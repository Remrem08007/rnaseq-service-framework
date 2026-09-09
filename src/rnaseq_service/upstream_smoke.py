"""Plan genuine pinned nf-core smoke runs using public upstream test data."""

from __future__ import annotations

import json
import os
import shlex
from datetime import datetime, timezone
from pathlib import Path

from .hpc_config import MODULE_TOKEN, SLURM_TIME, TOKEN
from .plan import CONTAINER_ENGINES, PROXY_VARIABLES, sha256_file
from .workflow_lock import WorkflowLockError, load_workflow_lock


NETWORK_MODES = {"direct", "proxy"}
TEST_PROFILES = {
    "rnaseq": "test",
    "differential": "test_rnaseq_deseq2_gsea",
}


class UpstreamSmokeError(ValueError):
    """Raised when a real upstream smoke plan cannot be trusted."""


def _control(path: Path) -> dict[str, object]:
    resolved = path.resolve(strict=True)
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256_file(resolved),
    }


def _stage_argv(
    *,
    workflow: str,
    revision: str,
    profile: str,
    container_engine: str,
    root: Path,
    stage: str,
) -> list[str]:
    execution = root / stage / "execution"
    return [
        "nextflow",
        "-log",
        str(execution / "nextflow.log"),
        "run",
        workflow,
        "-r",
        revision,
        "-profile",
        f"{profile},{container_engine}",
        "-work-dir",
        str(root / stage / "work"),
        "-resume",
        "-with-report",
        str(execution / "report.html"),
        "-with-trace",
        str(execution / "trace.tsv"),
        "-with-timeline",
        str(execution / "timeline.html"),
        "-with-dag",
        str(execution / "dag.html"),
        "--outdir",
        str(root / stage / "results"),
    ]


def create_upstream_smoke_plan(
    *,
    workflow_lock: Path,
    output: Path,
    run_root: Path,
    network_mode: str,
    container_engine: str = "apptainer",
) -> dict[str, object]:
    """Write a non-executing plan for official public-data test profiles."""

    if output.exists():
        raise UpstreamSmokeError(f"refusing to overwrite upstream smoke plan: {output}")
    if network_mode not in NETWORK_MODES:
        raise UpstreamSmokeError("upstream smoke network mode must be 'direct' or 'proxy'")
    if container_engine not in CONTAINER_ENGINES:
        raise UpstreamSmokeError(f"unsupported container engine: {container_engine!r}")
    try:
        lock = load_workflow_lock(workflow_lock)
    except WorkflowLockError as exc:
        raise UpstreamSmokeError(str(exc)) from exc
    root = run_root.resolve()
    stages = [
        {
            "id": "rnaseq",
            "workflow": {"name": lock.rnaseq.name, "revision": lock.rnaseq.revision},
            "test_profile": TEST_PROFILES["rnaseq"],
            "command_argv": _stage_argv(
                workflow=lock.rnaseq.name,
                revision=lock.rnaseq.revision,
                profile=TEST_PROFILES["rnaseq"],
                container_engine=container_engine,
                root=root,
                stage="rnaseq",
            ),
        },
        {
            "id": "differential",
            "workflow": {
                "name": lock.differential.name,
                "revision": lock.differential.revision,
            },
            "test_profile": TEST_PROFILES["differential"],
            "command_argv": _stage_argv(
                workflow=lock.differential.name,
                revision=lock.differential.revision,
                profile=TEST_PROFILES["differential"],
                container_engine=container_engine,
                root=root,
                stage="differential",
            ),
        },
    ]
    for stage in stages:
        stage["command_preview"] = shlex.join(stage["command_argv"])
    plan: dict[str, object] = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "stage": "upstream_smoke_planned_not_executed",
        "purpose": "real_pinned_upstream_test_profiles",
        "workflow_lock": _control(workflow_lock),
        "runtime": {
            "nextflow_version": lock.runtime.nextflow_version,
            "nf_core_tools_version": lock.runtime.nf_core_tools_version,
        },
        "execution": {
            "run_root": str(root),
            "network_mode": network_mode,
            "container_engine": container_engine,
            "proxy_environment_present": [
                name for name in PROXY_VARIABLES if os.environ.get(name)
            ],
            "proxy_values_recorded": False,
            "resume_required": True,
        },
        "stages": stages,
        "uses_public_upstream_test_data": True,
        "scientific_execution_expected": True,
        "contains_client_data": False,
        "contains_secrets": False,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        descriptor = os.open(output, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    except FileExistsError as exc:
        raise UpstreamSmokeError(f"refusing to overwrite upstream smoke plan: {output}") from exc
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(plan, handle, indent=2, sort_keys=True)
            handle.write("\n")
    except Exception:
        output.unlink(missing_ok=True)
        raise
    return plan


def _read_plan(path: Path) -> tuple[Path, dict[str, object]]:
    try:
        resolved = path.resolve(strict=True)
        plan = json.loads(resolved.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise UpstreamSmokeError(f"could not read upstream smoke plan: {exc}") from exc
    if (
        not isinstance(plan, dict)
        or plan.get("schema_version") != 1
        or plan.get("stage") != "upstream_smoke_planned_not_executed"
        or plan.get("uses_public_upstream_test_data") is not True
        or plan.get("contains_client_data") is not False
        or plan.get("contains_secrets") is not False
    ):
        raise UpstreamSmokeError("upstream smoke plan is not supported")
    lock = plan.get("workflow_lock")
    if not isinstance(lock, dict):
        raise UpstreamSmokeError("upstream smoke plan has no workflow-lock control")
    lock_path = Path(str(lock.get("path", ""))).resolve(strict=True)
    if lock_path.stat().st_size != lock.get("size_bytes"):
        raise UpstreamSmokeError("workflow lock size changed after smoke planning")
    if sha256_file(lock_path) != lock.get("sha256"):
        raise UpstreamSmokeError("workflow lock checksum changed after smoke planning")
    stages = plan.get("stages")
    if (
        not isinstance(stages, list)
        or [stage.get("id") for stage in stages if isinstance(stage, dict)]
        != ["rnaseq", "differential"]
    ):
        raise UpstreamSmokeError("upstream smoke plan has invalid stages")
    for stage in stages:
        argv = stage.get("command_argv")
        if (
            not isinstance(argv, list)
            or not argv
            or not all(isinstance(value, str) for value in argv)
            or argv[:4:3] != ["nextflow", "run"]
            or "-resume" not in argv
        ):
            raise UpstreamSmokeError(f"upstream smoke stage {stage.get('id')!r} is invalid")
        if stage.get("command_preview") != shlex.join(argv):
            raise UpstreamSmokeError(f"upstream smoke stage {stage.get('id')!r} preview changed")
    return resolved, plan


def create_upstream_smoke_launcher(
    *,
    run_plan: Path,
    output: Path,
    account: str,
    partition: str,
    wall_time: str,
    cpus: int,
    memory_gb: int,
    modules: list[str],
    purge_modules: bool = True,
) -> dict[str, object]:
    """Render, but never submit, a two-stage SLURM smoke launcher."""

    plan_path, plan = _read_plan(run_plan)
    if output.exists():
        raise UpstreamSmokeError(f"refusing to overwrite upstream smoke launcher: {output}")
    for label, value in (("account", account), ("partition", partition)):
        if TOKEN.fullmatch(value) is None:
            raise UpstreamSmokeError(f"invalid SLURM {label}: {value!r}")
    if SLURM_TIME.fullmatch(wall_time) is None:
        raise UpstreamSmokeError(f"invalid SLURM wall time: {wall_time!r}")
    if isinstance(cpus, bool) or not isinstance(cpus, int) or cpus < 1:
        raise UpstreamSmokeError("SLURM CPUs must be a positive integer")
    if isinstance(memory_gb, bool) or not isinstance(memory_gb, int) or memory_gb < 1:
        raise UpstreamSmokeError("SLURM memory must be a positive integer")
    if not modules or any(MODULE_TOKEN.fullmatch(module) is None for module in modules):
        raise UpstreamSmokeError("module names must be non-empty safe tokens")

    resolved = output.resolve()
    log_dir = resolved.parent / "logs"
    run_root = Path(str(plan["execution"]["run_root"]))
    lines = [
        "#!/usr/bin/env bash",
        "#SBATCH --job-name=rnaseq-upstream-smoke",
        f"#SBATCH --account={account}",
        f"#SBATCH --partition={partition}",
        f"#SBATCH --cpus-per-task={cpus}",
        f"#SBATCH --mem={memory_gb}G",
        f"#SBATCH --time={wall_time}",
        "#SBATCH --export=ALL",
        f"#SBATCH --output={log_dir / 'smoke-%j.out'}",
        f"#SBATCH --error={log_dir / 'smoke-%j.err'}",
        "",
        "set -euo pipefail",
        "umask 077",
        f"mkdir -p {shlex.quote(str(log_dir))}",
        f"mkdir -p {shlex.quote(str(run_root / 'rnaseq' / 'execution'))}",
        f"mkdir -p {shlex.quote(str(run_root / 'differential' / 'execution'))}",
    ]
    if purge_modules:
        lines.append("module purge")
    lines.extend(f"module load {shlex.quote(module)}" for module in modules)
    lines.extend(
        [
            "",
            "run_stage() {",
            "    local label=\"$1\"",
            "    local ordinal=\"$2\"",
            "    shift 2",
            "    local started child ticker status",
            "    started=$(date +%s)",
            "    echo \"[upstream-smoke] stage=${ordinal}/2 name=${label} starting\" >&2",
            "    \"$@\" &",
            "    child=$!",
            "    (",
            "        while kill -0 \"$child\" 2>/dev/null; do",
            "            sleep 60",
            "            if kill -0 \"$child\" 2>/dev/null; then",
            "                echo \"[upstream-smoke] stage=${ordinal}/2 name=${label} running elapsed_seconds=$(($(date +%s)-started))\" >&2",
            "            fi",
            "        done",
            "    ) &",
            "    ticker=$!",
            "    set +e",
            "    wait \"$child\"",
            "    status=$?",
            "    set -e",
            "    kill \"$ticker\" 2>/dev/null || true",
            "    wait \"$ticker\" 2>/dev/null || true",
            "    if (( status != 0 )); then",
            "        echo \"[upstream-smoke] stage=${ordinal}/2 name=${label} failed exit=${status}\" >&2",
            "        return \"$status\"",
            "    fi",
            "    echo \"[upstream-smoke] stage=${ordinal}/2 name=${label} complete elapsed_seconds=$(($(date +%s)-started))\" >&2",
            "}",
            "",
        ]
    )
    for ordinal, stage in enumerate(plan["stages"], start=1):
        lines.append(
            f"run_stage {shlex.quote(str(stage['id']))} {ordinal} "
            + shlex.join(stage["command_argv"])
        )
    lines.extend(["", "echo '[upstream-smoke] COMPLETE' >&2", ""])
    rendered = "\n".join(lines)
    resolved.parent.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    try:
        descriptor = os.open(resolved, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o700)
    except FileExistsError as exc:
        raise UpstreamSmokeError(
            f"refusing to overwrite upstream smoke launcher: {output}"
        ) from exc
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(rendered)
    except Exception:
        resolved.unlink(missing_ok=True)
        raise
    return {
        "schema_version": 1,
        "prepared": True,
        "submitted": False,
        "launcher": str(resolved),
        "launcher_sha256": sha256_file(resolved),
        "run_plan": str(plan_path),
        "run_plan_sha256": sha256_file(plan_path),
        "heartbeat_seconds": 60,
        "uses_public_upstream_test_data": True,
        "contains_client_data": False,
        "contains_secrets": False,
    }
