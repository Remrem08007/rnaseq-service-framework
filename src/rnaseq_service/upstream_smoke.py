"""Plan genuine pinned nf-core smoke runs using public upstream test data."""

from __future__ import annotations

import csv
import json
import os
import shlex
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

from .hpc_config import MODULE_TOKEN, SLURM_TIME, TOKEN
from .plan import CONTAINER_ENGINES, PROXY_VARIABLES, sha256_file
from .resources import ResourceSummaryError, summarize_trace
from .workflow_lock import WorkflowLockError, load_workflow_lock


NETWORK_MODES = {"direct", "proxy"}
TEST_PROFILES = {
    "rnaseq": "test",
    "differential": "test_rnaseq_deseq2_gsea",
}
HashProgress = Callable[[str, int, int], None]


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
    if run_root.exists():
        raise UpstreamSmokeError(
            f"upstream smoke run root must not already exist: {run_root}"
        )
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
    try:
        lock_path = Path(str(lock.get("path", ""))).resolve(strict=True)
    except OSError as exc:
        raise UpstreamSmokeError("workflow lock is unavailable after smoke planning") from exc
    if lock_path.stat().st_size != lock.get("size_bytes"):
        raise UpstreamSmokeError("workflow lock size changed after smoke planning")
    if sha256_file(lock_path) != lock.get("sha256"):
        raise UpstreamSmokeError("workflow lock checksum changed after smoke planning")
    try:
        workflow_lock = load_workflow_lock(lock_path)
    except WorkflowLockError as exc:
        raise UpstreamSmokeError(str(exc)) from exc
    execution = plan.get("execution")
    if (
        not isinstance(execution, dict)
        or execution.get("network_mode") not in NETWORK_MODES
        or execution.get("container_engine") not in CONTAINER_ENGINES
        or not isinstance(execution.get("run_root"), str)
    ):
        raise UpstreamSmokeError("upstream smoke plan has an invalid execution record")
    if plan.get("runtime") != {
        "nextflow_version": workflow_lock.runtime.nextflow_version,
        "nf_core_tools_version": workflow_lock.runtime.nf_core_tools_version,
    }:
        raise UpstreamSmokeError("upstream smoke runtime pin changed after planning")
    stages = plan.get("stages")
    if (
        not isinstance(stages, list)
        or not all(isinstance(stage, dict) for stage in stages)
        or [stage.get("id") for stage in stages if isinstance(stage, dict)]
        != ["rnaseq", "differential"]
    ):
        raise UpstreamSmokeError("upstream smoke plan has invalid stages")
    expected_workflows = {
        "rnaseq": (workflow_lock.rnaseq.name, workflow_lock.rnaseq.revision),
        "differential": (
            workflow_lock.differential.name,
            workflow_lock.differential.revision,
        ),
    }
    run_root = Path(str(execution["run_root"])).resolve()
    for stage in stages:
        stage_id = str(stage.get("id"))
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
        workflow_name, revision = expected_workflows[stage_id]
        if stage.get("workflow") != {"name": workflow_name, "revision": revision}:
            raise UpstreamSmokeError(f"upstream smoke stage {stage_id!r} workflow pin changed")
        if stage.get("test_profile") != TEST_PROFILES[stage_id]:
            raise UpstreamSmokeError(f"upstream smoke stage {stage_id!r} test profile changed")
        expected_options = {
            "-log": run_root / stage_id / "execution" / "nextflow.log",
            "-r": revision,
            "-profile": (
                f"{TEST_PROFILES[stage_id]},{execution['container_engine']}"
            ),
            "-work-dir": run_root / stage_id / "work",
            "-with-report": run_root / stage_id / "execution" / "report.html",
            "-with-trace": run_root / stage_id / "execution" / "trace.tsv",
            "-with-timeline": run_root / stage_id / "execution" / "timeline.html",
            "-with-dag": run_root / stage_id / "execution" / "dag.html",
            "--outdir": run_root / stage_id / "results",
        }
        if argv[4] != workflow_name:
            raise UpstreamSmokeError(f"upstream smoke stage {stage_id!r} command workflow changed")
        for option, expected in expected_options.items():
            if argv.count(option) != 1:
                raise UpstreamSmokeError(
                    f"upstream smoke stage {stage_id!r} has invalid {option} option"
                )
            index = argv.index(option)
            if index + 1 >= len(argv) or str(argv[index + 1]) != str(expected):
                raise UpstreamSmokeError(
                    f"upstream smoke stage {stage_id!r} {option} value changed"
                )
    return resolved, plan


def _required_file(path: Path, role: str) -> Path:
    if not path.is_file():
        raise UpstreamSmokeError(f"missing required {role}: {path}")
    if path.stat().st_size == 0:
        raise UpstreamSmokeError(f"empty required {role}: {path}")
    return path.resolve()


def _matches(root: Path, pattern: str, role: str) -> list[Path]:
    matches = sorted(path.resolve() for path in root.glob(pattern) if path.is_file())
    if not matches:
        raise UpstreamSmokeError(f"missing required {role} matching {root / pattern}")
    for path in matches:
        _required_file(path, role)
        try:
            path.relative_to(root.resolve())
        except ValueError as exc:
            raise UpstreamSmokeError(f"{role} resolves outside its stage output: {path}") from exc
    return matches


def _one_match(root: Path, pattern: str, role: str) -> Path:
    matches = _matches(root, pattern, role)
    if len(matches) != 1:
        raise UpstreamSmokeError(
            f"expected exactly one {role} matching {root / pattern}; found {len(matches)}"
        )
    return matches[0]


def _read_params(path: Path, expected_outdir: Path) -> dict[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise UpstreamSmokeError(f"could not read pipeline params.json: {exc}") from exc
    if not isinstance(payload, dict):
        raise UpstreamSmokeError("pipeline params.json must contain a JSON object")
    value = payload.get("outdir")
    if not isinstance(value, str) or Path(value).resolve() != expected_outdir.resolve():
        raise UpstreamSmokeError("pipeline params.json 'outdir' conflicts with the run plan")
    return payload


def _hash_artifact(
    *,
    stage: str,
    role: str,
    path: Path,
    run_root: Path,
    progress: HashProgress | None,
) -> dict[str, object]:
    import hashlib

    try:
        relative_path = path.relative_to(run_root)
    except ValueError as exc:
        raise UpstreamSmokeError(
            f"{stage} artifact resolves outside the smoke run root: {path}"
        ) from exc
    total = path.stat().st_size
    completed = 0
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
            completed += len(chunk)
            if progress is not None:
                progress(f"{stage}:{role}", completed, total)
    if progress is not None and total == 0:
        progress(f"{stage}:{role}", 0, 0)
    return {
        "role": role,
        "path": str(path),
        "relative_path": str(relative_path),
        "size_bytes": total,
        "sha256": digest.hexdigest(),
    }


def _stage_artifacts(stage_id: str, stage_root: Path) -> list[tuple[str, Path]]:
    execution = stage_root / "execution"
    results = stage_root / "results"
    artifacts: list[tuple[str, Path]] = [
        (
            "nextflow_report",
            _required_file(execution / "report.html", "Nextflow report"),
        ),
        ("nextflow_trace", _required_file(execution / "trace.tsv", "Nextflow trace")),
        (
            "nextflow_timeline",
            _required_file(execution / "timeline.html", "Nextflow timeline"),
        ),
        ("nextflow_dag", _required_file(execution / "dag.html", "Nextflow DAG")),
        (
            "pipeline_params",
            _required_file(results / "pipeline_info" / "params.json", "pipeline parameters"),
        ),
        (
            "software_versions",
            _one_match(
                results / "pipeline_info",
                "*software*versions.yml",
                "software versions file",
            ),
        ),
    ]
    if stage_id == "rnaseq":
        artifacts.extend(
            [
                (
                    "multiqc_report",
                    _one_match(results, "multiqc/*/multiqc_report.html", "MultiQC report"),
                ),
                (
                    "gene_counts",
                    _matches(
                        results,
                        "*/salmon.merged.gene_counts.tsv",
                        "merged Salmon gene-count matrix",
                    )[0],
                ),
                (
                    "gene_tpm",
                    _matches(
                        results,
                        "*/salmon.merged.gene_tpm.tsv",
                        "merged Salmon gene-TPM matrix",
                    )[0],
                ),
            ]
        )
    else:
        artifacts.extend(
            [
                (
                    "analysis_report",
                    _one_match(
                        results,
                        "report/rnaseq_deseq2_gsea/*_differentialabundance_report.html",
                        "differential-abundance report",
                    ),
                ),
                (
                    "normalised_counts",
                    _one_match(
                        results,
                        "tables/processed_abundance/rnaseq_deseq2_gsea/*.normalised_counts.tsv",
                        "normalised-count matrix",
                    ),
                ),
                (
                    "vst_counts",
                    _one_match(
                        results,
                        "tables/processed_abundance/rnaseq_deseq2_gsea/*.vst.tsv",
                        "variance-stabilised matrix",
                    ),
                ),
                (
                    "deseq2_results",
                    _matches(
                        results,
                        "tables/differential/rnaseq_deseq2_gsea/*.deseq2.results.tsv",
                        "DESeq2 result table",
                    )[0],
                ),
                (
                    "deseq2_filtered",
                    _matches(
                        results,
                        "tables/differential/rnaseq_deseq2_gsea/*.deseq2.results_filtered.tsv",
                        "filtered DESeq2 result table",
                    )[0],
                ),
                (
                    "volcano_plot",
                    _matches(
                        results,
                        "plots/differential/rnaseq_deseq2_gsea/*/png/volcano.png",
                        "volcano plot",
                    )[0],
                ),
                (
                    "gsea_report",
                    _matches(results, "report/gsea/**/*.Gsea.rpt", "GSEA report")[0],
                ),
            ]
        )
    return artifacts


def inspect_upstream_smoke_completion(
    run_plan: Path,
    *,
    progress: HashProgress | None = None,
) -> dict[str, object]:
    """Validate both genuine upstream test-profile runs without writing a receipt."""

    plan_path, plan = _read_plan(run_plan)
    run_root = Path(str(plan["execution"]["run_root"])).resolve()
    runtime_path = _required_file(
        run_root / "runtime" / "versions.tsv", "runtime-version evidence"
    )
    try:
        with runtime_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["tool", "version"]:
                raise UpstreamSmokeError("runtime-version evidence has invalid columns")
            runtime_rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as exc:
        raise UpstreamSmokeError(f"could not read runtime-version evidence: {exc}") from exc
    observed = {row["tool"]: row["version"] for row in runtime_rows}
    expected_nextflow = str(plan["runtime"]["nextflow_version"])
    expected_engine = str(plan["execution"]["container_engine"])
    if observed.get("nextflow") != expected_nextflow:
        raise UpstreamSmokeError(
            "observed Nextflow version does not match the workflow lock"
        )
    if observed.get("container_engine") != expected_engine:
        raise UpstreamSmokeError("observed container engine does not match the run plan")
    if not observed.get("container_runtime"):
        raise UpstreamSmokeError("runtime-version evidence has no container runtime version")
    runtime_artifact = _hash_artifact(
        stage="runtime",
        role="versions",
        path=runtime_path,
        run_root=run_root,
        progress=progress,
    )
    stage_records: list[dict[str, object]] = []
    for stage in plan["stages"]:
        stage_id = str(stage["id"])
        stage_root = run_root / stage_id
        artifacts = _stage_artifacts(stage_id, stage_root)
        by_role = dict(artifacts)
        try:
            trace_summary = summarize_trace(by_role["nextflow_trace"])
        except ResourceSummaryError as exc:
            raise UpstreamSmokeError(str(exc)) from exc
        if not trace_summary["all_tasks_successful"]:
            raise UpstreamSmokeError(
                f"{stage_id} trace is not complete and successful: "
                f"failed={trace_summary['failed_task_count']}, "
                f"nonterminal_or_unknown={trace_summary['other_task_count']}"
            )
        _read_params(by_role["pipeline_params"], stage_root / "results")
        inventory = [
            _hash_artifact(
                stage=stage_id,
                role=role,
                path=path,
                run_root=run_root,
                progress=progress,
            )
            for role, path in artifacts
        ]
        stage_records.append(
            {
                "id": stage_id,
                "workflow": stage["workflow"],
                "test_profile": stage["test_profile"],
                "task_summary": {
                    "task_count": trace_summary["task_count"],
                    "process_count": trace_summary["process_count"],
                    "failed_task_count": trace_summary["failed_task_count"],
                    "other_task_count": trace_summary["other_task_count"],
                },
                "artifacts": inventory,
            }
        )
    return {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "stage": "upstream_smoke_complete",
        "run_plan": _control(plan_path),
        "workflow_lock": plan["workflow_lock"],
        "runtime": plan["runtime"],
        "runtime_verification": {
            "expected_nextflow_version": expected_nextflow,
            "observed_nextflow_version": observed["nextflow"],
            "container_engine": observed["container_engine"],
            "container_runtime": observed["container_runtime"],
            "artifact": runtime_artifact,
        },
        "execution": plan["execution"],
        "stages": stage_records,
        "controls_verified": True,
        "scientific_execution_performed": True,
        "uses_public_upstream_test_data": True,
        "contains_client_data": False,
        "contains_secrets": False,
        "validation_scope": "official_upstream_test_profiles_not_client_study_validation",
    }


def create_upstream_smoke_receipt(
    *,
    run_plan: Path,
    output: Path,
    progress: HashProgress | None = None,
) -> dict[str, object]:
    """Exclusively seal a successful two-stage upstream smoke run."""

    if output.exists():
        raise UpstreamSmokeError(f"refusing to overwrite completion receipt: {output}")
    payload = inspect_upstream_smoke_completion(run_plan, progress=progress)
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        descriptor = os.open(output, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    except FileExistsError as exc:
        raise UpstreamSmokeError(f"refusing to overwrite completion receipt: {output}") from exc
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
    except Exception:
        output.unlink(missing_ok=True)
        raise
    return payload


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
        f"mkdir -p {shlex.quote(str(run_root / 'runtime'))}",
    ]
    if purge_modules:
        lines.append("module purge")
    lines.extend(f"module load {shlex.quote(module)}" for module in modules)
    lines.extend(
        [
            "",
            f"expected_nextflow={shlex.quote(str(plan['runtime']['nextflow_version']))}",
            f"container_engine={shlex.quote(str(plan['execution']['container_engine']))}",
            "actual_nextflow=$(nextflow -version 2>&1 | awk '/version/ && !found {for (i=1; i<=NF; i++) if ($i ~ /^[0-9]+\\.[0-9]+\\.[0-9]+(-edge)?$/) {print $i; found=1; break}}')",
            "if [[ \"$actual_nextflow\" != \"$expected_nextflow\" ]]; then",
            "    echo \"[upstream-smoke] Nextflow version mismatch: expected=${expected_nextflow} observed=${actual_nextflow:-unavailable}\" >&2",
            "    exit 64",
            "fi",
            "if ! command -v \"$container_engine\" >/dev/null 2>&1; then",
            "    echo \"[upstream-smoke] container engine unavailable: ${container_engine}\" >&2",
            "    exit 64",
            "fi",
            "container_runtime=$(\"$container_engine\" --version 2>&1 | sed -n '1p' | tr '\\t\\r\\n' '   ')",
            "if [[ -z \"$container_runtime\" ]]; then",
            "    echo \"[upstream-smoke] container runtime version is unavailable\" >&2",
            "    exit 64",
            "fi",
            f"runtime_evidence={shlex.quote(str(run_root / 'runtime' / 'versions.tsv'))}",
            "printf 'tool\\tversion\\nnextflow\\t%s\\ncontainer_engine\\t%s\\ncontainer_runtime\\t%s\\n' \\",
            "    \"$actual_nextflow\" \"$container_engine\" \"$container_runtime\" > \"$runtime_evidence\"",
            "chmod 600 \"$runtime_evidence\"",
            "echo \"[upstream-smoke] runtime verified nextflow=${actual_nextflow} engine=${container_engine}\" >&2",
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
