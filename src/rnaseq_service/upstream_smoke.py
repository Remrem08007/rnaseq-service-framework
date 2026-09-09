"""Plan genuine pinned nf-core smoke runs using public upstream test data."""

from __future__ import annotations

import json
import os
import shlex
from datetime import datetime, timezone
from pathlib import Path

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

