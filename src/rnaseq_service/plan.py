"""Create immutable, inspectable RNA-seq primary-processing run plans."""

from __future__ import annotations

import hashlib
import json
import os
import shlex
import re
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

from .bundle import BundleError, verify_bundle
from .inputs import InputManifestError, verify_input_manifest
from .preflight import PreflightReport, run_preflight
from .workflow_lock import WorkflowLockError, load_workflow_lock


NETWORK_MODES = {"direct", "proxy", "offline"}
CONTAINER_ENGINES = {"apptainer", "docker", "singularity"}
EXECUTORS = {"local", "slurm"}
PROXY_VARIABLES = (
    "HTTP_PROXY",
    "HTTPS_PROXY",
    "NO_PROXY",
    "http_proxy",
    "https_proxy",
    "no_proxy",
)
REFERENCE_KEY = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


class PlanError(ValueError):
    """Raised when a safe run plan cannot be created."""


class PreflightPlanError(PlanError):
    def __init__(self, report: PreflightReport):
        super().__init__("intake preflight failed")
        self.report = report


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _control_file(path: Path) -> dict[str, object]:
    resolved = path.resolve(strict=True)
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256_file(resolved),
    }


def build_rnaseq_argv(
    *,
    workflow_name: str,
    revision: str | None,
    samplesheet: Path,
    outdir: Path,
    workdir: Path,
    container_engine: str,
    infrastructure_config: Path | None,
    genome: str | None,
) -> list[str]:
    """Return an argv vector; callers do not need to construct a shell string."""

    argv = [
        "nextflow",
        "run",
        workflow_name,
    ]
    if revision is not None:
        argv.extend(["-r", revision])
    argv.extend(
        ["-profile", container_engine, "-work-dir", str(workdir)]
    )
    if infrastructure_config is not None:
        argv.extend(["-c", str(infrastructure_config)])
    argv.extend(
        [
            "-resume",
            "-with-report",
            str(outdir / "execution" / "report.html"),
            "-with-trace",
            str(outdir / "execution" / "trace.tsv"),
            "-with-timeline",
            str(outdir / "execution" / "timeline.html"),
            "-with-dag",
            str(outdir / "execution" / "dag.html"),
            "--input",
            str(samplesheet),
            "--outdir",
            str(outdir / "rnaseq"),
        ]
    )
    if genome is not None:
        argv.extend(["--genome", genome])
    return argv


def create_rnaseq_plan(
    *,
    samplesheet: Path,
    design: Path,
    contrasts: Path,
    workflow_lock: Path,
    output: Path,
    outdir: Path,
    workdir: Path,
    network_mode: str,
    container_engine: str,
    executor: str,
    infrastructure_config: Path | None = None,
    offline_manifest: Path | None = None,
    input_manifest: Path | None = None,
    genome: str | None = None,
    min_replicates: int = 2,
) -> dict[str, object]:
    """Validate inputs and exclusively write a reproducible plan JSON."""

    if network_mode not in NETWORK_MODES:
        raise PlanError(f"unsupported network mode: {network_mode!r}")
    if container_engine not in CONTAINER_ENGINES:
        raise PlanError(f"unsupported container engine: {container_engine!r}")
    if executor not in EXECUTORS:
        raise PlanError(f"unsupported executor: {executor!r}")
    if executor == "slurm" and infrastructure_config is None:
        raise PlanError("--infrastructure-config is required for the slurm executor")
    if executor == "local" and infrastructure_config is not None:
        raise PlanError("--infrastructure-config is only accepted with the slurm executor")
    if network_mode == "offline" and offline_manifest is None:
        raise PlanError("--offline-manifest is required in offline mode")
    if network_mode != "offline" and offline_manifest is not None:
        raise PlanError("--offline-manifest is only accepted in offline mode")
    if network_mode == "offline" and container_engine == "docker":
        raise PlanError("offline bundles currently support apptainer/singularity only")
    if genome is not None and REFERENCE_KEY.fullmatch(genome) is None:
        raise PlanError(f"invalid iGenomes reference key: {genome!r}")

    preflight = run_preflight(
        samplesheet,
        design,
        contrasts,
        check_files=True,
        min_replicates=min_replicates,
    )
    if not preflight.valid:
        raise PreflightPlanError(preflight)

    try:
        lock = load_workflow_lock(workflow_lock)
    except WorkflowLockError as exc:
        raise PlanError(str(exc)) from exc

    bundle: dict[str, object] | None = None
    offline_manifest_resolved: Path | None = None
    if offline_manifest is not None:
        try:
            bundle = verify_bundle(
                offline_manifest,
                expected_workflow_lock=workflow_lock,
            )
        except (BundleError, OSError) as exc:
            raise PlanError(f"offline bundle verification failed: {exc}") from exc
        offline_manifest_resolved = offline_manifest.resolve(strict=True)

    samplesheet_resolved = samplesheet.resolve(strict=True)
    input_dataset: dict[str, object] | None = None
    input_manifest_resolved: Path | None = None
    if input_manifest is not None:
        try:
            input_dataset = verify_input_manifest(
                input_manifest,
                expected_samplesheet=samplesheet_resolved,
            )
        except (InputManifestError, OSError) as exc:
            raise PlanError(f"input manifest verification failed: {exc}") from exc
        input_manifest_resolved = input_manifest.resolve(strict=True)
    outdir_resolved = outdir.resolve()
    workdir_resolved = workdir.resolve()
    infra_resolved = (
        infrastructure_config.resolve(strict=True)
        if infrastructure_config is not None
        else None
    )
    workflow_source = lock.rnaseq.name
    workflow_revision: str | None = lock.rnaseq.revision
    if bundle is not None and offline_manifest_resolved is not None:
        components = bundle["components"]
        workflow_source = str(
            (offline_manifest_resolved.parent / components["rnaseq_workflow"]).resolve()
        )
        workflow_revision = None
    argv = build_rnaseq_argv(
        workflow_name=workflow_source,
        revision=workflow_revision,
        samplesheet=samplesheet_resolved,
        outdir=outdir_resolved,
        workdir=workdir_resolved,
        container_engine=container_engine,
        infrastructure_config=infra_resolved,
        genome=genome,
    )

    controls = {
        "samplesheet": _control_file(samplesheet),
        "design": _control_file(design),
        "contrasts": _control_file(contrasts),
        "workflow_lock": _control_file(workflow_lock),
    }
    if infra_resolved is not None:
        controls["infrastructure_config"] = _control_file(infra_resolved)
    if offline_manifest_resolved is not None:
        controls["offline_manifest"] = _control_file(offline_manifest_resolved)
    if input_manifest_resolved is not None:
        controls["input_manifest"] = _control_file(input_manifest_resolved)

    required_environment: dict[str, str] = {}
    if bundle is not None and offline_manifest_resolved is not None:
        components = bundle["components"]
        required_environment = {
            "NXF_OFFLINE": "true",
            "NXF_SINGULARITY_CACHEDIR": str(
                (
                    offline_manifest_resolved.parent / components["container_root"]
                ).resolve()
            ),
            "NXF_PLUGINS_DIR": str(
                (
                    offline_manifest_resolved.parent / components["plugin_root"]
                ).resolve()
            ),
        }

    plan: dict[str, object] = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "stage": "planned_not_executed",
        "workflow": asdict(lock.rnaseq),
        "runtime": asdict(lock.runtime),
        "workflow_lock_last_reviewed": lock.last_reviewed,
        "control_files": controls,
        "intake_policy": {"min_replicates": min_replicates},
        "preflight": asdict(preflight),
        "execution": {
            "executor": executor,
            "container_engine": container_engine,
            "network_mode": network_mode,
            "outdir": str(outdir_resolved),
            "workdir": str(workdir_resolved),
            "proxy_environment_present": [
                name for name in PROXY_VARIABLES if os.environ.get(name)
            ],
            "proxy_values_recorded": False,
            "offline_bundle_verified": bundle is not None,
            "offline_bundle": bundle,
            "required_environment": required_environment,
        },
        "command_argv": argv,
        "command_preview": shlex.join(argv),
        "contains_client_results": False,
        "contains_secrets": False,
    }
    if genome is not None:
        plan["reference"] = {"mode": "igenomes", "genome": genome}
    if input_dataset is not None:
        plan["input_dataset"] = {
            "manifest": str(input_manifest_resolved),
            "n_fastq_files": input_dataset["n_fastq_files"],
            "total_bytes": input_dataset["total_bytes"],
            "metadata_verified_at_planning": True,
        }

    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        descriptor = os.open(output, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    except FileExistsError as exc:
        raise PlanError(f"refusing to overwrite existing plan: {output}") from exc
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(plan, handle, indent=2, sort_keys=True)
            handle.write("\n")
    except Exception:
        output.unlink(missing_ok=True)
        raise
    return plan
