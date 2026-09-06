"""Prepare, submit, and inspect immutable SLURM controller jobs."""

from __future__ import annotations

import json
import os
import re
import shlex
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from .hpc_config import HPCConfigError, load_hpc_settings, render_nextflow_config
from .plan import sha256_file


JOB_ID = re.compile(r"^(\d+)(?:;[A-Za-z0-9_.-]+)?$")
SAFE_ENVIRONMENT = {"NXF_OFFLINE", "NXF_SINGULARITY_CACHEDIR", "NXF_PLUGINS_DIR"}


class HPCRunError(ValueError):
    """Raised when an HPC launch or status operation cannot be trusted."""


def _read_plan(path: Path) -> dict[str, object]:
    try:
        payload = json.loads(path.resolve(strict=True).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise HPCRunError(f"could not read run plan: {exc}") from exc
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise HPCRunError("unsupported run plan schema")
    if payload.get("stage") != "planned_not_executed":
        raise HPCRunError("run plan is not in planned_not_executed state")
    execution = payload.get("execution")
    if not isinstance(execution, dict) or execution.get("executor") != "slurm":
        raise HPCRunError("run plan was not created for the slurm executor")
    argv = payload.get("command_argv")
    if not isinstance(argv, list) or not argv or not all(isinstance(v, str) for v in argv):
        raise HPCRunError("run plan command_argv is invalid")
    if argv[0:2] != ["nextflow", "run"] or "-resume" not in argv:
        raise HPCRunError("run plan is not a resumable Nextflow command")
    controls = payload.get("control_files")
    if not isinstance(controls, dict) or "infrastructure_config" not in controls:
        raise HPCRunError("run plan has no infrastructure config control")
    for label, record in controls.items():
        if not isinstance(record, dict):
            raise HPCRunError(f"invalid control file record: {label}")
        control_path = Path(str(record.get("path", ""))).resolve(strict=True)
        if control_path.stat().st_size != record.get("size_bytes"):
            raise HPCRunError(f"control file size changed after planning: {label}")
        if sha256_file(control_path) != record.get("sha256"):
            raise HPCRunError(f"control file checksum changed after planning: {label}")
    return payload


def _inside(root: Path, child: Path, label: str) -> None:
    try:
        child.resolve().relative_to(root.resolve())
    except ValueError as exc:
        raise HPCRunError(f"{label} must be inside configured work_root") from exc


def prepare_launcher(
    *,
    run_plan: Path,
    settings_path: Path,
    output: Path,
) -> dict[str, object]:
    """Validate all controls and exclusively render a controller sbatch script."""

    plan = _read_plan(run_plan)
    try:
        settings = load_hpc_settings(settings_path)
    except HPCConfigError as exc:
        raise HPCRunError(str(exc)) from exc
    controls = plan["control_files"]
    config_path = Path(str(controls["infrastructure_config"]["path"]))
    expected_config = render_nextflow_config(settings)
    if config_path.read_text(encoding="utf-8") != expected_config:
        raise HPCRunError("infrastructure config does not match the supplied settings")

    execution = plan["execution"]
    _inside(
        Path(settings.paths.work_root),
        Path(str(execution["workdir"])),
        "planned work directory",
    )
    required_environment = execution.get("required_environment", {})
    if not isinstance(required_environment, dict):
        raise HPCRunError("required_environment is invalid")
    if not set(required_environment).issubset(SAFE_ENVIRONMENT):
        raise HPCRunError("run plan requests an unsupported environment variable")
    if execution.get("network_mode") == "offline":
        if required_environment.get("NXF_OFFLINE") != "true":
            raise HPCRunError("offline run plan does not enforce NXF_OFFLINE=true")

    output_resolved = output.resolve()
    log_dir = output_resolved.parent / "logs"
    cluster = settings.cluster
    lines = [
        "#!/usr/bin/env bash",
        f"#SBATCH --job-name={cluster.job_name}",
        f"#SBATCH --account={cluster.account}",
        f"#SBATCH --partition={cluster.partition}",
        f"#SBATCH --cpus-per-task={cluster.launcher_cpus}",
        f"#SBATCH --mem={cluster.launcher_memory_gb}G",
        f"#SBATCH --time={cluster.launcher_time}",
        "#SBATCH --export=ALL",
        f"#SBATCH --output={log_dir / 'controller-%j.out'}",
        f"#SBATCH --error={log_dir / 'controller-%j.err'}",
        "",
        "set -euo pipefail",
        "umask 077",
        f"mkdir -p {shlex.quote(str(log_dir))}",
    ]
    if settings.software.purge_modules:
        lines.append("module purge")
    for module in settings.software.modules:
        lines.append(f"module load {shlex.quote(module)}")
    for name, value in sorted(required_environment.items()):
        if not isinstance(value, str) or "\x00" in value or "\n" in value:
            raise HPCRunError(f"unsafe environment value for {name}")
        lines.append(f"export {name}={shlex.quote(value)}")
    lines.extend(
        [
            f"cd {shlex.quote(str(output_resolved.parent))}",
            "echo '[controller] starting Nextflow' >&2",
            shlex.join(plan["command_argv"]),
            "echo '[controller] Nextflow completed' >&2",
            "",
        ]
    )
    rendered = "\n".join(lines)
    output_resolved.parent.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    try:
        descriptor = os.open(
            output_resolved, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o700
        )
    except FileExistsError as exc:
        raise HPCRunError(f"refusing to overwrite existing launcher: {output}") from exc
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(rendered)
    except Exception:
        output_resolved.unlink(missing_ok=True)
        raise
    return {
        "prepared": True,
        "launcher": str(output_resolved),
        "launcher_sha256": sha256_file(output_resolved),
        "run_plan_sha256": sha256_file(run_plan.resolve(strict=True)),
        "settings_sha256": sha256_file(settings_path.resolve(strict=True)),
        "resumable": True,
        "contains_secrets": False,
    }


def _write_receipt(path: Path, payload: dict[str, object], *, exclusive: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if exclusive:
        try:
            descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        except FileExistsError as exc:
            raise HPCRunError(f"refusing existing submission receipt: {path}") from exc
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        return
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.chmod(temporary, 0o600)
    temporary.replace(path)


def submit_launcher(*, launcher: Path, receipt: Path) -> dict[str, object]:
    """Reserve a receipt, submit once with sbatch, and record the scheduler ID."""

    launcher_resolved = launcher.resolve(strict=True)
    if not launcher_resolved.is_file():
        raise HPCRunError(f"launcher is not a file: {launcher}")
    payload: dict[str, object] = {
        "schema_version": 1,
        "status": "reserved",
        "reserved_at": datetime.now(timezone.utc).isoformat(),
        "launcher": str(launcher_resolved),
        "launcher_sha256": sha256_file(launcher_resolved),
        "contains_secrets": False,
    }
    receipt_resolved = receipt.resolve()
    _write_receipt(receipt_resolved, payload, exclusive=True)
    try:
        completed = subprocess.run(
            ["sbatch", "--parsable", str(launcher_resolved)],
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0:
            payload.update(
                status="submission_failed",
                returncode=completed.returncode,
                error=completed.stderr.strip(),
            )
            raise HPCRunError(f"sbatch failed with exit {completed.returncode}")
        match = JOB_ID.fullmatch(completed.stdout.strip())
        if match is None:
            payload.update(status="submission_failed", error="invalid sbatch job id")
            raise HPCRunError("sbatch returned an invalid job identifier")
        payload.update(
            status="submitted",
            submitted_at=datetime.now(timezone.utc).isoformat(),
            job_id=match.group(1),
        )
    except OSError as exc:
        payload.update(status="submission_failed", error=str(exc))
        raise HPCRunError(f"could not execute sbatch: {exc}") from exc
    finally:
        _write_receipt(receipt_resolved, payload, exclusive=False)
    return payload


def query_job(receipt: Path) -> dict[str, object]:
    """Query squeue first, then sacct for a completed or vanished job."""

    try:
        payload = json.loads(receipt.resolve(strict=True).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise HPCRunError(f"could not read submission receipt: {exc}") from exc
    if not isinstance(payload, dict) or payload.get("status") != "submitted":
        raise HPCRunError("receipt does not contain a submitted job")
    job_id = str(payload.get("job_id", ""))
    if JOB_ID.fullmatch(job_id) is None:
        raise HPCRunError("receipt contains an invalid job identifier")

    queued = subprocess.run(
        ["squeue", "--noheader", "--jobs", job_id, "--format", "%T|%M|%l|%R"],
        check=False,
        capture_output=True,
        text=True,
    )
    if queued.returncode == 0 and queued.stdout.strip():
        state, elapsed, limit, reason = queued.stdout.strip().split("|", 3)
        return {
            "job_id": job_id,
            "source": "squeue",
            "state": state,
            "elapsed": elapsed,
            "time_limit": limit,
            "location_or_reason": reason,
            "terminal": False,
        }

    accounting = subprocess.run(
        [
            "sacct",
            "--noheader",
            "--parsable2",
            "--jobs",
            job_id,
            "--format",
            "JobIDRaw,State,Elapsed,TotalCPU,MaxRSS,ExitCode",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if accounting.returncode != 0:
        raise HPCRunError("job is absent from squeue and sacct failed")
    rows = [line.split("|") for line in accounting.stdout.splitlines() if line.strip()]
    root_rows = [row for row in rows if len(row) >= 6 and row[0] == job_id]
    if not root_rows:
        raise HPCRunError("job is absent from both squeue and sacct")
    row = root_rows[0]
    state = row[1].split()[0].rstrip("+")
    terminal_states = {
        "BOOT_FAIL", "CANCELLED", "COMPLETED", "DEADLINE", "FAILED",
        "NODE_FAIL", "OUT_OF_MEMORY", "PREEMPTED", "TIMEOUT",
    }
    return {
        "job_id": job_id,
        "source": "sacct",
        "state": state,
        "elapsed": row[2],
        "total_cpu": row[3],
        "max_rss": row[4],
        "exit_code": row[5],
        "terminal": state in terminal_states,
        "successful": state == "COMPLETED" and row[5].startswith("0:0"),
    }
