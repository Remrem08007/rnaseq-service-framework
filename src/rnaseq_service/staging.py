"""Plan and execute pinned nf-core workflow/container staging."""

from __future__ import annotations

import json
import os
import queue
import re
import shutil
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.request
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

from .plan import PROXY_VARIABLES, sha256_file
from .workflow_lock import WorkflowLockError, WorkflowSpec, load_workflow_lock


STAGING_NETWORK_MODES = {"direct", "proxy"}
STAGING_ENDPOINTS = (
    ("github_api", "https://api.github.com/"),
    ("quay_registry", "https://quay.io/v2/"),
    ("docker_registry", "https://registry-1.docker.io/v2/"),
)


class StagingError(ValueError):
    """Raised when staging cannot proceed safely."""


def build_download_argv(
    spec: WorkflowSpec,
    *,
    outdir: Path,
    parallel_downloads: int,
) -> list[str]:
    if parallel_downloads < 1:
        raise StagingError("parallel downloads must be at least 1")
    short_name = spec.name.removeprefix("nf-core/")
    return [
        "nf-core",
        "pipelines",
        "download",
        short_name,
        "--revision",
        spec.revision,
        "--outdir",
        str(outdir),
        "--container-system",
        "singularity",
        "--container-cache-utilisation",
        "amend",
        "--compress",
        "none",
        "--parallel-downloads",
        str(parallel_downloads),
        "--download-configuration",
        "yes",
    ]


def _exclusive_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    except FileExistsError as exc:
        raise StagingError(f"refusing to overwrite existing file: {path}") from exc
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
    except Exception:
        path.unlink(missing_ok=True)
        raise


def create_staging_plan(
    *,
    workflow_lock: Path,
    bundle_dir: Path,
    output: Path,
    network_mode: str,
    parallel_downloads: int = 4,
) -> dict[str, object]:
    """Write an inspectable staging plan without contacting the network."""

    if network_mode not in STAGING_NETWORK_MODES:
        raise StagingError("staging network mode must be direct or proxy")
    if parallel_downloads < 1:
        raise StagingError("parallel downloads must be at least 1")
    present = [name for name in PROXY_VARIABLES if os.environ.get(name)]
    if network_mode == "proxy" and not present:
        raise StagingError("proxy mode requires a standard proxy environment variable")
    try:
        lock = load_workflow_lock(workflow_lock)
    except WorkflowLockError as exc:
        raise StagingError(str(exc)) from exc

    root = bundle_dir.resolve()
    pipelines = root / "pipelines"
    containers = root / "containers"
    plugins = root / "plugins"
    commands = []
    for spec in (lock.rnaseq, lock.differential):
        short_name = spec.name.removeprefix("nf-core/")
        target = pipelines / f"{short_name}-{spec.revision}"
        commands.append(
            {
                "workflow": asdict(spec),
                "output_root": str(target),
                "argv": build_download_argv(
                    spec,
                    outdir=target,
                    parallel_downloads=parallel_downloads,
                ),
            }
        )
    lock_resolved = workflow_lock.resolve(strict=True)
    plan: dict[str, object] = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "stage": "planned_not_executed",
        "network_mode": network_mode,
        "bundle_dir": str(root),
        "container_root": str(containers),
        "plugin_root": str(plugins),
        "parallel_downloads": parallel_downloads,
        "workflow_lock": {
            "path": str(lock_resolved),
            "sha256": sha256_file(lock_resolved),
        },
        "runtime": asdict(lock.runtime),
        "commands": commands,
        "runtime_environment": {
            "NXF_SINGULARITY_CACHEDIR": str(containers),
            "NXF_PLUGINS_DIR": str(plugins),
            "proxy_environment_present": present,
            "proxy_values_recorded": False,
        },
        "contains_client_data": False,
        "contains_secrets": False,
    }
    _exclusive_json(output, plan)
    return plan


def _redact(text: str, secrets: list[str]) -> str:
    rendered = text
    for value in sorted(set(secrets), key=len, reverse=True):
        if value:
            rendered = rendered.replace(value, "<redacted-proxy>")
    return rendered


def check_staging_environment(
    plan_path: Path,
    *,
    timeout_seconds: int = 10,
) -> dict[str, object]:
    """Check tools, exact versions, proxy presence, and required endpoints."""

    if timeout_seconds < 1:
        raise StagingError("network timeout must be at least 1 second")
    resolved_plan = plan_path.resolve(strict=True)
    try:
        plan = json.loads(resolved_plan.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise StagingError(f"could not read staging plan: {exc}") from exc
    if not isinstance(plan, dict) or plan.get("schema_version") != 1:
        raise StagingError("unsupported staging plan schema")
    network_mode = plan.get("network_mode")
    if network_mode not in STAGING_NETWORK_MODES:
        raise StagingError("invalid staging network mode")
    lock_record = plan.get("workflow_lock")
    if not isinstance(lock_record, dict):
        raise StagingError("staging plan has no workflow lock")
    lock_path = Path(str(lock_record.get("path", ""))).resolve(strict=True)
    if sha256_file(lock_path) != lock_record.get("sha256"):
        raise StagingError("workflow lock changed after planning")
    lock = load_workflow_lock(lock_path)

    proxy_names = [name for name in PROXY_VARIABLES if os.environ.get(name)]
    secrets = [os.environ[name] for name in proxy_names]
    checks: list[dict[str, object]] = []
    if network_mode == "proxy" and not proxy_names:
        checks.append(
            {
                "name": "proxy_environment",
                "ok": False,
                "detail": "no standard proxy variable is present",
            }
        )
    else:
        checks.append(
            {
                "name": "proxy_environment",
                "ok": True,
                "present_names": proxy_names,
                "values_recorded": False,
            }
        )

    environment = os.environ.copy()
    expected_versions = {
        "nf-core": lock.runtime.nf_core_tools_version,
        "nextflow": lock.runtime.nextflow_version,
    }
    version_argv = {"nf-core": ["nf-core", "--version"], "nextflow": ["nextflow", "-version"]}
    for executable, expected in expected_versions.items():
        available = shutil.which(executable, path=environment.get("PATH")) is not None
        detail: dict[str, object] = {
            "name": f"runtime:{executable}",
            "ok": False,
            "expected_version": expected,
        }
        if available:
            try:
                output = subprocess.check_output(
                    version_argv[executable],
                    stderr=subprocess.STDOUT,
                    text=True,
                    env=environment,
                )
                detail["ok"] = re.search(
                    rf"(?<![0-9.]){re.escape(expected)}(?![0-9.])", output
                ) is not None
                if not detail["ok"]:
                    detail["detail"] = _redact(output.strip(), secrets)
            except (OSError, subprocess.CalledProcessError) as exc:
                detail["detail"] = _redact(str(exc), secrets)
        else:
            detail["detail"] = "not found on PATH"
        checks.append(detail)

    container_available = any(
        shutil.which(name, path=environment.get("PATH"))
        for name in ("apptainer", "singularity")
    )
    checks.append(
        {
            "name": "runtime:apptainer_or_singularity",
            "ok": container_available,
            "detail": "available" if container_available else "not found on PATH",
        }
    )

    for name, url in STAGING_ENDPOINTS:
        request = urllib.request.Request(url, method="HEAD")
        check: dict[str, object] = {"name": f"endpoint:{name}", "url": url, "ok": False}
        try:
            with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
                check["status"] = response.status
                check["ok"] = True
        except urllib.error.HTTPError as exc:
            # Authentication failures still prove that the registry is reachable.
            check["status"] = exc.code
            check["ok"] = True
        except (urllib.error.URLError, OSError) as exc:
            check["detail"] = _redact(str(exc), secrets)
        checks.append(check)

    return {
        "ready": all(bool(check["ok"]) for check in checks),
        "network_mode": network_mode,
        "checks": checks,
        "proxy_values_recorded": False,
        "contains_secrets": False,
    }


def _stream_command(
    argv: list[str],
    *,
    environment: dict[str, str],
    secrets: list[str],
    step: int,
    steps: int,
    heartbeat_seconds: int,
) -> tuple[int, float]:
    started = time.monotonic()
    process = subprocess.Popen(
        argv,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        env=environment,
    )
    lines: queue.Queue[str | None] = queue.Queue()

    def read_output() -> None:
        assert process.stdout is not None
        for line in process.stdout:
            lines.put(line)
        lines.put(None)

    reader = threading.Thread(target=read_output, daemon=True)
    reader.start()
    last_heartbeat = started
    output_closed = False
    while process.poll() is None or not output_closed:
        try:
            line = lines.get(timeout=1)
            if line is None:
                output_closed = True
            else:
                print(_redact(line, secrets), end="", flush=True)
        except queue.Empty:
            pass
        now = time.monotonic()
        if process.poll() is None and now - last_heartbeat >= heartbeat_seconds:
            elapsed = int(now - started)
            print(
                f"[stage {step}/{steps}] running; elapsed={elapsed // 60:02d}:{elapsed % 60:02d}",
                file=sys.stderr,
                flush=True,
            )
            last_heartbeat = now
    reader.join()
    return process.wait(), time.monotonic() - started


def _discover_workflow_dir(output_root: Path) -> Path:
    """Return the unique shallowest downloaded Nextflow workflow directory."""

    candidates = sorted(
        path.parent
        for path in output_root.rglob("main.nf")
        if (path.parent / "nextflow.config").is_file()
    )
    if not candidates:
        raise StagingError(f"no downloaded workflow found under {output_root}")

    depths = {
        candidate: len(candidate.relative_to(output_root).parts)
        for candidate in candidates
    }
    shallowest_depth = min(depths.values())
    shallowest = [
        candidate
        for candidate, depth in depths.items()
        if depth == shallowest_depth
    ]
    if len(shallowest) != 1:
        raise StagingError(
            f"expected one top-level downloaded workflow under {output_root}, "
            f"found {len(shallowest)} at depth {shallowest_depth}"
        )
    return shallowest[0].resolve()


def run_staging_plan(
    plan_path: Path,
    *,
    heartbeat_seconds: int = 30,
) -> dict[str, object]:
    """Execute a validated plan once, with redacted output and a receipt."""

    if heartbeat_seconds < 1:
        raise StagingError("heartbeat interval must be at least 1 second")
    resolved_plan = plan_path.resolve(strict=True)
    try:
        plan = json.loads(resolved_plan.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise StagingError(f"could not read staging plan: {exc}") from exc
    if not isinstance(plan, dict) or plan.get("schema_version") != 1:
        raise StagingError("unsupported staging plan schema")
    if plan.get("stage") != "planned_not_executed":
        raise StagingError("staging plan is not executable")
    network_mode = plan.get("network_mode")
    if network_mode not in STAGING_NETWORK_MODES:
        raise StagingError("invalid staging network mode")

    lock_record = plan.get("workflow_lock")
    if not isinstance(lock_record, dict):
        raise StagingError("staging plan has no workflow lock")
    lock_path = Path(str(lock_record.get("path", ""))).resolve(strict=True)
    if sha256_file(lock_path) != lock_record.get("sha256"):
        raise StagingError("workflow lock changed after planning")
    lock = load_workflow_lock(lock_path)

    proxy_names = [name for name in PROXY_VARIABLES if os.environ.get(name)]
    if network_mode == "proxy" and not proxy_names:
        raise StagingError("proxy settings are no longer present")
    secrets = [os.environ[name] for name in proxy_names]
    environment = os.environ.copy()
    container_root = Path(str(plan.get("container_root", ""))).resolve()
    plugin_root = Path(str(plan.get("plugin_root", ""))).resolve()
    container_root.mkdir(parents=True, exist_ok=True)
    plugin_root.mkdir(parents=True, exist_ok=True)
    environment["NXF_SINGULARITY_CACHEDIR"] = str(container_root)
    environment["NXF_PLUGINS_DIR"] = str(plugin_root)

    for executable in ("nf-core", "nextflow"):
        if shutil.which(executable, path=environment.get("PATH")) is None:
            raise StagingError(f"required executable is not available: {executable}")
    if not any(shutil.which(name, path=environment.get("PATH")) for name in ("apptainer", "singularity")):
        raise StagingError("apptainer or singularity is required for container staging")
    version_commands = {
        "nf-core/tools": (["nf-core", "--version"], lock.runtime.nf_core_tools_version),
        "Nextflow": (["nextflow", "-version"], lock.runtime.nextflow_version),
    }
    for label, (version_argv, expected_version) in version_commands.items():
        try:
            version_output = subprocess.check_output(
                version_argv,
                stderr=subprocess.STDOUT,
                text=True,
                env=environment,
            )
        except subprocess.CalledProcessError as exc:
            raise StagingError(f"could not determine {label} version") from exc
        if re.search(
            rf"(?<![0-9.]){re.escape(expected_version)}(?![0-9.])",
            version_output,
        ) is None:
            raise StagingError(
                f"{label} must be exactly {expected_version}; observed: "
                f"{_redact(version_output.strip(), secrets)}"
            )

    commands = plan.get("commands")
    if not isinstance(commands, list) or len(commands) != 2:
        raise StagingError("staging plan must contain exactly two workflow downloads")
    expected_specs = (lock.rnaseq, lock.differential)
    validated: list[list[str]] = []
    for record, spec in zip(commands, expected_specs, strict=True):
        if not isinstance(record, dict):
            raise StagingError("invalid staging command record")
        expected = build_download_argv(
            spec,
            outdir=Path(str(record.get("output_root", ""))),
            parallel_downloads=int(plan.get("parallel_downloads", 0)),
        )
        if record.get("workflow") != asdict(spec) or record.get("argv") != expected:
            raise StagingError("staging command does not match the pinned workflow lock")
        if Path(str(record.get("output_root", ""))).exists():
            raise StagingError(
                f"refusing existing staging target: {record.get('output_root')}"
            )
        validated.append(expected)

    receipt_path = resolved_plan.with_name(resolved_plan.stem + ".receipt.json")
    receipt: dict[str, object] = {
        "schema_version": 1,
        "started_at": datetime.now(timezone.utc).isoformat(),
        "status": "reserved",
        "plan_sha256": sha256_file(resolved_plan),
        "steps": [],
        "proxy_values_recorded": False,
        "contains_secrets": False,
    }
    _exclusive_json(receipt_path, receipt)

    try:
        for index, argv in enumerate(validated, start=1):
            print(f"[stage {index}/{len(validated)}] starting {argv[3]}", file=sys.stderr)
            returncode, wall_seconds = _stream_command(
                argv,
                environment=environment,
                secrets=secrets,
                step=index,
                steps=len(validated),
                heartbeat_seconds=heartbeat_seconds,
            )
            receipt["steps"].append(
                {
                    "workflow": argv[3],
                    "returncode": returncode,
                    "wall_seconds": wall_seconds,
                }
            )
            if returncode != 0:
                raise StagingError(
                    f"workflow staging failed for {argv[3]} with exit {returncode}"
                )
            output_root = Path(str(commands[index - 1]["output_root"]))
            receipt["steps"][-1]["workflow_dir"] = str(
                _discover_workflow_dir(output_root)
            )
        receipt["status"] = "complete"
    except Exception:
        receipt["status"] = "failed"
        raise
    finally:
        receipt["finished_at"] = datetime.now(timezone.utc).isoformat()
        temporary = receipt_path.with_suffix(receipt_path.suffix + ".tmp")
        temporary.write_text(
            json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.chmod(temporary, 0o600)
        temporary.replace(receipt_path)
    return receipt
