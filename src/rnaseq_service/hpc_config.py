"""Validate external HPC settings and render a generic Nextflow configuration."""

from __future__ import annotations

import os
import re
import tomllib
from dataclasses import asdict, dataclass
from pathlib import Path


TOKEN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
MODULE_TOKEN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.+/-]*$")
SLURM_TIME = re.compile(r"^(?:\d+-)?\d{1,3}:[0-5]\d:[0-5]\d$")


class HPCConfigError(ValueError):
    """Raised when infrastructure settings are missing or unsafe."""


@dataclass(frozen=True)
class ClusterSettings:
    account: str
    partition: str
    job_name: str
    launcher_cpus: int
    launcher_memory_gb: int
    launcher_time: str
    max_task_cpus: int
    max_task_memory_gb: int
    max_task_time_hours: int
    queue_size: int
    per_cpu_memory: bool


@dataclass(frozen=True)
class PathSettings:
    work_root: str
    container_cache: str
    reference_cache: str | None


@dataclass(frozen=True)
class SoftwareSettings:
    purge_modules: bool
    modules: tuple[str, ...]


@dataclass(frozen=True)
class HPCSettings:
    cluster: ClusterSettings
    paths: PathSettings
    software: SoftwareSettings


def _table(payload: dict[str, object], name: str, path: Path) -> dict[str, object]:
    value = payload.get(name)
    if not isinstance(value, dict):
        raise HPCConfigError(f"{path}: [{name}] table is required")
    return value


def _string(table: dict[str, object], key: str, path: Path) -> str:
    value = table.get(key)
    if not isinstance(value, str) or not value.strip():
        raise HPCConfigError(f"{path}: {key!r} must be a non-empty string")
    return value.strip()


def _positive_int(table: dict[str, object], key: str, path: Path) -> int:
    value = table.get(key)
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise HPCConfigError(f"{path}: {key!r} must be a positive integer")
    return value


def _absolute_path(table: dict[str, object], key: str, path: Path) -> str:
    value = Path(_string(table, key, path)).expanduser()
    if not value.is_absolute():
        raise HPCConfigError(f"{path}: {key!r} must be an absolute path")
    return str(value.resolve())


def load_hpc_settings(path: Path) -> HPCSettings:
    """Load a private TOML contract; no institution defaults are assumed."""

    try:
        with path.open("rb") as handle:
            payload = tomllib.load(handle)
    except FileNotFoundError as exc:
        raise HPCConfigError(f"{path}: settings file does not exist") from exc
    except (OSError, tomllib.TOMLDecodeError) as exc:
        raise HPCConfigError(f"{path}: could not read settings: {exc}") from exc

    cluster = _table(payload, "cluster", path)
    paths = _table(payload, "paths", path)
    software = _table(payload, "software", path)
    account = _string(cluster, "account", path)
    partition = _string(cluster, "partition", path)
    job_name = _string(cluster, "job_name", path)
    for label, value in (
        ("account", account),
        ("partition", partition),
        ("job_name", job_name),
    ):
        if TOKEN.fullmatch(value) is None:
            raise HPCConfigError(f"{path}: unsafe {label}: {value!r}")
    launcher_time = _string(cluster, "launcher_time", path)
    if SLURM_TIME.fullmatch(launcher_time) is None:
        raise HPCConfigError(
            f"{path}: launcher_time must use [days-]hours:minutes:seconds"
        )
    per_cpu_memory = cluster.get("per_cpu_memory")
    if not isinstance(per_cpu_memory, bool):
        raise HPCConfigError(f"{path}: per_cpu_memory must be true or false")

    raw_modules = software.get("modules")
    if not isinstance(raw_modules, list) or not raw_modules:
        raise HPCConfigError(f"{path}: software.modules must be a non-empty list")
    modules: list[str] = []
    for value in raw_modules:
        if not isinstance(value, str) or MODULE_TOKEN.fullmatch(value) is None:
            raise HPCConfigError(f"{path}: unsafe module name: {value!r}")
        modules.append(value)
    purge_modules = software.get("purge_modules")
    if not isinstance(purge_modules, bool):
        raise HPCConfigError(f"{path}: purge_modules must be true or false")

    reference_value = paths.get("reference_cache")
    reference_cache: str | None
    if reference_value is None:
        reference_cache = None
    else:
        reference_cache = _absolute_path(paths, "reference_cache", path)

    return HPCSettings(
        cluster=ClusterSettings(
            account=account,
            partition=partition,
            job_name=job_name,
            launcher_cpus=_positive_int(cluster, "launcher_cpus", path),
            launcher_memory_gb=_positive_int(cluster, "launcher_memory_gb", path),
            launcher_time=launcher_time,
            max_task_cpus=_positive_int(cluster, "max_task_cpus", path),
            max_task_memory_gb=_positive_int(cluster, "max_task_memory_gb", path),
            max_task_time_hours=_positive_int(cluster, "max_task_time_hours", path),
            queue_size=_positive_int(cluster, "queue_size", path),
            per_cpu_memory=per_cpu_memory,
        ),
        paths=PathSettings(
            work_root=_absolute_path(paths, "work_root", path),
            container_cache=_absolute_path(paths, "container_cache", path),
            reference_cache=reference_cache,
        ),
        software=SoftwareSettings(
            purge_modules=purge_modules,
            modules=tuple(modules),
        ),
    )


def _groovy_string(value: str) -> str:
    return "'" + value.replace("\\", "\\\\").replace("'", "\\'") + "'"


def render_nextflow_config(settings: HPCSettings) -> str:
    """Render only validated literals into a portable SLURM/Apptainer config."""

    cluster = settings.cluster
    paths = settings.paths
    lines = [
        "// Generated by rnaseq-service-framework; do not edit in place.",
        "process {",
        "    executor = 'slurm'",
        f"    queue = {_groovy_string(cluster.partition)}",
        f"    clusterOptions = {_groovy_string('--account=' + cluster.account)}",
        "}",
        "",
        "executor {",
        f"    queueSize = {cluster.queue_size}",
        f"    perCpuMemAllocation = {str(cluster.per_cpu_memory).lower()}",
        "}",
        "",
        "apptainer {",
        "    enabled = true",
        "    autoMounts = true",
        f"    cacheDir = {_groovy_string(paths.container_cache)}",
        "}",
        "",
        "params {",
        f"    max_cpus = {cluster.max_task_cpus}",
        f"    max_memory = {cluster.max_task_memory_gb}.GB",
        f"    max_time = {cluster.max_task_time_hours}.h",
    ]
    if paths.reference_cache is not None:
        lines.append(f"    igenomes_base = {_groovy_string(paths.reference_cache)}")
    lines.extend(["}", ""])
    return "\n".join(lines)


def write_nextflow_config(settings_path: Path, output: Path) -> dict[str, object]:
    """Exclusively write a generated config and return its public-safe summary."""

    settings = load_hpc_settings(settings_path)
    rendered = render_nextflow_config(settings)
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        descriptor = os.open(output, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    except FileExistsError as exc:
        raise HPCConfigError(f"refusing to overwrite existing config: {output}") from exc
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(rendered)
    except Exception:
        output.unlink(missing_ok=True)
        raise
    return {
        "created": True,
        "output": str(output.resolve()),
        "settings": asdict(settings),
        "contains_secrets": False,
    }
