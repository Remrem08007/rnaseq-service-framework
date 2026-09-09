from __future__ import annotations

import stat
from pathlib import Path

import pytest

from rnaseq_service.hpc_config import (
    HPCConfigError,
    load_hpc_settings,
    render_nextflow_config,
    write_nextflow_config,
)


def settings_file(tmp_path: Path, *, account: str = "project-123") -> Path:
    path = tmp_path / "infrastructure.toml"
    path.write_text(
        f"""[cluster]
account = "{account}"
partition = "compute"
job_name = "rnaseq-service"
launcher_cpus = 2
launcher_memory_gb = 8
launcher_time = "12:00:00"
max_task_cpus = 32
max_task_memory_gb = 128
max_task_time_hours = 24
queue_size = 100
per_cpu_memory = false

[paths]
work_root = "{tmp_path / 'work'}"
container_cache = "{tmp_path / 'containers'}"
reference_cache = "{tmp_path / 'references'}"

[software]
purge_modules = true
modules = ["StdEnv/2023", "nextflow/26.04.4", "apptainer/1.3.5"]
""",
        encoding="utf-8",
    )
    return path


def test_rendered_config_is_generic_slurm_apptainer(tmp_path: Path) -> None:
    settings = load_hpc_settings(settings_file(tmp_path))

    rendered = render_nextflow_config(settings)

    assert "executor = 'slurm'" in rendered
    assert "queue = 'compute'" in rendered
    assert "clusterOptions = '--account=project-123'" in rendered
    assert "perCpuMemAllocation = false" in rendered
    assert "apptainer" in rendered
    assert "max_memory = 128.GB" in rendered
    assert "max_time = 24.h" in rendered
    assert str(tmp_path / "references") in rendered


def test_config_is_private_and_never_overwritten(tmp_path: Path) -> None:
    source = settings_file(tmp_path)
    output = tmp_path / "rendered" / "nextflow.config"

    write_nextflow_config(source, output)

    assert stat.S_IMODE(output.stat().st_mode) == 0o600
    with pytest.raises(HPCConfigError, match="refusing to overwrite"):
        write_nextflow_config(source, output)


def test_shell_like_cluster_value_is_rejected(tmp_path: Path) -> None:
    source = settings_file(tmp_path, account="project;touch-bad")

    with pytest.raises(HPCConfigError, match="unsafe account"):
        load_hpc_settings(source)


def test_relative_shared_path_is_rejected(tmp_path: Path) -> None:
    source = settings_file(tmp_path)
    source.write_text(
        source.read_text(encoding="utf-8").replace(
            f'work_root = "{tmp_path / "work"}"', 'work_root = "relative/work"'
        ),
        encoding="utf-8",
    )

    with pytest.raises(HPCConfigError, match="absolute path"):
        load_hpc_settings(source)
