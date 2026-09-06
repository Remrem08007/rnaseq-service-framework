from __future__ import annotations

import csv
import json
import stat
import subprocess
from pathlib import Path

import pytest

from rnaseq_service.hpc_config import write_nextflow_config
from rnaseq_service.hpc_run import HPCRunError, prepare_launcher, query_job, submit_launcher
from rnaseq_service.inputs import create_input_manifest
from rnaseq_service.plan import create_rnaseq_plan
from rnaseq_service.primary import PrimaryDeliveryError, prepare_safe_restart


ROOT = Path(__file__).parents[1]


def write_csv(path: Path, headers: list[str], rows: list[list[str]]) -> Path:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        writer.writerows(rows)
    return path


def create_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    samples = []
    design_rows = []
    for condition in ("control", "treated"):
        for replicate in (1, 2):
            sample = f"{condition}_{replicate}"
            r1 = tmp_path / f"{sample}_R1.fastq.gz"
            r2 = tmp_path / f"{sample}_R2.fastq.gz"
            r1.write_bytes(b"reads")
            r2.write_bytes(b"reads")
            samples.append([sample, str(r1), str(r2), "auto"])
            design_rows.append([sample, condition])
    sample_sheet = write_csv(
        tmp_path / "samplesheet.csv",
        ["sample", "fastq_1", "fastq_2", "strandedness"],
        samples,
    )
    design = write_csv(tmp_path / "design.csv", ["sample", "condition"], design_rows)
    contrasts = write_csv(
        tmp_path / "contrasts.csv",
        ["contrast_id", "variable", "reference", "target"],
        [["treated_vs_control", "condition", "control", "treated"]],
    )
    return sample_sheet, design, contrasts


def make_settings(tmp_path: Path) -> Path:
    settings = tmp_path / "infrastructure.toml"
    settings.write_text(
        f"""[cluster]
account = "project-123"
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
work_root = "{tmp_path / 'shared-work'}"
container_cache = "{tmp_path / 'containers'}"
[software]
purge_modules = true
modules = ["nextflow/26.04.6", "apptainer/1.3.5"]
""",
        encoding="utf-8",
    )
    return settings


def make_run_plan(tmp_path: Path) -> tuple[Path, Path]:
    settings = make_settings(tmp_path)
    config = tmp_path / "generated" / "nextflow.config"
    write_nextflow_config(settings, config)
    samplesheet, design, contrasts = create_inputs(tmp_path)
    input_manifest = tmp_path / "input-manifest.json"
    create_input_manifest(samplesheet=samplesheet, output=input_manifest)
    plan_path = tmp_path / "plans" / "run.json"
    create_rnaseq_plan(
        samplesheet=samplesheet,
        design=design,
        contrasts=contrasts,
        workflow_lock=ROOT / "config" / "workflows.toml",
        output=plan_path,
        outdir=tmp_path / "results",
        workdir=tmp_path / "shared-work" / "study-1",
        network_mode="direct",
        container_engine="apptainer",
        executor="slurm",
        infrastructure_config=config,
        input_manifest=input_manifest,
    )
    return plan_path, settings


def test_prepare_launcher_verifies_plan_and_renders_resume(tmp_path: Path) -> None:
    plan, settings = make_run_plan(tmp_path)
    launcher = tmp_path / "launch" / "controller.sbatch"

    result = prepare_launcher(run_plan=plan, settings_path=settings, output=launcher)
    rendered = launcher.read_text(encoding="utf-8")

    assert result["resumable"] is True
    assert "#SBATCH --account=project-123" in rendered
    assert "#SBATCH --partition=compute" in rendered
    assert "module load nextflow/26.04.6" in rendered
    assert "nextflow run" in rendered
    assert "-resume" in rendered
    assert stat.S_IMODE(launcher.stat().st_mode) == 0o700
    assert (launcher.parent / "logs").is_dir()


def test_prepare_rejects_changed_control_file(tmp_path: Path) -> None:
    plan, settings = make_run_plan(tmp_path)
    payload = json.loads(plan.read_text(encoding="utf-8"))
    config = Path(payload["control_files"]["infrastructure_config"]["path"])
    config.write_text(config.read_text(encoding="utf-8") + "// changed\n", encoding="utf-8")

    with pytest.raises(HPCRunError, match="changed after planning"):
        prepare_launcher(
            run_plan=plan,
            settings_path=settings,
            output=tmp_path / "controller.sbatch",
        )


def test_submit_reserves_receipt_and_records_job_id(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    launcher = tmp_path / "controller.sbatch"
    launcher.write_text("#!/bin/bash\ntrue\n", encoding="utf-8")
    receipt = tmp_path / "submission.json"
    monkeypatch.setattr(
        "rnaseq_service.hpc_run.subprocess.run",
        lambda *args, **kwargs: subprocess.CompletedProcess(args[0], 0, "123456;cluster\n", ""),
    )

    result = submit_launcher(launcher=launcher, receipt=receipt)

    assert result["status"] == "submitted"
    assert result["job_id"] == "123456"
    assert stat.S_IMODE(receipt.stat().st_mode) == 0o600
    with pytest.raises(HPCRunError, match="existing submission receipt"):
        submit_launcher(launcher=launcher, receipt=receipt)


def test_status_uses_squeue_for_active_job(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    receipt = tmp_path / "submission.json"
    receipt.write_text(
        json.dumps({"status": "submitted", "job_id": "123456"}), encoding="utf-8"
    )
    monkeypatch.setattr(
        "rnaseq_service.hpc_run.subprocess.run",
        lambda *args, **kwargs: subprocess.CompletedProcess(
            args[0], 0, "RUNNING|01:02|12:00:00|node01\n", ""
        ),
    )

    result = query_job(receipt)

    assert result["source"] == "squeue"
    assert result["state"] == "RUNNING"
    assert result["terminal"] is False


def test_status_falls_back_to_sacct_for_completed_job(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    receipt = tmp_path / "submission.json"
    receipt.write_text(
        json.dumps({"status": "submitted", "job_id": "123456"}), encoding="utf-8"
    )

    def fake_run(argv, **kwargs):
        if argv[0] == "squeue":
            return subprocess.CompletedProcess(argv, 0, "", "")
        return subprocess.CompletedProcess(
            argv, 0, "123456|COMPLETED|02:00:00|05:00:00|8G|0:0\n", ""
        )

    monkeypatch.setattr("rnaseq_service.hpc_run.subprocess.run", fake_run)

    result = query_job(receipt)

    assert result["source"] == "sacct"
    assert result["terminal"] is True
    assert result["successful"] is True


def restart_fixture(tmp_path: Path) -> tuple[Path, Path, Path]:
    plan, settings = make_run_plan(tmp_path)
    payload = json.loads(plan.read_text(encoding="utf-8"))
    Path(payload["execution"]["workdir"]).mkdir(parents=True)
    launcher = tmp_path / "launch-1" / "controller.sbatch"
    launch = prepare_launcher(run_plan=plan, settings_path=settings, output=launcher)
    receipt = tmp_path / "launch-1" / "submission.json"
    receipt.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "status": "submitted",
                "job_id": "123456",
                "launcher": str(launcher.resolve()),
                "launcher_sha256": launch["launcher_sha256"],
            }
        ),
        encoding="utf-8",
    )
    return plan, settings, receipt


def test_safe_restart_requires_terminal_failure_and_reuses_workdir(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    plan, settings, receipt = restart_fixture(tmp_path)

    def fake_run(argv, **kwargs):
        if argv[0] == "squeue":
            return subprocess.CompletedProcess(argv, 0, "", "")
        return subprocess.CompletedProcess(
            argv, 0, "123456|OUT_OF_MEMORY|01:00:00|02:00:00|8G|0:9\n", ""
        )

    monkeypatch.setattr("rnaseq_service.hpc_run.subprocess.run", fake_run)
    output = tmp_path / "launch-2" / "controller.sbatch"

    result = prepare_safe_restart(
        run_plan=plan,
        settings_path=settings,
        previous_receipt=receipt,
        output=output,
    )

    assert result["restart_prepared"] is True
    assert result["previous_scheduler_state"] == "OUT_OF_MEMORY"
    assert result["submitted"] is False
    assert "-resume" in output.read_text(encoding="utf-8")


def test_safe_restart_refuses_active_job(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    plan, settings, receipt = restart_fixture(tmp_path)
    monkeypatch.setattr(
        "rnaseq_service.hpc_run.subprocess.run",
        lambda argv, **kwargs: subprocess.CompletedProcess(
            argv, 0, "RUNNING|00:05|12:00:00|node01\n", ""
        ),
    )

    with pytest.raises(PrimaryDeliveryError, match="not terminal"):
        prepare_safe_restart(
            run_plan=plan,
            settings_path=settings,
            previous_receipt=receipt,
            output=tmp_path / "launch-2" / "controller.sbatch",
        )
