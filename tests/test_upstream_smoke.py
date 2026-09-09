from __future__ import annotations

import json
import stat
from pathlib import Path

import pytest

from rnaseq_service.upstream_smoke import (
    UpstreamSmokeError,
    create_upstream_smoke_launcher,
    create_upstream_smoke_plan,
)
from rnaseq_service.upstream_smoke_cli import main as smoke_main


ROOT = Path(__file__).parents[1]


def test_plan_pins_both_official_public_test_profiles(tmp_path: Path) -> None:
    output = tmp_path / "private" / "upstream-smoke.json"
    plan = create_upstream_smoke_plan(
        workflow_lock=ROOT / "config" / "workflows.toml",
        output=output,
        run_root=tmp_path / "runs",
        network_mode="direct",
    )

    assert plan["stage"] == "upstream_smoke_planned_not_executed"
    assert plan["uses_public_upstream_test_data"] is True
    assert plan["scientific_execution_expected"] is True
    assert plan["contains_client_data"] is False
    assert plan["contains_secrets"] is False
    assert [stage["id"] for stage in plan["stages"]] == ["rnaseq", "differential"]
    rnaseq, differential = plan["stages"]
    assert rnaseq["workflow"] == {"name": "nf-core/rnaseq", "revision": "3.26.0"}
    assert differential["workflow"] == {
        "name": "nf-core/differentialabundance",
        "revision": "2.0.0",
    }
    assert "test,apptainer" in rnaseq["command_argv"]
    assert "test_rnaseq_deseq2_gsea,apptainer" in differential["command_argv"]
    assert all("-resume" in stage["command_argv"] for stage in plan["stages"])
    assert stat.S_IMODE(output.stat().st_mode) == 0o600


def test_plan_rejects_offline_and_overwrite(tmp_path: Path) -> None:
    kwargs = {
        "workflow_lock": ROOT / "config" / "workflows.toml",
        "output": tmp_path / "plan.json",
        "run_root": tmp_path / "runs",
        "network_mode": "offline",
    }
    with pytest.raises(UpstreamSmokeError, match="direct.*proxy"):
        create_upstream_smoke_plan(**kwargs)
    kwargs["network_mode"] = "direct"
    create_upstream_smoke_plan(**kwargs)
    with pytest.raises(UpstreamSmokeError, match="refusing to overwrite"):
        create_upstream_smoke_plan(**kwargs)


def test_upstream_smoke_cli_plans_without_running(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    output = tmp_path / "plan.json"
    assert smoke_main(
        [
            "plan",
            "--workflow-lock",
            str(ROOT / "config" / "workflows.toml"),
            "--output",
            str(output),
            "--run-root",
            str(tmp_path / "runs"),
            "--network-mode",
            "direct",
        ]
    ) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary == {
        "scientific_execution_expected": True,
        "stages": ["rnaseq", "differential"],
        "status": "upstream_smoke_planned_not_executed",
        "uses_public_upstream_test_data": True,
    }


def test_prepare_renders_non_submitting_heartbeat_launcher(tmp_path: Path) -> None:
    plan_path = tmp_path / "private" / "plan.json"
    create_upstream_smoke_plan(
        workflow_lock=ROOT / "config" / "workflows.toml",
        output=plan_path,
        run_root=tmp_path / "runs",
        network_mode="direct",
    )
    launcher = tmp_path / "private" / "launch" / "smoke.sbatch"

    result = create_upstream_smoke_launcher(
        run_plan=plan_path,
        output=launcher,
        account="def-project",
        partition="compute",
        wall_time="04:00:00",
        cpus=8,
        memory_gb=32,
        modules=["StdEnv/2023", "nextflow/26.04.6", "apptainer/1.3.5"],
    )

    text = launcher.read_text(encoding="utf-8")
    assert result["prepared"] is True
    assert result["submitted"] is False
    assert result["heartbeat_seconds"] == 60
    assert "#SBATCH --account=def-project" in text
    assert "stage=${ordinal}/2" in text
    assert "test,apptainer" in text
    assert "test_rnaseq_deseq2_gsea,apptainer" in text
    assert text.count("run_stage ") == 2
    assert "sbatch" not in text
    assert stat.S_IMODE(launcher.stat().st_mode) == 0o700


def test_prepare_rejects_changed_plan_control_and_overwrite(tmp_path: Path) -> None:
    plan_path = tmp_path / "plan.json"
    lock = tmp_path / "workflows.toml"
    lock.write_text(
        (ROOT / "config" / "workflows.toml").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    create_upstream_smoke_plan(
        workflow_lock=lock,
        output=plan_path,
        run_root=tmp_path / "runs",
        network_mode="direct",
    )
    kwargs = {
        "run_plan": plan_path,
        "output": tmp_path / "smoke.sbatch",
        "account": "def-project",
        "partition": "compute",
        "wall_time": "04:00:00",
        "cpus": 8,
        "memory_gb": 32,
        "modules": ["nextflow", "apptainer"],
    }
    create_upstream_smoke_launcher(**kwargs)
    with pytest.raises(UpstreamSmokeError, match="refusing to overwrite"):
        create_upstream_smoke_launcher(**kwargs)

    original = lock.read_text(encoding="utf-8")
    try:
        lock.write_text(original + "\n", encoding="utf-8")
        with pytest.raises(UpstreamSmokeError, match="workflow lock .* changed"):
            create_upstream_smoke_launcher(**{**kwargs, "output": tmp_path / "other.sbatch"})
    finally:
        lock.write_text(original, encoding="utf-8")


def test_prepare_cli(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    plan_path = tmp_path / "plan.json"
    create_upstream_smoke_plan(
        workflow_lock=ROOT / "config" / "workflows.toml",
        output=plan_path,
        run_root=tmp_path / "runs",
        network_mode="proxy",
    )
    launcher = tmp_path / "smoke.sbatch"
    assert smoke_main(
        [
            "prepare",
            "--run-plan",
            str(plan_path),
            "--output",
            str(launcher),
            "--account",
            "def-project",
            "--partition",
            "compute",
            "--time",
            "04:00:00",
            "--module",
            "nextflow",
            "--module",
            "apptainer",
        ]
    ) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["prepared"] is True
    assert summary["submitted"] is False
