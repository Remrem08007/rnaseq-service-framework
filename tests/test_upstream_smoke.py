from __future__ import annotations

import json
import stat
from pathlib import Path

import pytest

from rnaseq_service.upstream_smoke import UpstreamSmokeError, create_upstream_smoke_plan
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
