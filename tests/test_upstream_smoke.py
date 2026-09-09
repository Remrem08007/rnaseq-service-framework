from __future__ import annotations

import json
import stat
from pathlib import Path

import pytest

from rnaseq_service.upstream_smoke import (
    UpstreamSmokeError,
    create_upstream_smoke_launcher,
    create_upstream_smoke_plan,
    create_upstream_smoke_receipt,
    inspect_upstream_smoke_completion,
)
from rnaseq_service.upstream_smoke_cli import main as smoke_main


ROOT = Path(__file__).parents[1]


def _file(path: Path, content: str = "evidence\n") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _complete_fixture(plan: dict[str, object]) -> None:
    run_root = Path(str(plan["execution"]["run_root"]))
    trace = (
        "task_id\tname\tstatus\trealtime\t%cpu\tpeak_rss\tduration\n"
        "1\tTASK\tCOMPLETED\t1s\t100%\t1 MB\t1s\n"
    )
    for stage_id in ("rnaseq", "differential"):
        stage_root = run_root / stage_id
        execution = stage_root / "execution"
        results = stage_root / "results"
        _file(execution / "trace.tsv", trace)
        for name in ("report.html", "timeline.html", "dag.html"):
            _file(execution / name)
        _file(
            results / "pipeline_info" / "params.json",
            json.dumps({"outdir": str(results.resolve())}),
        )
        version_name = (
            "nf_core_rnaseq_software_mqc_versions.yml"
            if stage_id == "rnaseq"
            else "nf_core_differentialabundance_software_versions.yml"
        )
        _file(results / "pipeline_info" / version_name)
    rnaseq = run_root / "rnaseq" / "results"
    _file(rnaseq / "multiqc" / "star_salmon" / "multiqc_report.html")
    _file(rnaseq / "star_salmon" / "salmon.merged.gene_counts.tsv")
    _file(rnaseq / "star_salmon" / "salmon.merged.gene_tpm.tsv")
    differential = run_root / "differential" / "results"
    prefix = "treatment_mCherry_hND6__SRP254919"
    _file(
        differential
        / "report"
        / "rnaseq_deseq2_gsea"
        / "SRP254919_differentialabundance_report.html"
    )
    _file(
        differential
        / "tables"
        / "processed_abundance"
        / "rnaseq_deseq2_gsea"
        / "all.normalised_counts.tsv"
    )
    _file(
        differential
        / "tables"
        / "processed_abundance"
        / "rnaseq_deseq2_gsea"
        / "all.vst.tsv"
    )
    _file(
        differential
        / "tables"
        / "differential"
        / "rnaseq_deseq2_gsea"
        / f"{prefix}.deseq2.results.tsv"
    )
    _file(
        differential
        / "tables"
        / "differential"
        / "rnaseq_deseq2_gsea"
        / f"{prefix}.deseq2.results_filtered.tsv"
    )
    _file(
        differential
        / "plots"
        / "differential"
        / "rnaseq_deseq2_gsea"
        / prefix
        / "png"
        / "volcano.png"
    )
    _file(
        differential
        / "report"
        / "gsea"
        / "rnaseq_deseq2_gsea"
        / prefix
        / "hallmarks"
        / f"{prefix}.hallmarks.Gsea.rpt"
    )


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


def test_completion_requires_both_successful_scientific_stages(tmp_path: Path) -> None:
    plan_path = tmp_path / "plan.json"
    plan = create_upstream_smoke_plan(
        workflow_lock=ROOT / "config" / "workflows.toml",
        output=plan_path,
        run_root=tmp_path / "runs",
        network_mode="direct",
    )
    _complete_fixture(plan)
    progress: list[tuple[str, int, int]] = []

    result = inspect_upstream_smoke_completion(
        plan_path,
        progress=lambda role, current, total: progress.append((role, current, total)),
    )

    assert result["stage"] == "upstream_smoke_complete"
    assert result["scientific_execution_performed"] is True
    assert result["uses_public_upstream_test_data"] is True
    assert result["contains_client_data"] is False
    assert [stage["id"] for stage in result["stages"]] == ["rnaseq", "differential"]
    assert all(stage["task_summary"]["task_count"] == 1 for stage in result["stages"])
    assert {artifact["role"] for artifact in result["stages"][0]["artifacts"]} >= {
        "multiqc_report",
        "gene_counts",
        "gene_tpm",
    }
    assert {artifact["role"] for artifact in result["stages"][1]["artifacts"]} >= {
        "deseq2_results",
        "volcano_plot",
        "gsea_report",
    }
    assert progress


def test_completion_receipt_is_immutable_and_private(tmp_path: Path) -> None:
    plan_path = tmp_path / "plan.json"
    plan = create_upstream_smoke_plan(
        workflow_lock=ROOT / "config" / "workflows.toml",
        output=plan_path,
        run_root=tmp_path / "runs",
        network_mode="proxy",
    )
    _complete_fixture(plan)
    receipt = tmp_path / "private" / "upstream-smoke.complete.json"

    create_upstream_smoke_receipt(run_plan=plan_path, output=receipt)

    assert stat.S_IMODE(receipt.stat().st_mode) == 0o600
    with pytest.raises(UpstreamSmokeError, match="refusing to overwrite"):
        create_upstream_smoke_receipt(run_plan=plan_path, output=receipt)


def test_completion_rejects_nonterminal_trace_and_missing_science(tmp_path: Path) -> None:
    plan_path = tmp_path / "plan.json"
    plan = create_upstream_smoke_plan(
        workflow_lock=ROOT / "config" / "workflows.toml",
        output=plan_path,
        run_root=tmp_path / "runs",
        network_mode="direct",
    )
    _complete_fixture(plan)
    trace = Path(str(plan["execution"]["run_root"])) / "rnaseq" / "execution" / "trace.tsv"
    trace.write_text(
        "task_id\tname\tstatus\trealtime\t%cpu\tpeak_rss\tduration\n"
        "1\tTASK\tRUNNING\t1s\t100%\t1 MB\t1s\n",
        encoding="utf-8",
    )
    with pytest.raises(UpstreamSmokeError, match="not complete and successful"):
        inspect_upstream_smoke_completion(plan_path)

    trace.write_text(
        "task_id\tname\tstatus\trealtime\t%cpu\tpeak_rss\tduration\n"
        "1\tTASK\tCOMPLETED\t1s\t100%\t1 MB\t1s\n",
        encoding="utf-8",
    )
    gsea = next(
        (Path(str(plan["execution"]["run_root"])) / "differential" / "results").glob(
            "report/gsea/**/*.Gsea.rpt"
        )
    )
    gsea.unlink()
    with pytest.raises(UpstreamSmokeError, match="GSEA report"):
        inspect_upstream_smoke_completion(plan_path)


def test_complete_cli_reports_real_upstream_scope(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    plan_path = tmp_path / "plan.json"
    plan = create_upstream_smoke_plan(
        workflow_lock=ROOT / "config" / "workflows.toml",
        output=plan_path,
        run_root=tmp_path / "runs",
        network_mode="direct",
    )
    _complete_fixture(plan)
    output = tmp_path / "receipt.json"

    assert smoke_main(
        ["complete", "--run-plan", str(plan_path), "--output", str(output), "--quiet"]
    ) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["status"] == "upstream_smoke_complete"
    assert summary["scientific_execution_performed"] is True
    assert [stage["id"] for stage in summary["stages"]] == ["rnaseq", "differential"]
