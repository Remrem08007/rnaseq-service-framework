from __future__ import annotations

import csv
import json
import stat
from pathlib import Path

import pytest

from rnaseq_service.differential import (
    DifferentialHandoffError,
    create_differential_plan,
    create_differential_receipt,
)
from rnaseq_service.differential_cli import main as differential_main
from rnaseq_service.qc import evaluate_qc, finalize_qc
from test_qc import qc_fixture, write_reviewed_decisions


def accepted_fixture(tmp_path: Path) -> Path:
    completion, policy = qc_fixture(tmp_path)
    assessment_dir = tmp_path / "assessment"
    evaluate_qc(completion_receipt=completion, policy_path=policy, outdir=assessment_dir)
    decisions = tmp_path / "reviewed.tsv"
    write_reviewed_decisions(decisions)
    accepted_dir = tmp_path / "accepted"
    finalize_qc(
        assessment_path=assessment_dir / "qc_assessment.json",
        decisions_path=decisions,
        outdir=accepted_dir,
    )
    return accepted_dir / "qc_acceptance.json"


def plan_kwargs(tmp_path: Path) -> dict[str, object]:
    return {
        "qc_acceptance": accepted_fixture(tmp_path),
        "handoff_dir": tmp_path / "handoff",
        "output": tmp_path / "private" / "differential-plan.json",
        "outdir": tmp_path / "results" / "differential",
        "workdir": tmp_path / "work" / "differential",
        "study_name": "study001",
        "network_mode": "direct",
        "container_engine": "apptainer",
        "executor": "local",
        "seed": 20260908,
    }


def completed_differential_fixture(tmp_path: Path) -> tuple[Path, Path]:
    kwargs = plan_kwargs(tmp_path)
    plan = create_differential_plan(**kwargs)
    run_root = Path(kwargs["outdir"])
    execution = run_root / "execution"
    execution.mkdir(parents=True)
    for name in ("report.html", "timeline.html", "dag.html"):
        (execution / name).write_text(f"<html>{name}</html>\n", encoding="utf-8")
    (execution / "trace.tsv").write_text(
        "name\tstatus\tduration\trealtime\t%cpu\tpeak_rss\n"
        "DESEQ2_NORM\tCOMPLETED\t2s\t1s\t90%\t100 MB\n"
        "DESEQ2_DIFFERENTIAL\tCACHED\t3s\t2s\t95%\t200 MB\n",
        encoding="utf-8",
    )
    pipeline = run_root / "differentialabundance"
    info = pipeline / "pipeline_info"
    info.mkdir(parents=True)
    (info / "nf_core_differentialabundance_software_versions.yml").write_text(
        "DESeq2: 1.34.0\n", encoding="utf-8"
    )
    controls = plan["control_files"]
    (info / "params.json").write_text(
        json.dumps(
            {
                "input": controls["observations"]["path"],
                "matrix": controls["gene_counts"]["path"],
                "feature_length_matrix": controls["gene_lengths"]["path"],
                "contrasts": controls["contrasts"]["path"],
                "outdir": str(pipeline.resolve()),
                "study_name": "study001",
                "study_type": "rnaseq",
                "observations_id_col": "sample",
                "features_id_col": "gene_id",
                "genome": "GRCh38",
                "differential_method": "deseq2",
                "deseq2_vs_method": "vst",
                "functional_method": "none",
                "skip_reports": False,
            }
        ),
        encoding="utf-8",
    )
    (info / "samplesheet.valid.csv").write_text(
        Path(controls["observations"]["path"]).read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    report = pipeline / "report" / "study001"
    report.mkdir(parents=True)
    (report / "study001_quarto.html").write_text("<html>report</html>\n", encoding="utf-8")
    (report / "study001.zip").write_bytes(b"synthetic-editable-report")
    processed = pipeline / "tables" / "processed_abundance" / "study001"
    processed.mkdir(parents=True)
    sample_header = "gene_id\t" + "\t".join(plan["accepted_samples"]) + "\n"
    for name in ("all.normalised_counts.tsv", "all.vst.tsv"):
        (processed / name).write_text(sample_header + "ENSG000001\t1\t2\t3\t4\n", encoding="utf-8")
    differential = pipeline / "tables" / "differential" / "study001"
    differential.mkdir(parents=True)
    result_header = "gene_id\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n"
    result_row = "ENSG000001\t10\t2\t0.5\t4\t0.001\t0.01\n"
    for name in (
        "treated_vs_control.deseq2.results.tsv",
        "treated_vs_control.deseq2.results_filtered.tsv",
        "treated_vs_control_deseq2.annotated.tsv",
    ):
        (differential / name).write_text(result_header + result_row, encoding="utf-8")
    volcano = pipeline / "plots" / "differential" / "study001" / "treated_vs_control" / "png"
    volcano.mkdir(parents=True)
    (volcano / "volcano.png").write_bytes(b"synthetic-png")
    qc = pipeline / "plots" / "qc" / "study001"
    qc.mkdir(parents=True)
    (qc / "treated_vs_control.deseq2.dispersion.png").write_bytes(b"synthetic-png")
    exploratory = pipeline / "plots" / "exploratory" / "study001" / "condition" / "png"
    exploratory.mkdir(parents=True)
    (exploratory / "pca2d.png").write_bytes(b"synthetic-png")
    (exploratory / "mad_correlation.png").write_bytes(b"synthetic-png")
    return Path(kwargs["output"]), run_root


def test_plan_subsets_accepted_inputs_and_pins_workflow(tmp_path: Path) -> None:
    kwargs = plan_kwargs(tmp_path)
    progress: list[tuple[str, int, int]] = []

    plan = create_differential_plan(
        **kwargs,
        progress=lambda label, current, total: progress.append((label, current, total)),
    )

    assert plan["stage"] == "planned_not_executed"
    assert plan["analysis"] == "differential_abundance"
    assert plan["workflow"] == {
        "name": "nf-core/differentialabundance",
        "revision": "2.0.0",
    }
    assert plan["accepted_sample_count"] == 4
    assert plan["feature_count"] == 1
    assert plan["contrast_count"] == 1
    assert "--feature_length_matrix" in plan["command_argv"]
    assert plan["command_argv"][0:5] == [
        "nextflow", "run", "nf-core/differentialabundance", "-r", "2.0.0"
    ]
    assert progress[-1][0:2] == ("gene_lengths", progress[-1][2])
    assert stat.S_IMODE(Path(kwargs["output"]).stat().st_mode) == 0o600
    with (Path(kwargs["handoff_dir"]) / "gene_counts.accepted.tsv").open() as handle:
        header = handle.readline().rstrip().split("\t")
    assert header == ["gene_id", "control_1", "control_2", "treated_1", "treated_2"]
    with (Path(kwargs["handoff_dir"]) / "contrasts.accepted.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows == [
        {
            "id": "treated_vs_control",
            "variable": "condition",
            "reference": "control",
            "target": "treated",
            "blocking": "",
        }
    ]


def test_changed_source_matrix_is_rejected(tmp_path: Path) -> None:
    kwargs = plan_kwargs(tmp_path)
    acceptance = json.loads(Path(kwargs["qc_acceptance"]).read_text(encoding="utf-8"))
    assessment = json.loads(Path(acceptance["qc_assessment"]["path"]).read_text(encoding="utf-8"))
    completion = json.loads(Path(assessment["primary_completion"]["path"]).read_text(encoding="utf-8"))
    counts = next(item for item in completion["artifacts"] if item["role"] == "gene_counts")
    Path(counts["path"]).write_text("changed\n", encoding="utf-8")

    with pytest.raises(DifferentialHandoffError, match="gene-count matrix size changed"):
        create_differential_plan(**kwargs)


def test_count_and_length_features_must_match(tmp_path: Path) -> None:
    kwargs = plan_kwargs(tmp_path)
    acceptance = json.loads(Path(kwargs["qc_acceptance"]).read_text(encoding="utf-8"))
    assessment = json.loads(Path(acceptance["qc_assessment"]["path"]).read_text(encoding="utf-8"))
    completion_path = Path(assessment["primary_completion"]["path"])
    completion = json.loads(completion_path.read_text(encoding="utf-8"))
    lengths = next(item for item in completion["artifacts"] if item["role"] == "gene_lengths")
    length_path = Path(lengths["path"])
    length_path.write_text(length_path.read_text().replace("ENSG000001", "ENSG999999"))
    lengths["size_bytes"] = length_path.stat().st_size
    from rnaseq_service.plan import sha256_file

    lengths["sha256"] = sha256_file(length_path)
    completion_path.write_text(json.dumps(completion), encoding="utf-8")
    assessment["primary_completion"]["size_bytes"] = completion_path.stat().st_size
    assessment["primary_completion"]["sha256"] = sha256_file(completion_path)
    assessment_path = Path(acceptance["qc_assessment"]["path"])
    assessment_path.write_text(json.dumps(assessment), encoding="utf-8")
    acceptance["qc_assessment"]["size_bytes"] = assessment_path.stat().st_size
    acceptance["qc_assessment"]["sha256"] = sha256_file(assessment_path)
    Path(kwargs["qc_acceptance"]).write_text(json.dumps(acceptance), encoding="utf-8")

    with pytest.raises(DifferentialHandoffError, match="different feature IDs"):
        create_differential_plan(**kwargs)


def test_plan_and_handoff_are_non_overwriting(tmp_path: Path) -> None:
    kwargs = plan_kwargs(tmp_path)
    create_differential_plan(**kwargs)

    with pytest.raises(DifferentialHandoffError, match="refusing to overwrite"):
        create_differential_plan(**kwargs)


def test_differential_cli_creates_plan_with_progress_suppressed(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    kwargs = plan_kwargs(tmp_path)

    assert differential_main(
        [
            "plan",
            "--qc-acceptance",
            str(kwargs["qc_acceptance"]),
            "--handoff-dir",
            str(kwargs["handoff_dir"]),
            "--output",
            str(kwargs["output"]),
            "--outdir",
            str(kwargs["outdir"]),
            "--workdir",
            str(kwargs["workdir"]),
            "--study-name",
            str(kwargs["study_name"]),
            "--network-mode",
            "direct",
            "--quiet",
        ]
    ) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["status"] == "planned_not_executed"
    assert summary["accepted_sample_count"] == 4
    assert summary["feature_count"] == 1


def test_completed_differential_run_is_checksum_bound(tmp_path: Path) -> None:
    plan, _ = completed_differential_fixture(tmp_path)
    output = tmp_path / "private" / "differential-completion.json"
    progress: list[tuple[str, int, int]] = []

    result = create_differential_receipt(
        run_plan=plan,
        output=output,
        progress=lambda role, current, total: progress.append((role, current, total)),
    )

    assert result["stage"] == "differential_complete_interpretation_pending"
    assert result["controls_verified"] is True
    assert result["accepted_sample_count"] == 4
    assert result["feature_count"] == 1
    assert result["contrast_count"] == 1
    assert result["contrasts"] == [
        {"id": "treated_vs_control", "result_rows": 1, "filtered_rows": 1}
    ]
    assert result["interpretation_status"] == "pending_human_review"
    assert result["automatic_biological_claims"] is False
    assert {record["role"] for record in result["artifacts"]} >= {
        "analysis_report",
        "analysis_bundle",
        "normalised_counts",
        "vst_counts",
        "deseq2_results",
        "deseq2_filtered",
        "deseq2_annotated",
        "volcano_plot",
        "pca_plot",
        "mad_plot",
    }
    assert progress
    assert stat.S_IMODE(output.stat().st_mode) == 0o600
    with pytest.raises(DifferentialHandoffError, match="refusing to overwrite"):
        create_differential_receipt(run_plan=plan, output=output)


def test_failed_differential_trace_is_not_completion(tmp_path: Path) -> None:
    plan, run_root = completed_differential_fixture(tmp_path)
    trace = run_root / "execution" / "trace.tsv"
    trace.write_text(trace.read_text().replace("COMPLETED", "FAILED"), encoding="utf-8")

    with pytest.raises(DifferentialHandoffError, match="not complete and successful"):
        create_differential_receipt(run_plan=plan, output=tmp_path / "completion.json")


def test_missing_contrast_result_is_rejected(tmp_path: Path) -> None:
    plan, run_root = completed_differential_fixture(tmp_path)
    result = next(
        (run_root / "differentialabundance" / "tables" / "differential").rglob(
            "treated_vs_control.deseq2.results.tsv"
        )
    )
    result.unlink()

    with pytest.raises(DifferentialHandoffError, match="DESeq2 results for treated_vs_control"):
        create_differential_receipt(run_plan=plan, output=tmp_path / "completion.json")


def test_completion_cli_reports_sealed_stage(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    plan, _ = completed_differential_fixture(tmp_path)
    output = tmp_path / "private" / "differential-completion-cli.json"

    assert differential_main(
        ["complete", "--run-plan", str(plan), "--output", str(output), "--quiet"]
    ) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["status"] == "differential_complete_interpretation_pending"
    assert summary["contrast_count"] == 1
