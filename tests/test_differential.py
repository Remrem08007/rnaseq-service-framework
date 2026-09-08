from __future__ import annotations

import csv
import json
import stat
from pathlib import Path

import pytest

from rnaseq_service.differential import DifferentialHandoffError, create_differential_plan
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
