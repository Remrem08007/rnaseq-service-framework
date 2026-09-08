from __future__ import annotations

import csv
import json
import stat
from pathlib import Path

import pytest

from rnaseq_service.primary import create_primary_receipt
from rnaseq_service.qc import QCError, evaluate_qc
from test_primary import completed_fixture


def write_policy(path: Path, *, value_column: str = "mapped", required: bool = True) -> Path:
    path.write_text(
        "schema_version = 1\n"
        "[study]\n"
        "robust_z_threshold = 3.0\n"
        "min_samples_for_outliers = 4\n"
        "[[metrics]]\n"
        'id = "mapping_percent"\n'
        'source_file = "multiqc_general_stats.txt"\n'
        'sample_column = "Sample"\n'
        f'value_column = "{value_column}"\n'
        f"required = {str(required).lower()}\n"
        "lower = 70.0\n"
        "upper = 100.0\n"
        "outlier = true\n",
        encoding="utf-8",
    )
    return path


def qc_fixture(tmp_path: Path) -> tuple[Path, Path]:
    plan, run_root = completed_fixture(tmp_path)
    table = (
        run_root
        / "rnaseq"
        / "multiqc"
        / "star_salmon"
        / "multiqc_report_data"
        / "multiqc_general_stats.txt"
    )
    with table.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Sample", "mapped"])
        writer.writerows(
            [
                ["control_1", "96"],
                ["control_2", "95"],
                ["treated_1", "94"],
                ["treated_2", "40"],
            ]
        )
    completion = tmp_path / "private" / "completion.json"
    create_primary_receipt(run_plan=plan, output=completion)
    return completion, write_policy(tmp_path / "policy.toml")


def test_qc_assessment_flags_but_never_excludes(tmp_path: Path) -> None:
    completion, policy = qc_fixture(tmp_path)
    outdir = tmp_path / "qc"
    progress: list[tuple[str, int, int]] = []

    result = evaluate_qc(
        completion_receipt=completion,
        policy_path=policy,
        outdir=outdir,
        progress=lambda metric, current, total: progress.append((metric, current, total)),
    )

    assert result["status"] == "review_required"
    assert result["flagged_samples"] == ["treated_2"]
    flags = result["samples"][3]["flags"]
    assert {flag["kind"] for flag in flags} == {"below_lower", "robust_outlier"}
    assert result["automatic_exclusions"] == 0
    assert progress[-1] == ("mapping_percent", 1, 1)
    assert stat.S_IMODE((outdir / "qc_assessment.json").stat().st_mode) == 0o600
    with (outdir / "qc_decisions.tsv").open(encoding="utf-8") as handle:
        decisions = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["decision"] for row in decisions} == {"review"}
    with pytest.raises(QCError, match="refusing to overwrite"):
        evaluate_qc(completion_receipt=completion, policy_path=policy, outdir=outdir)


def test_missing_required_column_writes_stop_state(tmp_path: Path) -> None:
    completion, _ = qc_fixture(tmp_path)
    policy = write_policy(tmp_path / "missing.toml", value_column="not_present")

    result = evaluate_qc(
        completion_receipt=completion,
        policy_path=policy,
        outdir=tmp_path / "qc-stop",
    )

    assert result["status"] == "stop_missing_evidence"
    assert result["missing_required_sample_count"] == 4
    assert result["evidence_issues"][0]["required"] is True


def test_changed_multiqc_evidence_is_rejected(tmp_path: Path) -> None:
    completion, policy = qc_fixture(tmp_path)
    receipt = json.loads(completion.read_text(encoding="utf-8"))
    record = next(
        item for item in receipt["artifacts"] if item["role"] == "multiqc_general_stats"
    )
    Path(record["path"]).write_text("changed\n", encoding="utf-8")

    with pytest.raises(QCError, match="size changed|checksum changed"):
        evaluate_qc(
            completion_receipt=completion,
            policy_path=policy,
            outdir=tmp_path / "qc",
        )


def test_source_path_cannot_escape_multiqc_directory(tmp_path: Path) -> None:
    completion, _ = qc_fixture(tmp_path)
    policy = write_policy(tmp_path / "escape.toml")
    text = policy.read_text(encoding="utf-8").replace(
        'source_file = "multiqc_general_stats.txt"', 'source_file = "../secret.tsv"'
    )
    policy.write_text(text, encoding="utf-8")

    with pytest.raises(QCError, match="escapes the MultiQC directory"):
        evaluate_qc(
            completion_receipt=completion,
            policy_path=policy,
            outdir=tmp_path / "qc",
        )


def test_categorical_strand_disagreement_is_flagged(tmp_path: Path) -> None:
    completion, _ = qc_fixture(tmp_path)
    receipt = json.loads(completion.read_text(encoding="utf-8"))
    root = Path(receipt["multiqc"]["data_directory"])
    (root / "strand.tsv").write_text(
        "Sample\tstatus\n"
        "control_1\tpass\ncontrol_2\tpass\ntreated_1\tpass\ntreated_2\tfail\n",
        encoding="utf-8",
    )
    run_plan = Path(receipt["run_plan"]["path"])
    completion.unlink()
    create_primary_receipt(run_plan=run_plan, output=completion)
    policy = tmp_path / "strand-policy.toml"
    policy.write_text(
        "schema_version = 1\n"
        "[[metrics]]\n"
        'id = "strand_agreement"\n'
        'kind = "categorical"\n'
        'source_file = "strand.tsv"\n'
        'sample_column = "Sample"\n'
        'value_column = "status"\n'
        'allowed_values = ["pass"]\n'
        "required = true\n"
        "outlier = false\n",
        encoding="utf-8",
    )

    result = evaluate_qc(
        completion_receipt=completion,
        policy_path=policy,
        outdir=tmp_path / "strand-qc",
    )

    assert result["flagged_samples"] == ["treated_2"]
    assert result["samples"][3]["flags"][0]["kind"] == "unexpected_category"
