from __future__ import annotations

import csv
import json
from pathlib import Path

from rnaseq_service.cli import main
from rnaseq_service.preflight import PreflightReport, run_preflight, validate_samplesheet


FIXTURES = Path(__file__).parent / "fixtures"


def write_csv(path: Path, headers: list[str], rows: list[list[str]]) -> Path:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        writer.writerows(rows)
    return path


def test_valid_fixture_passes() -> None:
    report = run_preflight(
        FIXTURES / "samplesheet.csv",
        FIXTURES / "design.csv",
        FIXTURES / "contrasts.csv",
    )

    assert report.valid
    assert report.errors == []
    assert report.summary == {
        "design_variables": ["condition"],
        "n_contrasts": 1,
        "n_fastq_files": 8,
        "n_samples": 4,
        "n_samplesheet_rows": 4,
    }


def test_missing_samplesheet_column_fails(tmp_path: Path) -> None:
    samplesheet = write_csv(
        tmp_path / "samples.csv",
        ["sample", "fastq_1", "strandedness"],
        [["sample_1", "sample_1.fastq.gz", "auto"]],
    )
    report = PreflightReport()

    validate_samplesheet(samplesheet, report)

    assert not report.valid
    assert any("missing required columns: fastq_2" in error for error in report.errors)


def test_invalid_strandedness_and_sample_id_fail(tmp_path: Path) -> None:
    samplesheet = write_csv(
        tmp_path / "samples.csv",
        ["sample", "fastq_1", "fastq_2", "strandedness"],
        [["bad sample", "r1.fastq.gz", "r2.fastq.gz", "sometimes"]],
    )
    report = PreflightReport()

    validate_samplesheet(samplesheet, report)

    assert not report.valid
    assert any("invalid sample identifier" in error for error in report.errors)
    assert any("invalid strandedness" in error for error in report.errors)


def test_fastq_cannot_be_reused(tmp_path: Path) -> None:
    samplesheet = write_csv(
        tmp_path / "samples.csv",
        ["sample", "fastq_1", "fastq_2", "strandedness"],
        [
            ["sample_1", "shared.fastq.gz", "", "unstranded"],
            ["sample_2", "shared.fastq.gz", "", "unstranded"],
        ],
    )
    report = PreflightReport()

    validate_samplesheet(samplesheet, report)

    assert not report.valid
    assert any("already assigned" in error for error in report.errors)


def test_sample_cannot_mix_read_layouts(tmp_path: Path) -> None:
    samplesheet = write_csv(
        tmp_path / "samples.csv",
        ["sample", "fastq_1", "fastq_2", "strandedness"],
        [
            ["sample_1", "lane1_R1.fastq.gz", "lane1_R2.fastq.gz", "auto"],
            ["sample_1", "lane2_R1.fastq.gz", "", "auto"],
        ],
    )
    report = PreflightReport()

    validate_samplesheet(samplesheet, report)

    assert not report.valid
    assert any("mixes single-end and paired-end" in error for error in report.errors)


def test_sample_lanes_must_agree_on_strandedness(tmp_path: Path) -> None:
    samplesheet = write_csv(
        tmp_path / "samples.csv",
        ["sample", "fastq_1", "fastq_2", "strandedness"],
        [
            ["sample_1", "lane1_R1.fastq.gz", "lane1_R2.fastq.gz", "forward"],
            ["sample_1", "lane2_R1.fastq.gz", "lane2_R2.fastq.gz", "reverse"],
        ],
    )
    report = PreflightReport()

    validate_samplesheet(samplesheet, report)

    assert not report.valid
    assert any("inconsistent strandedness" in error for error in report.errors)


def test_design_must_match_samples(tmp_path: Path) -> None:
    design = write_csv(
        tmp_path / "design.csv",
        ["sample", "condition"],
        [["sample_1", "control"], ["unknown", "treated"]],
    )
    report = run_preflight(
        FIXTURES / "samplesheet.csv",
        design,
        FIXTURES / "contrasts.csv",
    )

    assert not report.valid
    assert any("samples missing from design" in error for error in report.errors)
    assert any("design contains unknown samples" in error for error in report.errors)


def test_contrast_levels_and_replicates_are_checked(tmp_path: Path) -> None:
    contrasts = write_csv(
        tmp_path / "contrasts.csv",
        ["contrast_id", "variable", "reference", "target"],
        [
            ["missing_level", "condition", "control", "absent"],
            ["same_level", "condition", "control", "control"],
        ],
    )
    report = run_preflight(
        FIXTURES / "samplesheet.csv",
        FIXTURES / "design.csv",
        contrasts,
    )

    assert not report.valid
    assert any("target level 'absent' is absent" in error for error in report.errors)
    assert any("reference and target must differ" in error for error in report.errors)


def test_check_files_resolves_relative_to_samplesheet(tmp_path: Path) -> None:
    reads = tmp_path / "reads"
    reads.mkdir()
    (reads / "sample.fastq.gz").touch()
    samplesheet = write_csv(
        tmp_path / "samples.csv",
        ["sample", "fastq_1", "fastq_2", "strandedness"],
        [["sample", "reads/sample.fastq.gz", "", "unstranded"]],
    )
    report = PreflightReport()

    validate_samplesheet(samplesheet, report, check_files=True)

    assert report.valid


def test_cli_writes_json_and_returns_two_on_failure(tmp_path: Path) -> None:
    output = tmp_path / "report.json"
    exit_code = main(
        [
            "--samplesheet",
            str(FIXTURES / "samplesheet.csv"),
            "--design",
            str(FIXTURES / "design.csv"),
            "--contrasts",
            str(tmp_path / "missing.csv"),
            "--json-out",
            str(output),
        ]
    )

    assert exit_code == 2
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["valid"] is False
    assert payload["errors"]
