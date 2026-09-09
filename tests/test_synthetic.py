from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import pytest

from rnaseq_service.preflight import run_preflight
from rnaseq_service.synthetic import PAIR_COUNTS, SyntheticStudyError, create_synthetic_study
from rnaseq_service.synthetic_cli import main as synthetic_main


def tree_hashes(root: Path) -> dict[str, str]:
    return {
        str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def test_synthetic_study_is_deterministic_valid_and_directional(tmp_path: Path) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    progress: list[tuple[str, int, int]] = []

    truth = create_synthetic_study(
        output_dir=first,
        progress=lambda sample, current, total: progress.append((sample, current, total)),
    )
    create_synthetic_study(output_dir=second)

    assert tree_hashes(first) == tree_hashes(second)
    assert truth["dataset"] == "deterministic_synthetic_not_client_data"
    assert truth["contrast"]["expected_positive_log2_fold_change"] == [
        "gene_treated_marker"
    ]
    assert truth["contrast"]["expected_negative_log2_fold_change"] == [
        "gene_control_marker"
    ]
    assert progress[-1] == ("treated_3", 6, 6)
    report = run_preflight(
        first / "intake" / "samplesheet.csv",
        first / "intake" / "design.csv",
        first / "intake" / "contrasts.csv",
        check_files=True,
    )
    assert report.valid, report.errors
    assert report.summary["n_samples"] == 6
    assert report.summary["n_fastq_files"] == 12


def test_fastqs_have_matching_pairs_and_truth_counts(tmp_path: Path) -> None:
    root = tmp_path / "study"
    create_synthetic_study(output_dir=root)

    for sample, counts in PAIR_COUNTS.items():
        with gzip.open(root / "fastq" / f"{sample}_R1.fastq.gz", "rt") as handle:
            read_1 = handle.read().splitlines()
        with gzip.open(root / "fastq" / f"{sample}_R2.fastq.gz", "rt") as handle:
            read_2 = handle.read().splitlines()
        assert len(read_1) == len(read_2) == sum(counts) * 4
        assert all(len(read_1[index]) == 50 for index in range(1, len(read_1), 4))
        assert all(len(read_2[index]) == 50 for index in range(1, len(read_2), 4))
        assert [line.removesuffix("/1") for line in read_1[::4]] == [
            line.removesuffix("/2") for line in read_2[::4]
        ]


def test_synthetic_study_is_non_overwriting(tmp_path: Path) -> None:
    root = tmp_path / "study"
    create_synthetic_study(output_dir=root)

    with pytest.raises(SyntheticStudyError, match="refusing to overwrite"):
        create_synthetic_study(output_dir=root)


def test_synthetic_cli(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    root = tmp_path / "study"

    assert synthetic_main(["--output-dir", str(root), "--quiet"]) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["n_samples"] == 6
    assert summary["n_fastq_files"] == 12
    assert Path(summary["truth"]) == root / "truth.json"
