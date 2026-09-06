from __future__ import annotations

import csv
import hashlib
import json
import stat
from pathlib import Path

import pytest

from rnaseq_service.inputs import (
    InputManifestError,
    create_input_manifest,
    verify_input_manifest,
)


def samplesheet_fixture(tmp_path: Path) -> tuple[Path, list[Path]]:
    reads = tmp_path / "reads"
    reads.mkdir()
    paths = [reads / "sample_R1.fastq.gz", reads / "sample_R2.fastq.gz"]
    paths[0].write_bytes(b"ACGT\n")
    paths[1].write_bytes(b"TGCA\n")
    samplesheet = tmp_path / "samplesheet.csv"
    with samplesheet.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample", "fastq_1", "fastq_2", "strandedness"])
        writer.writerow(["sample", str(paths[0]), str(paths[1]), "auto"])
    return samplesheet, paths


def test_manifest_hashes_fastqs_with_progress_and_private_mode(tmp_path: Path) -> None:
    samplesheet, paths = samplesheet_fixture(tmp_path)
    output = tmp_path / "private" / "inputs.json"
    progress: list[tuple[str, int, int]] = []

    manifest = create_input_manifest(
        samplesheet=samplesheet,
        output=output,
        progress=lambda name, current, total: progress.append((name, current, total)),
    )

    assert manifest["n_fastq_files"] == 2
    assert manifest["total_bytes"] == 10
    assert manifest["files"][0]["sha256"] == hashlib.sha256(b"ACGT\n").hexdigest()
    assert progress[-1][1:] == (10, 10)
    assert stat.S_IMODE(output.stat().st_mode) == 0o600
    assert verify_input_manifest(output, expected_samplesheet=samplesheet) == manifest
    with pytest.raises(InputManifestError, match="refusing to overwrite"):
        create_input_manifest(samplesheet=samplesheet, output=output)


def test_manifest_rejects_changed_fastq_metadata(tmp_path: Path) -> None:
    samplesheet, paths = samplesheet_fixture(tmp_path)
    output = tmp_path / "inputs.json"
    create_input_manifest(samplesheet=samplesheet, output=output)
    paths[0].write_bytes(b"changed\n")

    with pytest.raises(InputManifestError, match="size changed"):
        verify_input_manifest(output, expected_samplesheet=samplesheet)


def test_manifest_can_rehash_fastqs(tmp_path: Path) -> None:
    samplesheet, paths = samplesheet_fixture(tmp_path)
    output = tmp_path / "inputs.json"
    create_input_manifest(samplesheet=samplesheet, output=output)
    payload = json.loads(output.read_text(encoding="utf-8"))
    payload["files"][0]["sha256"] = "0" * 64
    output.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(InputManifestError, match="checksum changed"):
        verify_input_manifest(output, expected_samplesheet=samplesheet, verify_hashes=True)


def test_manifest_must_match_expected_samplesheet(tmp_path: Path) -> None:
    samplesheet, _ = samplesheet_fixture(tmp_path)
    output = tmp_path / "inputs.json"
    create_input_manifest(samplesheet=samplesheet, output=output)
    other = tmp_path / "other.csv"
    other.write_text(samplesheet.read_text(encoding="utf-8"), encoding="utf-8")

    with pytest.raises(InputManifestError, match="different samplesheet"):
        verify_input_manifest(output, expected_samplesheet=other)
