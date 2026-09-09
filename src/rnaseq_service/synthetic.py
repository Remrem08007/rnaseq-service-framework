"""Generate a deterministic, non-client paired-end RNA-seq validation study."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import random
import shutil
from pathlib import Path
from typing import Callable


SyntheticProgress = Callable[[str, int, int], None]
GENES = (
    "gene_housekeeping",
    "gene_control_marker",
    "gene_treated_marker",
    "gene_neutral",
)
PAIR_COUNTS = {
    "control_1": (12, 18, 3, 7),
    "control_2": (12, 17, 4, 7),
    "control_3": (12, 19, 3, 7),
    "treated_1": (12, 3, 18, 7),
    "treated_2": (12, 4, 17, 7),
    "treated_3": (12, 3, 19, 7),
}
DNA = "ACGT"


class SyntheticStudyError(ValueError):
    """Raised when a synthetic study cannot be created safely."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _write_gzip(path: Path, text: str) -> None:
    with path.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            compressed.write(text.encode("ascii"))


def _transcripts(seed: int) -> dict[str, str]:
    rng = random.Random(seed)
    return {
        gene: "".join(rng.choice(DNA) for _ in range(500))
        for gene in GENES
    }


def _fastq_pair(
    sample: str,
    counts: tuple[int, int, int, int],
    transcripts: dict[str, str],
    seed: int,
) -> tuple[str, str]:
    read_length = 50
    fragment_length = 140
    rng = random.Random(f"{seed}:{sample}")
    read_1: list[str] = []
    read_2: list[str] = []
    for gene, count in zip(GENES, counts, strict=True):
        transcript = transcripts[gene]
        for pair_index in range(1, count + 1):
            start = rng.randrange(0, len(transcript) - fragment_length + 1)
            r1 = transcript[start : start + read_length]
            r2_source = transcript[
                start + fragment_length - read_length : start + fragment_length
            ]
            name = f"synthetic:{sample}:{gene}:{pair_index}:{start}"
            read_1.extend((f"@{name}/1", r1, "+", "I" * read_length))
            read_2.extend(
                (f"@{name}/2", _reverse_complement(r2_source), "+", "I" * read_length)
            )
    return "\n".join(read_1) + "\n", "\n".join(read_2) + "\n"


def create_synthetic_study(
    *,
    output_dir: Path,
    seed: int = 20260909,
    progress: SyntheticProgress | None = None,
) -> dict[str, object]:
    """Create an immutable tiny study, mini-reference, and expected truth."""

    if output_dir.exists():
        raise SyntheticStudyError(f"refusing to overwrite synthetic study: {output_dir}")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise SyntheticStudyError("seed must be a non-negative integer")
    root = output_dir.resolve()
    try:
        fastq_dir = root / "fastq"
        intake_dir = root / "intake"
        reference_dir = root / "reference"
        fastq_dir.mkdir(parents=True)
        intake_dir.mkdir()
        reference_dir.mkdir()
        transcripts = _transcripts(seed)

        fasta = reference_dir / "synthetic.fa"
        gtf = reference_dir / "synthetic.gtf"
        with fasta.open("x", encoding="ascii") as handle:
            for index, gene in enumerate(GENES, start=1):
                handle.write(f">chrSynthetic{index}\n{transcripts[gene]}\n")
        with gtf.open("x", encoding="ascii") as handle:
            for index, gene in enumerate(GENES, start=1):
                attributes = f'gene_id "{gene}"; transcript_id "{gene}.t1"; gene_name "{gene}";'
                handle.write(
                    f"chrSynthetic{index}\tsynthetic\texon\t1\t500\t.\t+\t.\t{attributes}\n"
                )

        sample_rows: list[list[str]] = []
        design_rows: list[list[str]] = []
        fastq_records: list[dict[str, object]] = []
        total_samples = len(PAIR_COUNTS)
        for completed, (sample, counts) in enumerate(PAIR_COUNTS.items(), start=1):
            r1_text, r2_text = _fastq_pair(sample, counts, transcripts, seed)
            r1 = fastq_dir / f"{sample}_R1.fastq.gz"
            r2 = fastq_dir / f"{sample}_R2.fastq.gz"
            _write_gzip(r1, r1_text)
            _write_gzip(r2, r2_text)
            sample_rows.append(
                [sample, f"../fastq/{r1.name}", f"../fastq/{r2.name}", "unstranded"]
            )
            condition = sample.split("_", 1)[0]
            design_rows.append([sample, condition])
            fastq_records.extend(
                {
                    "sample": sample,
                    "read": read,
                    "relative_path": str(path.relative_to(root)),
                    "size_bytes": path.stat().st_size,
                    "sha256": _sha256(path),
                }
                for read, path in (("R1", r1), ("R2", r2))
            )
            if progress is not None:
                progress(sample, completed, total_samples)

        def write_csv(path: Path, header: list[str], rows: list[list[str]]) -> None:
            with path.open("x", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, lineterminator="\n")
                writer.writerow(header)
                writer.writerows(rows)

        samplesheet = intake_dir / "samplesheet.csv"
        design = intake_dir / "design.csv"
        contrasts = intake_dir / "contrasts.csv"
        write_csv(samplesheet, ["sample", "fastq_1", "fastq_2", "strandedness"], sample_rows)
        write_csv(design, ["sample", "condition"], design_rows)
        write_csv(
            contrasts,
            ["contrast_id", "variable", "reference", "target"],
            [["treated_vs_control", "condition", "control", "treated"]],
        )
        truth: dict[str, object] = {
            "schema_version": 1,
            "dataset": "deterministic_synthetic_not_client_data",
            "seed": seed,
            "paired_end": True,
            "read_length": 50,
            "fragment_length": 140,
            "samples": list(PAIR_COUNTS),
            "pair_counts": {
                sample: dict(zip(GENES, counts, strict=True))
                for sample, counts in PAIR_COUNTS.items()
            },
            "contrast": {
                "id": "treated_vs_control",
                "reference": "control",
                "target": "treated",
                "expected_positive_log2_fold_change": ["gene_treated_marker"],
                "expected_negative_log2_fold_change": ["gene_control_marker"],
                "expected_near_zero_log2_fold_change": ["gene_housekeeping", "gene_neutral"],
            },
            "files": sorted(fastq_records, key=lambda record: str(record["relative_path"])),
            "reference": {
                "fasta": str(fasta.relative_to(root)),
                "fasta_sha256": _sha256(fasta),
                "gtf": str(gtf.relative_to(root)),
                "gtf_sha256": _sha256(gtf),
            },
        }
        truth_path = root / "truth.json"
        with truth_path.open("x", encoding="utf-8") as handle:
            json.dump(truth, handle, indent=2, sort_keys=True)
            handle.write("\n")
        return truth
    except Exception:
        shutil.rmtree(root, ignore_errors=True)
        raise

