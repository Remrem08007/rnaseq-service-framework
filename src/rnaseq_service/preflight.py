"""Validate RNA-seq intake artifacts before launching external workflows."""

from __future__ import annotations

import csv
import json
import os
import re
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable


SAMPLE_COLUMNS = ("sample", "fastq_1", "fastq_2", "strandedness")
CONTRAST_COLUMNS = ("contrast_id", "variable", "reference", "target")
SAMPLE_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
STRANDEDNESS_VALUES = {"auto", "forward", "reverse", "unstranded"}


@dataclass
class PreflightReport:
    """Machine-readable preflight outcome."""

    valid: bool = True
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    summary: dict[str, object] = field(default_factory=dict)

    def error(self, message: str) -> None:
        self.errors.append(message)
        self.valid = False

    def warn(self, message: str) -> None:
        self.warnings.append(message)

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2, sort_keys=True) + "\n"


def _read_csv(path: Path, report: PreflightReport) -> tuple[list[str], list[dict[str, str]]]:
    try:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                report.error(f"{path}: missing CSV header")
                return [], []
            headers = [name.strip() for name in reader.fieldnames]
            if any(not name for name in headers):
                report.error(f"{path}: header contains an empty column name")
            if len(headers) != len(set(headers)):
                report.error(f"{path}: header contains duplicate column names")
            rows: list[dict[str, str]] = []
            for line_number, raw in enumerate(reader, start=2):
                if None in raw:
                    report.error(f"{path}:{line_number}: row has more values than the header")
                    continue
                row = {str(key).strip(): (value or "").strip() for key, value in raw.items()}
                if not any(row.values()):
                    report.warn(f"{path}:{line_number}: ignored empty row")
                    continue
                row["__line__"] = str(line_number)
                rows.append(row)
            return headers, rows
    except FileNotFoundError:
        report.error(f"{path}: file does not exist")
    except (OSError, UnicodeError, csv.Error) as exc:
        report.error(f"{path}: could not read CSV: {exc}")
    return [], []


def _missing_columns(headers: Iterable[str], required: Iterable[str]) -> list[str]:
    available = set(headers)
    return [column for column in required if column not in available]


def _looks_like_fastq(path: str) -> bool:
    return path.lower().endswith(FASTQ_SUFFIXES)


def _resolve_input_path(value: str, csv_path: Path) -> Path:
    expanded = Path(os.path.expandvars(os.path.expanduser(value)))
    if expanded.is_absolute():
        return expanded
    return csv_path.parent / expanded


def validate_samplesheet(
    path: Path,
    report: PreflightReport,
    *,
    check_files: bool = False,
) -> set[str]:
    """Validate the nf-core/rnaseq four-column samplesheet contract."""

    headers, rows = _read_csv(path, report)
    missing = _missing_columns(headers, SAMPLE_COLUMNS)
    if missing:
        report.error(f"{path}: missing required columns: {', '.join(missing)}")
        return set()
    if not rows:
        report.error(f"{path}: samplesheet has no data rows")
        return set()

    samples: set[str] = set()
    layouts: dict[str, set[str]] = defaultdict(set)
    strandedness_by_sample: dict[str, set[str]] = defaultdict(set)
    seen_rows: set[tuple[str, str, str, str]] = set()
    fastq_owners: dict[str, tuple[str, int]] = {}

    for row in rows:
        line = int(row["__line__"])
        sample = row["sample"]
        fastq_1 = row["fastq_1"]
        fastq_2 = row["fastq_2"]
        strandedness = row["strandedness"].lower()

        if not SAMPLE_ID_PATTERN.fullmatch(sample):
            report.error(
                f"{path}:{line}: invalid sample identifier {sample!r}; use letters, "
                "numbers, '.', '_' or '-' and do not start with punctuation"
            )
        else:
            samples.add(sample)

        if not fastq_1:
            report.error(f"{path}:{line}: fastq_1 is required")
        elif not _looks_like_fastq(fastq_1):
            report.error(f"{path}:{line}: fastq_1 is not a FASTQ path: {fastq_1!r}")

        if fastq_2 and not _looks_like_fastq(fastq_2):
            report.error(f"{path}:{line}: fastq_2 is not a FASTQ path: {fastq_2!r}")
        if fastq_1 and fastq_1 == fastq_2:
            report.error(f"{path}:{line}: fastq_1 and fastq_2 refer to the same path")

        if strandedness not in STRANDEDNESS_VALUES:
            allowed = ", ".join(sorted(STRANDEDNESS_VALUES))
            report.error(
                f"{path}:{line}: invalid strandedness {row['strandedness']!r}; "
                f"expected one of {allowed}"
            )

        layout = "paired" if fastq_2 else "single"
        layouts[sample].add(layout)
        strandedness_by_sample[sample].add(strandedness)
        signature = (sample, fastq_1, fastq_2, strandedness)
        if signature in seen_rows:
            report.error(f"{path}:{line}: duplicate samplesheet row")
        seen_rows.add(signature)

        for fastq in (fastq_1, fastq_2):
            if not fastq:
                continue
            previous = fastq_owners.get(fastq)
            if previous is not None:
                previous_sample, previous_line = previous
                report.error(
                    f"{path}:{line}: FASTQ {fastq!r} was already assigned to "
                    f"sample {previous_sample!r} on line {previous_line}"
                )
            else:
                fastq_owners[fastq] = (sample, line)
            if check_files:
                resolved = _resolve_input_path(fastq, path)
                if not resolved.is_file():
                    report.error(f"{path}:{line}: FASTQ does not exist: {resolved}")

    for sample, sample_layouts in sorted(layouts.items()):
        if len(sample_layouts) > 1:
            report.error(
                f"{path}: sample {sample!r} mixes single-end and paired-end rows"
            )
    for sample, values in sorted(strandedness_by_sample.items()):
        if len(values) > 1:
            rendered = ", ".join(sorted(values))
            report.error(
                f"{path}: sample {sample!r} has inconsistent strandedness "
                f"across rows: {rendered}"
            )

    report.summary.update(
        {
            "n_samples": len(samples),
            "n_samplesheet_rows": len(rows),
            "n_fastq_files": len(fastq_owners),
        }
    )
    return samples


def validate_design(path: Path, samples: set[str], report: PreflightReport) -> dict[str, dict[str, str]]:
    """Validate one-row-per-sample analysis metadata."""

    headers, rows = _read_csv(path, report)
    if "sample" not in headers:
        report.error(f"{path}: missing required column: sample")
        return {}
    variables = [header for header in headers if header != "sample"]
    if not variables:
        report.error(f"{path}: design must contain at least one variable column")
        return {}

    design: dict[str, dict[str, str]] = {}
    for row in rows:
        line = int(row["__line__"])
        sample = row["sample"]
        if not sample:
            report.error(f"{path}:{line}: sample is required")
            continue
        if sample in design:
            report.error(f"{path}:{line}: duplicate design row for sample {sample!r}")
            continue
        values = {variable: row.get(variable, "") for variable in variables}
        empty = [variable for variable, value in values.items() if not value]
        if empty:
            report.error(
                f"{path}:{line}: empty design values for: {', '.join(empty)}"
            )
        design[sample] = values

    design_samples = set(design)
    missing_from_design = sorted(samples - design_samples)
    unknown_in_design = sorted(design_samples - samples)
    if missing_from_design:
        report.error(
            f"{path}: samples missing from design: {', '.join(missing_from_design)}"
        )
    if unknown_in_design:
        report.error(
            f"{path}: design contains unknown samples: {', '.join(unknown_in_design)}"
        )

    report.summary["design_variables"] = variables
    return design


def validate_contrasts(
    path: Path,
    design: dict[str, dict[str, str]],
    report: PreflightReport,
    *,
    min_replicates: int = 2,
) -> None:
    """Validate requested categorical contrasts against the design table."""

    headers, rows = _read_csv(path, report)
    missing = _missing_columns(headers, CONTRAST_COLUMNS)
    if missing:
        report.error(f"{path}: missing required columns: {', '.join(missing)}")
        return
    if not rows:
        report.error(f"{path}: contrasts table has no data rows")
        return

    variables = set(next(iter(design.values()), {}))
    ids: set[str] = set()
    for row in rows:
        line = int(row["__line__"])
        contrast_id = row["contrast_id"]
        variable = row["variable"]
        reference = row["reference"]
        target = row["target"]

        if not SAMPLE_ID_PATTERN.fullmatch(contrast_id):
            report.error(f"{path}:{line}: invalid contrast_id {contrast_id!r}")
        elif contrast_id in ids:
            report.error(f"{path}:{line}: duplicate contrast_id {contrast_id!r}")
        ids.add(contrast_id)

        if variable not in variables:
            report.error(f"{path}:{line}: unknown design variable {variable!r}")
            continue
        if not reference or not target:
            report.error(f"{path}:{line}: reference and target are required")
            continue
        if reference == target:
            report.error(f"{path}:{line}: reference and target must differ")
            continue

        counts = Counter(values[variable] for values in design.values())
        for role, level in (("reference", reference), ("target", target)):
            if level not in counts:
                report.error(
                    f"{path}:{line}: {role} level {level!r} is absent from "
                    f"design variable {variable!r}"
                )
            elif counts[level] < min_replicates:
                report.error(
                    f"{path}:{line}: {role} level {level!r} has {counts[level]} "
                    f"replicate(s); minimum is {min_replicates}"
                )

    report.summary["n_contrasts"] = len(rows)


def run_preflight(
    samplesheet: Path,
    design: Path,
    contrasts: Path,
    *,
    check_files: bool = False,
    min_replicates: int = 2,
) -> PreflightReport:
    """Validate all intake artifacts and return one combined report."""

    report = PreflightReport()
    samples = validate_samplesheet(samplesheet, report, check_files=check_files)
    design_rows = validate_design(design, samples, report)
    validate_contrasts(
        contrasts,
        design_rows,
        report,
        min_replicates=min_replicates,
    )
    return report
