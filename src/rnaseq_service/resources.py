"""Summarize Nextflow task resource observations from a trace TSV."""

from __future__ import annotations

import csv
import json
import math
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Callable, Mapping


ProgressCallback = Callable[[int, int], None]
StorageProgressCallback = Callable[[str, int], None]
SIZE = re.compile(r"^([0-9]+(?:\.[0-9]+)?)\s*([KMGTPE]?i?B)$", re.IGNORECASE)
DURATION_PART = re.compile(r"([0-9]+(?:\.[0-9]+)?)\s*(ms|us|µs|ns|[smhd])", re.IGNORECASE)


class ResourceSummaryError(ValueError):
    """Raised when a trace cannot produce a trustworthy resource summary."""


def parse_size(value: str) -> int | None:
    text = value.strip()
    if text in {"", "-"}:
        return None
    match = SIZE.fullmatch(text)
    if match is None:
        raise ResourceSummaryError(f"unsupported memory value: {value!r}")
    number = float(match.group(1))
    unit = match.group(2).upper().replace("IB", "B")
    powers = {"B": 0, "KB": 1, "MB": 2, "GB": 3, "TB": 4, "PB": 5, "EB": 6}
    return round(number * (1024 ** powers[unit]))


def parse_duration(value: str) -> float | None:
    text = value.strip()
    if text in {"", "-"}:
        return None
    scales = {
        "ns": 1e-9,
        "us": 1e-6,
        "µs": 1e-6,
        "ms": 1e-3,
        "s": 1.0,
        "m": 60.0,
        "h": 3600.0,
        "d": 86400.0,
    }
    matches = list(DURATION_PART.finditer(text))
    if not matches or "".join(match.group(0) for match in matches).replace(" ", "") != text.replace(" ", ""):
        raise ResourceSummaryError(f"unsupported duration value: {value!r}")
    return sum(float(match.group(1)) * scales[match.group(2).lower()] for match in matches)


def parse_cpu(value: str) -> float | None:
    text = value.strip()
    if text in {"", "-"}:
        return None
    if text.endswith("%"):
        text = text[:-1]
    try:
        number = float(text)
    except ValueError as exc:
        raise ResourceSummaryError(f"unsupported CPU value: {value!r}") from exc
    if not math.isfinite(number) or number < 0:
        raise ResourceSummaryError(f"unsupported CPU value: {value!r}")
    return number


def _process_name(row: dict[str, str]) -> str:
    explicit = row.get("process", "").strip()
    if explicit:
        return explicit
    name = row["name"].strip()
    return re.sub(r"\s+\([^()]*(?:\([^()]*\)[^()]*)*\)$", "", name) or name


def summarize_trace(
    trace: Path,
    *,
    progress: ProgressCallback | None = None,
) -> dict[str, object]:
    """Aggregate task states and observed resources by process name."""

    resolved = trace.resolve(strict=True)
    with resolved.open("r", encoding="utf-8", newline="") as handle:
        total_rows = max(0, sum(1 for _ in handle) - 1)
    required = {"name", "status", "realtime", "%cpu", "peak_rss"}
    groups: dict[str, dict[str, object]] = defaultdict(
        lambda: {
            "task_count": 0,
            "completed": 0,
            "cached": 0,
            "failed": 0,
            "other": 0,
            "peak_rss_max_bytes": None,
            "realtime_sum_seconds": 0.0,
            "duration_sum_seconds": 0.0,
            "cpu_percent_sum": 0.0,
            "cpu_observations": 0,
        }
    )
    missing_metric_counts = {"peak_rss": 0, "realtime": 0, "%cpu": 0, "duration": 0}
    with resolved.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        headers = set(reader.fieldnames or [])
        missing = sorted(required - headers)
        if missing:
            raise ResourceSummaryError(f"trace is missing required columns: {missing}")
        for index, row in enumerate(reader, start=1):
            process = _process_name(row)
            group = groups[process]
            group["task_count"] += 1
            state = row["status"].strip().upper()
            if state == "COMPLETED":
                group["completed"] += 1
            elif state == "CACHED":
                group["cached"] += 1
            elif state in {"FAILED", "ABORTED"}:
                group["failed"] += 1
            else:
                group["other"] += 1

            peak = parse_size(row["peak_rss"])
            if peak is None:
                missing_metric_counts["peak_rss"] += 1
            elif group["peak_rss_max_bytes"] is None or peak > group["peak_rss_max_bytes"]:
                group["peak_rss_max_bytes"] = peak
            realtime = parse_duration(row["realtime"])
            if realtime is None:
                missing_metric_counts["realtime"] += 1
            else:
                group["realtime_sum_seconds"] += realtime
            duration = parse_duration(row.get("duration", ""))
            if duration is None:
                missing_metric_counts["duration"] += 1
            else:
                group["duration_sum_seconds"] += duration
            cpu = parse_cpu(row["%cpu"])
            if cpu is None:
                missing_metric_counts["%cpu"] += 1
            else:
                group["cpu_percent_sum"] += cpu
                group["cpu_observations"] += 1
            if progress is not None:
                progress(index, total_rows)

    process_rows: list[dict[str, object]] = []
    for process in sorted(groups):
        group = groups[process]
        observations = int(group.pop("cpu_observations"))
        cpu_sum = float(group.pop("cpu_percent_sum"))
        process_rows.append(
            {
                "process": process,
                **group,
                "cpu_percent_mean": cpu_sum / observations if observations else None,
            }
        )
    task_count = sum(int(row["task_count"]) for row in process_rows)
    failed = sum(int(row["failed"]) for row in process_rows)
    return {
        "schema_version": 2,
        "trace": str(resolved),
        "task_count": task_count,
        "process_count": len(process_rows),
        "failed_task_count": failed,
        "all_tasks_successful": failed == 0 and task_count > 0,
        "missing_metric_counts": missing_metric_counts,
        "processes": process_rows,
        "metrics_are_estimates": True,
    }


def measure_storage(
    paths: Mapping[str, Path],
    *,
    progress: StorageProgressCallback | None = None,
) -> list[dict[str, object]]:
    """Measure logical and allocated bytes without following symlinks."""

    observations: list[dict[str, object]] = []
    for label in sorted(paths):
        root = paths[label].resolve(strict=True)
        if not root.is_dir():
            raise ResourceSummaryError(f"storage path is not a directory: {root}")
        file_count = 0
        directory_count = 0
        symlink_count = 0
        logical_bytes = 0
        allocated_bytes = 0
        if progress is not None:
            progress(label, 0)
        try:
            for current, directories, files in os.walk(root, followlinks=False):
                directory_count += 1
                current_path = Path(current)
                for name in directories:
                    if (current_path / name).is_symlink():
                        symlink_count += 1
                for name in files:
                    entry = current_path / name
                    metadata = entry.lstat()
                    if entry.is_symlink():
                        symlink_count += 1
                        continue
                    if not entry.is_file():
                        continue
                    file_count += 1
                    logical_bytes += metadata.st_size
                    allocated_bytes += metadata.st_blocks * 512
                    if progress is not None:
                        progress(label, file_count)
        except OSError as exc:
            raise ResourceSummaryError(f"could not measure storage path {root}: {exc}") from exc
        if progress is not None:
            progress(label, file_count)
        observations.append(
            {
                "label": label,
                "path": str(root),
                "file_count": file_count,
                "directory_count": directory_count,
                "symlink_count_ignored": symlink_count,
                "logical_bytes": logical_bytes,
                "allocated_bytes": allocated_bytes,
            }
        )
    return observations


def _exclusive_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    except FileExistsError as exc:
        raise ResourceSummaryError(f"refusing to overwrite existing output: {path}") from exc
    with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
        handle.write(content)


def write_resource_summary(
    *,
    trace: Path,
    json_output: Path,
    tsv_output: Path,
    progress: ProgressCallback | None = None,
    storage_paths: Mapping[str, Path] | None = None,
    storage_progress: StorageProgressCallback | None = None,
) -> dict[str, object]:
    """Create immutable JSON and per-process TSV resource summaries."""

    if json_output.exists() or tsv_output.exists():
        existing = json_output if json_output.exists() else tsv_output
        raise ResourceSummaryError(f"refusing to overwrite existing output: {existing}")
    summary = summarize_trace(trace, progress=progress)
    summary["storage"] = measure_storage(
        storage_paths or {},
        progress=storage_progress,
    )
    _exclusive_text(json_output, json.dumps(summary, indent=2, sort_keys=True) + "\n")
    headers = [
        "process", "task_count", "completed", "cached", "failed", "other",
        "peak_rss_max_bytes", "realtime_sum_seconds", "duration_sum_seconds",
        "cpu_percent_mean",
    ]
    lines = ["\t".join(headers)]
    for row in summary["processes"]:
        lines.append("\t".join("" if row[key] is None else str(row[key]) for key in headers))
    try:
        _exclusive_text(tsv_output, "\n".join(lines) + "\n")
    except Exception:
        json_output.unlink(missing_ok=True)
        raise
    return summary
