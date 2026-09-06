from __future__ import annotations

import json
import stat
from pathlib import Path

import pytest

from rnaseq_service.resources import (
    ResourceSummaryError,
    parse_duration,
    parse_size,
    measure_storage,
    summarize_trace,
    write_resource_summary,
)


def trace_file(tmp_path: Path) -> Path:
    path = tmp_path / "trace.tsv"
    path.write_text(
        "name\tstatus\tduration\trealtime\t%cpu\tpeak_rss\n"
        "FASTQC (sample_1)\tCOMPLETED\t1m 30s\t1m 20s\t95.0%\t500 MB\n"
        "FASTQC (sample_2)\tCACHED\t2s\t1s\t20%\t100 MB\n"
        "STAR_ALIGN (sample_1)\tFAILED\t2h\t1h 30m\t350%\t8 GB\n",
        encoding="utf-8",
    )
    return path


def test_unit_parsers() -> None:
    assert parse_size("1.5 GB") == round(1.5 * 1024**3)
    assert parse_duration("1h 2m 3.5s") == 3723.5
    assert parse_duration("352ms") == pytest.approx(0.352)


def test_trace_summary_groups_processes_and_failures(tmp_path: Path) -> None:
    progress = []

    summary = summarize_trace(trace_file(tmp_path), progress=lambda row, total: progress.append((row, total)))

    assert summary["task_count"] == 3
    assert summary["process_count"] == 2
    assert summary["failed_task_count"] == 1
    assert summary["all_tasks_successful"] is False
    fastqc = next(row for row in summary["processes"] if row["process"] == "FASTQC")
    assert fastqc["task_count"] == 2
    assert fastqc["cached"] == 1
    assert fastqc["peak_rss_max_bytes"] == 500 * 1024**2
    assert progress[-1] == (3, 3)


def test_summary_outputs_are_private_and_immutable(tmp_path: Path) -> None:
    json_out = tmp_path / "summary" / "resources.json"
    tsv_out = tmp_path / "summary" / "resources.tsv"

    summary = write_resource_summary(
        trace=trace_file(tmp_path),
        json_output=json_out,
        tsv_output=tsv_out,
    )

    assert json.loads(json_out.read_text(encoding="utf-8"))["task_count"] == 3
    assert "STAR_ALIGN" in tsv_out.read_text(encoding="utf-8")
    assert summary["metrics_are_estimates"] is True
    assert summary["storage"] == []
    assert stat.S_IMODE(json_out.stat().st_mode) == 0o600
    assert stat.S_IMODE(tsv_out.stat().st_mode) == 0o600
    with pytest.raises(ResourceSummaryError, match="refusing to overwrite"):
        write_resource_summary(
            trace=trace_file(tmp_path),
            json_output=json_out,
            tsv_output=tsv_out,
        )


def test_missing_required_trace_column_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "bad.tsv"
    path.write_text("name\tstatus\nFASTQC (x)\tCOMPLETED\n", encoding="utf-8")

    with pytest.raises(ResourceSummaryError, match="missing required columns"):
        summarize_trace(path)


def test_storage_measurement_counts_bytes_and_ignores_symlinks(tmp_path: Path) -> None:
    root = tmp_path / "work"
    nested = root / "nested"
    nested.mkdir(parents=True)
    (root / "one.bin").write_bytes(b"1234")
    (nested / "two.bin").write_bytes(b"123456")
    (root / "link.bin").symlink_to(nested / "two.bin")
    progress: list[tuple[str, int]] = []

    observations = measure_storage(
        {"work": root},
        progress=lambda label, count: progress.append((label, count)),
    )

    assert observations[0]["file_count"] == 2
    assert observations[0]["directory_count"] == 2
    assert observations[0]["symlink_count_ignored"] == 1
    assert observations[0]["logical_bytes"] == 10
    assert observations[0]["allocated_bytes"] >= 0
    assert progress[-1] == ("work", 2)


def test_summary_includes_requested_storage_observations(tmp_path: Path) -> None:
    results = tmp_path / "results"
    results.mkdir()
    (results / "output.txt").write_text("result", encoding="utf-8")
    json_out = tmp_path / "summary" / "resources.json"
    tsv_out = tmp_path / "summary" / "resources.tsv"

    summary = write_resource_summary(
        trace=trace_file(tmp_path),
        json_output=json_out,
        tsv_output=tsv_out,
        storage_paths={"results": results},
    )

    assert summary["storage"][0]["label"] == "results"
    assert summary["storage"][0]["logical_bytes"] == 6
