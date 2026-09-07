"""Validate and seal a completed nf-core/rnaseq primary-processing run."""

from __future__ import annotations

import csv
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

from .plan import sha256_file
from .resources import ResourceSummaryError, measure_storage, summarize_trace


HashProgress = Callable[[str, int, int], None]
LogProgress = Callable[[str, int, bool], None]


class PrimaryDeliveryError(ValueError):
    """Raised when primary-processing completion cannot be trusted."""


def _read_verified_plan(path: Path) -> tuple[Path, dict[str, object]]:
    try:
        resolved = path.resolve(strict=True)
        payload = json.loads(resolved.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PrimaryDeliveryError(f"could not read run plan: {exc}") from exc
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise PrimaryDeliveryError("unsupported run plan schema")
    if payload.get("stage") != "planned_not_executed":
        raise PrimaryDeliveryError("run plan has an unsupported stage")
    workflow = payload.get("workflow")
    if not isinstance(workflow, dict) or workflow.get("name") != "nf-core/rnaseq":
        raise PrimaryDeliveryError("run plan is not for nf-core/rnaseq")
    reference = payload.get("reference")
    if (
        not isinstance(reference, dict)
        or reference.get("mode") != "igenomes"
        or not isinstance(reference.get("genome"), str)
    ):
        raise PrimaryDeliveryError("run plan has no supported reviewed reference genome")
    execution = payload.get("execution")
    if not isinstance(execution, dict) or not isinstance(execution.get("outdir"), str):
        raise PrimaryDeliveryError("run plan execution record is invalid")
    controls = payload.get("control_files")
    required_controls = {
        "samplesheet", "design", "contrasts", "workflow_lock", "input_manifest"
    }
    if not isinstance(controls, dict) or not required_controls.issubset(controls):
        raise PrimaryDeliveryError("run plan control-file records are incomplete")
    for label, record in controls.items():
        if not isinstance(record, dict):
            raise PrimaryDeliveryError(f"invalid control file record: {label}")
        try:
            control = Path(str(record.get("path", ""))).resolve(strict=True)
        except OSError as exc:
            raise PrimaryDeliveryError(f"control file is unavailable: {label}") from exc
        if not control.is_file() or control.stat().st_size != record.get("size_bytes"):
            raise PrimaryDeliveryError(f"control file size changed after planning: {label}")
        if sha256_file(control) != record.get("sha256"):
            raise PrimaryDeliveryError(f"control file checksum changed after planning: {label}")
    argv = payload.get("command_argv")
    if not isinstance(argv, list) or not all(isinstance(item, str) for item in argv):
        raise PrimaryDeliveryError("run plan command_argv is invalid")
    expected_output = str(Path(execution["outdir"]).resolve() / "rnaseq")
    try:
        output_index = argv.index("--outdir")
    except ValueError as exc:
        raise PrimaryDeliveryError("run plan has no nf-core output argument") from exc
    if output_index + 1 >= len(argv) or str(Path(argv[output_index + 1]).resolve()) != expected_output:
        raise PrimaryDeliveryError("run plan output argument conflicts with its execution record")
    return resolved, payload


def _required_file(path: Path, role: str) -> Path:
    if not path.is_file():
        raise PrimaryDeliveryError(f"missing required {role}: {path}")
    if path.stat().st_size == 0:
        raise PrimaryDeliveryError(f"empty required {role}: {path}")
    return path.resolve()


def _one_match(root: Path, pattern: str, role: str) -> Path:
    matches = sorted(path for path in root.glob(pattern) if path.is_file())
    if len(matches) != 1:
        raise PrimaryDeliveryError(
            f"expected exactly one {role} matching {root / pattern}; found {len(matches)}"
        )
    return _required_file(matches[0], role)


def _samples(samplesheet: Path) -> set[str]:
    try:
        with samplesheet.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None or "sample" not in reader.fieldnames:
                raise PrimaryDeliveryError("planned samplesheet has no sample column")
            return {row["sample"].strip() for row in reader if row.get("sample", "").strip()}
    except (OSError, csv.Error) as exc:
        raise PrimaryDeliveryError(f"could not read planned samplesheet: {exc}") from exc


def _check_matrix(path: Path, expected_samples: set[str], role: str) -> dict[str, int]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            header = handle.readline().rstrip("\r\n").split("\t")
            first_data_row = handle.readline()
    except (OSError, UnicodeError) as exc:
        raise PrimaryDeliveryError(f"could not read {role}: {exc}") from exc
    missing = sorted(expected_samples - set(header))
    if missing:
        raise PrimaryDeliveryError(
            f"{role} is missing {len(missing)} planned sample column(s): {', '.join(missing[:5])}"
        )
    if not first_data_row:
        raise PrimaryDeliveryError(f"{role} contains no feature rows")
    return {"sample_columns": len(expected_samples), "header_columns": len(header)}


def _check_params(
    path: Path,
    samplesheet: Path,
    rnaseq_root: Path,
    genome: str,
) -> None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PrimaryDeliveryError(f"could not read pipeline params.json: {exc}") from exc
    if not isinstance(payload, dict):
        raise PrimaryDeliveryError("pipeline params.json must contain a JSON object")
    for key, expected in (("input", samplesheet), ("outdir", rnaseq_root)):
        value = payload.get(key)
        if not isinstance(value, str) or Path(value).resolve() != expected.resolve():
            raise PrimaryDeliveryError(f"pipeline params.json {key!r} does not match the run plan")
    if payload.get("genome") != genome:
        raise PrimaryDeliveryError("pipeline params.json 'genome' does not match the run plan")


def _hash_file(path: Path, role: str, progress: HashProgress | None) -> str:
    import hashlib

    total = path.stat().st_size
    current = 0
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
            current += len(chunk)
            if progress is not None:
                progress(role, current, total)
    if progress is not None and total == 0:
        progress(role, 0, 0)
    return digest.hexdigest()


def inspect_primary_completion(
    run_plan: Path,
    *,
    progress: HashProgress | None = None,
) -> dict[str, object]:
    """Verify core outputs and return the receipt payload without writing it."""

    plan_path, plan = _read_verified_plan(run_plan)
    execution = plan["execution"]
    run_root = Path(str(execution["outdir"])).resolve()
    rnaseq_root = run_root / "rnaseq"
    samplesheet = Path(str(plan["control_files"]["samplesheet"]["path"])).resolve()
    trace = _required_file(run_root / "execution" / "trace.tsv", "Nextflow trace")
    try:
        trace_summary = summarize_trace(trace)
    except ResourceSummaryError as exc:
        raise PrimaryDeliveryError(str(exc)) from exc
    if not trace_summary["all_tasks_successful"]:
        raise PrimaryDeliveryError(
            "Nextflow trace is not complete and successful: "
            f"failed={trace_summary['failed_task_count']}, "
            f"nonterminal_or_unknown={trace_summary['other_task_count']}"
        )

    artifacts: list[tuple[str, Path]] = [
        ("nextflow_report", _required_file(run_root / "execution" / "report.html", "Nextflow report")),
        ("nextflow_trace", trace),
        ("nextflow_timeline", _required_file(run_root / "execution" / "timeline.html", "Nextflow timeline")),
        ("nextflow_dag", _required_file(run_root / "execution" / "dag.html", "Nextflow DAG")),
        ("software_versions", _required_file(rnaseq_root / "pipeline_info" / "software_versions.yml", "software versions")),
        ("pipeline_params", _required_file(rnaseq_root / "pipeline_info" / "params.json", "pipeline parameters")),
        ("validated_samplesheet", _required_file(rnaseq_root / "pipeline_info" / "samplesheet.valid.csv", "validated samplesheet")),
        ("gene_counts", _required_file(rnaseq_root / "star_salmon" / "salmon.merged.gene_counts.tsv", "gene-count matrix")),
        ("gene_tpm", _required_file(rnaseq_root / "star_salmon" / "salmon.merged.gene_tpm.tsv", "gene-TPM matrix")),
        ("multiqc_report", _one_match(rnaseq_root, "multiqc/*/multiqc_report.html", "MultiQC report")),
    ]
    by_role = dict(artifacts)
    multiqc_data = by_role["multiqc_report"].parent / "multiqc_data"
    if not multiqc_data.is_dir():
        raise PrimaryDeliveryError(f"missing required MultiQC data directory: {multiqc_data}")
    storage = measure_storage({"multiqc_data": multiqc_data})
    if storage[0]["file_count"] == 0:
        raise PrimaryDeliveryError("MultiQC data directory contains no files")

    expected_samples = _samples(samplesheet)
    if not expected_samples:
        raise PrimaryDeliveryError("planned samplesheet contains no samples")
    matrix_checks = {
        "gene_counts": _check_matrix(by_role["gene_counts"], expected_samples, "gene-count matrix"),
        "gene_tpm": _check_matrix(by_role["gene_tpm"], expected_samples, "gene-TPM matrix"),
    }
    _check_params(
        by_role["pipeline_params"],
        samplesheet,
        rnaseq_root,
        str(plan["reference"]["genome"]),
    )

    inventory = []
    for role, path in artifacts:
        inventory.append(
            {
                "role": role,
                "path": str(path),
                "relative_path": str(path.relative_to(run_root)),
                "size_bytes": path.stat().st_size,
                "sha256": _hash_file(path, role, progress),
            }
        )
    return {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "stage": "pipeline_complete_qc_pending",
        "run_plan": {
            "path": str(plan_path),
            "size_bytes": plan_path.stat().st_size,
            "sha256": sha256_file(plan_path),
        },
        "workflow": plan["workflow"],
        "runtime": plan["runtime"],
        "reference": plan["reference"],
        "controls_verified": True,
        "task_summary": {
            "task_count": trace_summary["task_count"],
            "process_count": trace_summary["process_count"],
            "failed_task_count": trace_summary["failed_task_count"],
            "other_task_count": trace_summary["other_task_count"],
        },
        "sample_count": len(expected_samples),
        "matrix_checks": matrix_checks,
        "artifacts": inventory,
        "multiqc": {
            "report": str(by_role["multiqc_report"]),
            "data_directory": str(multiqc_data.resolve()),
            "data_storage": storage[0],
        },
        "qc_status": "pending_review",
        "contains_client_results": True,
        "contains_secrets": False,
    }


def create_primary_receipt(
    *,
    run_plan: Path,
    output: Path,
    progress: HashProgress | None = None,
) -> dict[str, object]:
    """Exclusively write a completion receipt after all validation succeeds."""

    if output.exists():
        raise PrimaryDeliveryError(f"refusing to overwrite existing receipt: {output}")
    payload = inspect_primary_completion(run_plan, progress=progress)
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        descriptor = os.open(output, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    except FileExistsError as exc:
        raise PrimaryDeliveryError(f"refusing to overwrite existing receipt: {output}") from exc
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
    except Exception:
        output.unlink(missing_ok=True)
        raise
    return payload


DIAGNOSTIC_PATTERNS = {
    "out_of_memory": re.compile(r"out.of.memory|oom.kill|killed process|exceeded.*memory", re.I),
    "time_limit": re.compile(r"time.?limit|\btimeout\b", re.I),
    "storage": re.compile(r"no space left|disk quota|quota exceeded", re.I),
    "network_or_registry": re.compile(
        r"unknownhost|connection (?:timed out|refused)|ssl(?:error|exception)|registry.*error",
        re.I,
    ),
    "missing_or_denied_input": re.compile(
        r"no such file|file not found|permission denied", re.I
    ),
    "process_failure": re.compile(r"error\s*~|process .+ terminated|exit status", re.I),
}

DIAGNOSTIC_RECOMMENDATIONS = {
    "out_of_memory": "Review the failed process trace and approved task memory ceiling before retrying.",
    "time_limit": "Review the failed process duration and approved task time ceiling before retrying.",
    "storage": "Confirm quota and free space for the work and results roots before retrying.",
    "network_or_registry": "Verify the selected direct, proxy, or offline staging mode before retrying.",
    "missing_or_denied_input": "Restore the planned file/path permissions; do not replace a hashed control file.",
    "process_failure": "Inspect the named process work directory and tool log before retrying.",
}


def _read_submission_receipt(path: Path) -> tuple[Path, dict[str, object]]:
    try:
        resolved = path.resolve(strict=True)
        payload = json.loads(resolved.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PrimaryDeliveryError(f"could not read submission receipt: {exc}") from exc
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise PrimaryDeliveryError("unsupported submission receipt schema")
    return resolved, payload


def _diagnostic_sources(
    *,
    submission_receipt: Path | None,
    nextflow_log: Path | None,
) -> list[Path]:
    candidates: list[Path] = []
    if nextflow_log is not None:
        candidates.append(nextflow_log.resolve(strict=True))
    if submission_receipt is not None:
        _, receipt = _read_submission_receipt(submission_receipt)
        launcher_value = receipt.get("launcher")
        if not isinstance(launcher_value, str):
            raise PrimaryDeliveryError("submission receipt has no launcher path")
        launch_dir = Path(launcher_value).resolve().parent
        candidates.append(launch_dir / ".nextflow.log")
        candidates.extend(sorted((launch_dir / "logs").glob("controller-*.err")))
        candidates.extend(sorted((launch_dir / "logs").glob("controller-*.out")))
    unique: dict[str, Path] = {}
    for candidate in candidates:
        if candidate.is_file():
            unique[str(candidate.resolve())] = candidate.resolve()
    return list(unique.values())


def diagnose_primary_run(
    *,
    run_plan: Path,
    submission_receipt: Path | None = None,
    nextflow_log: Path | None = None,
    progress: LogProgress | None = None,
) -> dict[str, object]:
    """Classify an incomplete run without serializing raw log content."""

    plan_path, plan = _read_verified_plan(run_plan)
    execution = plan["execution"]
    run_root = Path(str(execution["outdir"])).resolve()
    trace = run_root / "execution" / "trace.tsv"
    trace_summary: dict[str, object] | None = None
    trace_error: str | None = None
    if trace.is_file():
        try:
            trace_summary = summarize_trace(trace)
        except ResourceSummaryError as exc:
            trace_error = str(exc)

    findings: dict[str, dict[str, object]] = {}
    sources = _diagnostic_sources(
        submission_receipt=submission_receipt,
        nextflow_log=nextflow_log,
    )
    for source in sources:
        counts = {category: 0 for category in DIAGNOSTIC_PATTERNS}
        line_numbers: dict[str, list[int]] = {category: [] for category in DIAGNOSTIC_PATTERNS}
        try:
            with source.open("r", encoding="utf-8", errors="replace") as handle:
                final_line = 0
                for final_line, line in enumerate(handle, start=1):
                    for category, pattern in DIAGNOSTIC_PATTERNS.items():
                        if pattern.search(line):
                            counts[category] += 1
                            if len(line_numbers[category]) < 20:
                                line_numbers[category].append(final_line)
                    if progress is not None:
                        progress(source.name, final_line, False)
                if progress is not None:
                    progress(source.name, final_line, True)
        except OSError as exc:
            raise PrimaryDeliveryError(f"could not scan diagnostic log {source}: {exc}") from exc
        for category, count in counts.items():
            if count:
                record = findings.setdefault(
                    category,
                    {
                        "match_count": 0,
                        "locations": [],
                        "recommendation": DIAGNOSTIC_RECOMMENDATIONS[category],
                    },
                )
                record["match_count"] += count
                record["locations"].append(
                    {"path": str(source), "line_numbers": line_numbers[category]}
                )

    if trace_summary is None:
        status = "incomplete_no_valid_trace"
    elif trace_summary["all_tasks_successful"]:
        status = "pipeline_trace_complete_run_completion_not_yet_sealed"
    else:
        status = "failed_or_incomplete"
    return {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "run_plan": {"path": str(plan_path), "sha256": sha256_file(plan_path)},
        "trace": str(trace) if trace.exists() else None,
        "trace_summary": trace_summary,
        "trace_error": trace_error,
        "log_sources": [str(source) for source in sources],
        "findings": findings,
        "raw_log_content_recorded": False,
        "contains_client_results": True,
        "contains_secrets": False,
    }


def prepare_safe_restart(
    *,
    run_plan: Path,
    settings_path: Path,
    previous_receipt: Path,
    output: Path,
) -> dict[str, object]:
    """Prepare, but never submit, a new launcher for a safely resumable failed run."""

    from .hpc_run import HPCRunError, prepare_launcher, query_job

    _, plan = _read_verified_plan(run_plan)
    execution = plan["execution"]
    if execution.get("executor") != "slurm":
        raise PrimaryDeliveryError("safe restart currently requires a SLURM run plan")
    receipt_path, receipt = _read_submission_receipt(previous_receipt)
    launcher_value = receipt.get("launcher")
    launcher_hash = receipt.get("launcher_sha256")
    if not isinstance(launcher_value, str) or not isinstance(launcher_hash, str):
        raise PrimaryDeliveryError("previous receipt has no launcher identity")
    try:
        launcher = Path(launcher_value).resolve(strict=True)
    except OSError as exc:
        raise PrimaryDeliveryError("previous launcher is unavailable") from exc
    if sha256_file(launcher) != launcher_hash:
        raise PrimaryDeliveryError("previous launcher checksum has changed")

    receipt_status = receipt.get("status")
    scheduler: dict[str, object] | None
    if receipt_status == "submitted":
        try:
            scheduler = query_job(receipt_path)
        except HPCRunError as exc:
            raise PrimaryDeliveryError(str(exc)) from exc
        if not scheduler.get("terminal"):
            raise PrimaryDeliveryError("previous scheduler job is not terminal")
        if scheduler.get("successful"):
            raise PrimaryDeliveryError("previous scheduler job completed successfully; do not restart it")
    elif receipt_status == "submission_failed":
        scheduler = None
    else:
        raise PrimaryDeliveryError(f"previous receipt status is not restartable: {receipt_status!r}")

    try:
        workdir = Path(str(execution["workdir"])).resolve(strict=True)
    except OSError as exc:
        raise PrimaryDeliveryError("planned work directory is unavailable for -resume") from exc
    if not workdir.is_dir():
        raise PrimaryDeliveryError("planned work path is not a directory")
    trace = Path(str(execution["outdir"])).resolve() / "execution" / "trace.tsv"
    if trace.is_file():
        try:
            summary = summarize_trace(trace)
        except ResourceSummaryError as exc:
            raise PrimaryDeliveryError(str(exc)) from exc
        if summary["all_tasks_successful"]:
            raise PrimaryDeliveryError("Nextflow trace is already complete; do not restart it")
    try:
        result = prepare_launcher(
            run_plan=run_plan,
            settings_path=settings_path,
            output=output,
        )
    except HPCRunError as exc:
        raise PrimaryDeliveryError(str(exc)) from exc
    return {
        **result,
        "restart_prepared": True,
        "restart_of_receipt": str(receipt_path),
        "restart_of_job_id": receipt.get("job_id"),
        "previous_scheduler_state": scheduler.get("state") if scheduler else "SUBMISSION_FAILED",
        "workdir_reused": str(workdir),
        "submitted": False,
    }
