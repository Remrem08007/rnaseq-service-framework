"""Validate and seal a completed nf-core/rnaseq primary-processing run."""

from __future__ import annotations

import csv
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

from .plan import sha256_file
from .resources import ResourceSummaryError, measure_storage, summarize_trace


HashProgress = Callable[[str, int, int], None]


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
    execution = payload.get("execution")
    if not isinstance(execution, dict) or not isinstance(execution.get("outdir"), str):
        raise PrimaryDeliveryError("run plan execution record is invalid")
    controls = payload.get("control_files")
    required_controls = {"samplesheet", "design", "contrasts", "workflow_lock"}
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


def _check_params(path: Path, samplesheet: Path, rnaseq_root: Path) -> None:
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
    _check_params(by_role["pipeline_params"], samplesheet, rnaseq_root)

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
