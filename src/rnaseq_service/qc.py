"""Evaluate configurable sample-QC rules against checksum-bound MultiQC evidence."""

from __future__ import annotations

import csv
import json
import math
import os
import re
import shutil
import statistics
import tomllib
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

from .plan import sha256_file


QCProgress = Callable[[str, int, int], None]
METRIC_ID = re.compile(r"^[a-z][a-z0-9_]*$")


class QCError(ValueError):
    """Raised when QC evidence or policy cannot be trusted."""


def _read_json(path: Path, label: str) -> tuple[Path, dict[str, object]]:
    try:
        resolved = path.resolve(strict=True)
        payload = json.loads(resolved.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise QCError(f"could not read {label}: {exc}") from exc
    if not isinstance(payload, dict):
        raise QCError(f"{label} must contain a JSON object")
    return resolved, payload


def _read_completion(path: Path) -> tuple[Path, dict[str, object], dict[str, dict[str, object]]]:
    resolved, payload = _read_json(path, "primary completion receipt")
    if payload.get("schema_version") != 2:
        raise QCError("primary completion receipt must use schema version 2")
    if payload.get("stage") != "pipeline_complete_qc_pending":
        raise QCError("primary completion receipt is not awaiting QC")
    if payload.get("qc_status") != "pending_review":
        raise QCError("primary completion receipt has an unsupported QC status")
    raw_artifacts = payload.get("artifacts")
    if not isinstance(raw_artifacts, list):
        raise QCError("primary completion receipt has no artifact inventory")
    artifacts: dict[str, dict[str, object]] = {}
    for record in raw_artifacts:
        if not isinstance(record, dict) or not isinstance(record.get("role"), str):
            raise QCError("primary completion receipt has an invalid artifact record")
        role = record["role"]
        if role in artifacts:
            raise QCError(f"duplicate artifact role in completion receipt: {role}")
        artifacts[role] = record
    for role in ("multiqc_data_json", "multiqc_general_stats", "validated_samplesheet"):
        if role not in artifacts:
            raise QCError(f"primary completion receipt is missing artifact role: {role}")
    return resolved, payload, artifacts


def _verified_artifact(record: dict[str, object], role: str) -> Path:
    try:
        path = Path(str(record["path"])).resolve(strict=True)
    except (KeyError, OSError) as exc:
        raise QCError(f"completion artifact is unavailable: {role}") from exc
    if not path.is_file() or path.stat().st_size != record.get("size_bytes"):
        raise QCError(f"completion artifact size changed: {role}")
    if sha256_file(path) != record.get("sha256"):
        raise QCError(f"completion artifact checksum changed: {role}")
    return path


def _load_policy(path: Path) -> tuple[Path, dict[str, object], list[dict[str, object]]]:
    try:
        resolved = path.resolve(strict=True)
        payload = tomllib.loads(resolved.read_text(encoding="utf-8"))
    except (OSError, tomllib.TOMLDecodeError) as exc:
        raise QCError(f"could not read QC policy: {exc}") from exc
    if payload.get("schema_version") != 1:
        raise QCError("QC policy must use schema_version = 1")
    study = payload.get("study", {})
    if not isinstance(study, dict):
        raise QCError("QC policy [study] must be a table")
    robust_z = study.get("robust_z_threshold", 4.5)
    minimum = study.get("min_samples_for_outliers", 4)
    if isinstance(robust_z, bool) or not isinstance(robust_z, (int, float)) or robust_z <= 0:
        raise QCError("study.robust_z_threshold must be positive")
    if isinstance(minimum, bool) or not isinstance(minimum, int) or minimum < 3:
        raise QCError("study.min_samples_for_outliers must be an integer >= 3")
    metrics = payload.get("metrics")
    if not isinstance(metrics, list) or not metrics:
        raise QCError("QC policy must define at least one [[metrics]] table")
    seen: set[str] = set()
    normalized: list[dict[str, object]] = []
    for index, metric in enumerate(metrics, start=1):
        if not isinstance(metric, dict):
            raise QCError(f"metrics entry {index} must be a table")
        metric_id = metric.get("id")
        if not isinstance(metric_id, str) or METRIC_ID.fullmatch(metric_id) is None:
            raise QCError(f"metrics entry {index} has an invalid id")
        if metric_id in seen:
            raise QCError(f"duplicate metric id: {metric_id}")
        seen.add(metric_id)
        for key in ("source_file", "sample_column", "value_column"):
            if not isinstance(metric.get(key), str) or not str(metric[key]).strip():
                raise QCError(f"metric {metric_id!r} requires {key}")
        required = metric.get("required", True)
        outlier = metric.get("outlier", True)
        if not isinstance(required, bool) or not isinstance(outlier, bool):
            raise QCError(f"metric {metric_id!r} required/outlier must be booleans")
        kind = metric.get("kind", "numeric")
        if kind not in {"numeric", "categorical"}:
            raise QCError(f"metric {metric_id!r} kind must be 'numeric' or 'categorical'")
        lower = metric.get("lower")
        upper = metric.get("upper")
        for key, value in (("lower", lower), ("upper", upper)):
            if value is not None and (isinstance(value, bool) or not isinstance(value, (int, float))):
                raise QCError(f"metric {metric_id!r} {key} must be numeric")
        if lower is not None and upper is not None and lower > upper:
            raise QCError(f"metric {metric_id!r} lower exceeds upper")
        allowed_values = metric.get("allowed_values")
        if kind == "numeric" and allowed_values is not None:
            raise QCError(f"numeric metric {metric_id!r} cannot define allowed_values")
        if kind == "categorical":
            if lower is not None or upper is not None or outlier:
                raise QCError(
                    f"categorical metric {metric_id!r} cannot define numeric limits or outlier = true"
                )
            if (
                not isinstance(allowed_values, list)
                or not allowed_values
                or not all(isinstance(value, str) and value for value in allowed_values)
            ):
                raise QCError(f"categorical metric {metric_id!r} requires allowed_values")
        normalized.append(dict(metric))
    return resolved, {
        "robust_z_threshold": float(robust_z),
        "min_samples_for_outliers": minimum,
    }, normalized


def _planned_samples(samplesheet: Path) -> list[str]:
    try:
        with samplesheet.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None or "sample" not in reader.fieldnames:
                raise QCError("validated samplesheet has no sample column")
            samples = [row["sample"].strip() for row in reader if row.get("sample", "").strip()]
    except (OSError, csv.Error) as exc:
        raise QCError(f"could not read validated samplesheet: {exc}") from exc
    if not samples or len(samples) != len(set(samples)):
        raise QCError("validated samplesheet must contain unique, non-empty samples")
    return samples


def _evidence_path(root: Path, relative: str) -> Path:
    candidate = (root / relative).resolve()
    try:
        candidate.relative_to(root.resolve())
    except ValueError as exc:
        raise QCError(f"QC policy source_file escapes the MultiQC directory: {relative!r}") from exc
    return candidate


def _bound_multiqc_files(
    completion: dict[str, object], root: Path
) -> dict[str, dict[str, object]]:
    multiqc = completion.get("multiqc")
    if not isinstance(multiqc, dict) or not isinstance(multiqc.get("files"), list):
        raise QCError("primary completion receipt has no checksum-bound MultiQC inventory")
    bound: dict[str, dict[str, object]] = {}
    for record in multiqc["files"]:
        if not isinstance(record, dict) or not isinstance(record.get("relative_path"), str):
            raise QCError("primary completion receipt has an invalid MultiQC file record")
        relative = record["relative_path"]
        path = _evidence_path(root, relative)
        if str(path.relative_to(root)) != relative or relative in bound:
            raise QCError("primary completion receipt has an invalid MultiQC relative path")
        bound[relative] = record
    return bound


def _read_metric(
    root: Path,
    metric: dict[str, object],
    planned: set[str],
    bound_files: dict[str, dict[str, object]],
) -> tuple[dict[str, float | str], dict[str, object], str | None]:
    metric_id = str(metric["id"])
    path = _evidence_path(root, str(metric["source_file"]))
    evidence: dict[str, object] = {"path": str(path), "exists": path.is_file()}
    relative = str(path.relative_to(root))
    bound = bound_files.get(relative)
    if bound is None:
        evidence["bound_by_completion"] = False
        return {}, evidence, f"{metric_id}: source file was not bound by the completion receipt"
    evidence["bound_by_completion"] = True
    if not path.is_file():
        return {}, evidence, f"{metric_id}: missing source file {metric['source_file']}"
    if path.stat().st_size != bound.get("size_bytes"):
        raise QCError(f"checksum-bound MultiQC evidence size changed: {relative}")
    if sha256_file(path) != bound.get("sha256"):
        raise QCError(f"checksum-bound MultiQC evidence checksum changed: {relative}")
    evidence["size_bytes"] = path.stat().st_size
    evidence["sha256"] = sha256_file(path)
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            columns = reader.fieldnames or []
            evidence["columns"] = columns
            required_columns = {str(metric["sample_column"]), str(metric["value_column"])}
            missing_columns = sorted(required_columns - set(columns))
            if missing_columns:
                return {}, evidence, f"{metric_id}: missing column(s) {', '.join(missing_columns)}"
            values: dict[str, float | str] = {}
            invalid: list[str] = []
            for line_number, row in enumerate(reader, start=2):
                sample = (row.get(str(metric["sample_column"])) or "").strip()
                if sample not in planned:
                    continue
                raw = (row.get(str(metric["value_column"])) or "").strip()
                if not raw:
                    continue
                if sample in values:
                    raise QCError(f"metric {metric_id!r} contains duplicate sample {sample!r}")
                if metric.get("kind", "numeric") == "numeric":
                    try:
                        value: float | str = float(raw)
                    except ValueError:
                        invalid.append(f"{sample}@{line_number}")
                        continue
                    if not math.isfinite(value):
                        invalid.append(f"{sample}@{line_number}")
                        continue
                else:
                    value = raw
                values[sample] = value
    except (OSError, csv.Error, UnicodeError) as exc:
        raise QCError(f"could not read metric {metric_id!r}: {exc}") from exc
    if invalid:
        return values, evidence, f"{metric_id}: non-numeric value(s) for {', '.join(invalid[:5])}"
    return values, evidence, None


def _robust_scores(values: dict[str, float], minimum: int) -> dict[str, float]:
    if len(values) < minimum:
        return {}
    center = statistics.median(values.values())
    deviations = [abs(value - center) for value in values.values()]
    mad = statistics.median(deviations)
    if mad == 0:
        return {}
    return {sample: 0.67448975 * (value - center) / mad for sample, value in values.items()}


def evaluate_qc(
    *,
    completion_receipt: Path,
    policy_path: Path,
    outdir: Path,
    progress: QCProgress | None = None,
) -> dict[str, object]:
    """Create a non-overwriting QC assessment and human decision template."""

    if outdir.exists():
        raise QCError(f"refusing to overwrite existing QC directory: {outdir}")
    completion_path, completion, artifacts = _read_completion(completion_receipt)
    policy_resolved, study, metrics = _load_policy(policy_path)
    general_stats = _verified_artifact(artifacts["multiqc_general_stats"], "multiqc_general_stats")
    _verified_artifact(artifacts["multiqc_data_json"], "multiqc_data_json")
    samplesheet = _verified_artifact(artifacts["validated_samplesheet"], "validated_samplesheet")
    multiqc_root = Path(str(completion["multiqc"]["data_directory"])).resolve(strict=True)
    if general_stats.parent != multiqc_root:
        raise QCError("completion receipt MultiQC directory conflicts with its artifact inventory")
    bound_files = _bound_multiqc_files(completion, multiqc_root)
    samples = _planned_samples(samplesheet)
    sample_set = set(samples)
    rows = {
        sample: {"sample": sample, "metrics": {}, "missing_metrics": [], "flags": []}
        for sample in samples
    }
    evidence: dict[str, dict[str, object]] = {}
    evidence_issues: list[dict[str, object]] = []
    for index, metric in enumerate(metrics, start=1):
        metric_id = str(metric["id"])
        if progress is not None:
            progress(metric_id, index - 1, len(metrics))
        values, record, issue = _read_metric(
            multiqc_root, metric, sample_set, bound_files
        )
        evidence[metric_id] = record
        if issue is not None:
            evidence_issues.append(
                {"metric": metric_id, "required": bool(metric.get("required", True)), "issue": issue}
            )
        for sample in samples:
            value = values.get(sample)
            rows[sample]["metrics"][metric_id] = value
            if value is None:
                rows[sample]["missing_metrics"].append(metric_id)
                continue
            lower = metric.get("lower")
            upper = metric.get("upper")
            if metric.get("kind", "numeric") == "categorical":
                allowed = set(metric["allowed_values"])
                if value not in allowed:
                    rows[sample]["flags"].append(
                        {
                            "metric": metric_id,
                            "kind": "unexpected_category",
                            "value": value,
                            "allowed_values": sorted(allowed),
                        }
                    )
                continue
            numeric_value = float(value)
            if lower is not None and numeric_value < float(lower):
                rows[sample]["flags"].append(
                    {"metric": metric_id, "kind": "below_lower", "value": numeric_value, "threshold": float(lower)}
                )
            if upper is not None and numeric_value > float(upper):
                rows[sample]["flags"].append(
                    {"metric": metric_id, "kind": "above_upper", "value": numeric_value, "threshold": float(upper)}
                )
        if bool(metric.get("outlier", True)):
            numeric_values = {sample: float(value) for sample, value in values.items()}
            scores = _robust_scores(numeric_values, int(study["min_samples_for_outliers"]))
            for sample, score in scores.items():
                if abs(score) > float(study["robust_z_threshold"]):
                    rows[sample]["flags"].append(
                        {
                            "metric": metric_id,
                            "kind": "robust_outlier",
                            "value": numeric_values[sample],
                            "robust_z": score,
                            "threshold": float(study["robust_z_threshold"]),
                        }
                    )
        if progress is not None:
            progress(metric_id, index, len(metrics))
    required_ids = {str(metric["id"]) for metric in metrics if bool(metric.get("required", True))}
    missing_required_samples = [
        sample for sample in samples if required_ids.intersection(rows[sample]["missing_metrics"])
    ]
    required_evidence_issues = [issue for issue in evidence_issues if issue["required"]]
    status = "stop_missing_evidence" if missing_required_samples or required_evidence_issues else "review_required"
    flagged = [sample for sample in samples if rows[sample]["flags"]]
    payload: dict[str, object] = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "stage": "qc_assessed_human_review_required",
        "status": status,
        "primary_completion": {
            "path": str(completion_path),
            "size_bytes": completion_path.stat().st_size,
            "sha256": sha256_file(completion_path),
        },
        "policy": {
            "path": str(policy_resolved),
            "size_bytes": policy_resolved.stat().st_size,
            "sha256": sha256_file(policy_resolved),
            "study": study,
            "metrics": metrics,
        },
        "evidence": evidence,
        "evidence_issues": evidence_issues,
        "sample_count": len(samples),
        "flagged_sample_count": len(flagged),
        "flagged_samples": flagged,
        "missing_required_sample_count": len(missing_required_samples),
        "missing_required_samples": missing_required_samples,
        "samples": [rows[sample] for sample in samples],
        "automatic_exclusions": 0,
        "human_decision_required": True,
    }
    outdir.mkdir(parents=True, mode=0o700)
    try:
        assessment = outdir / "qc_assessment.json"
        with assessment.open("x", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.chmod(assessment, 0o600)
        metric_table = outdir / "qc_sample_metrics.tsv"
        metric_ids = [str(metric["id"]) for metric in metrics]
        with metric_table.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["sample", *metric_ids, "missing_metrics", "flags"])
            for sample in samples:
                row = rows[sample]
                writer.writerow(
                    [
                        sample,
                        *["" if row["metrics"][metric_id] is None else row["metrics"][metric_id] for metric_id in metric_ids],
                        ";".join(row["missing_metrics"]),
                        ";".join(f"{flag['metric']}:{flag['kind']}" for flag in row["flags"]),
                    ]
                )
        os.chmod(metric_table, 0o600)
        decisions = outdir / "qc_decisions.tsv"
        with decisions.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["sample", "decision", "reason", "reviewer"])
            for sample in samples:
                writer.writerow([sample, "review", "", ""])
        os.chmod(decisions, 0o600)
    except Exception:
        shutil.rmtree(outdir)
        raise
    return payload


def _verified_record(path: Path, record: dict[str, object], label: str) -> None:
    if not path.is_file() or path.stat().st_size != record.get("size_bytes"):
        raise QCError(f"{label} size changed after QC assessment")
    if sha256_file(path) != record.get("sha256"):
        raise QCError(f"{label} checksum changed after QC assessment")


def _read_decisions(
    path: Path, expected_samples: list[str], flagged_samples: set[str]
) -> tuple[Path, list[dict[str, str]]]:
    try:
        resolved = path.resolve(strict=True)
        with resolved.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = ["sample", "decision", "reason", "reviewer"]
            if reader.fieldnames != required:
                raise QCError(
                    "QC decisions header must be exactly: sample, decision, reason, reviewer"
                )
            rows = [
                {key: (value or "").strip() for key, value in row.items()}
                for row in reader
                if any((value or "").strip() for value in row.values())
            ]
    except (OSError, csv.Error, UnicodeError) as exc:
        raise QCError(f"could not read QC decisions: {exc}") from exc
    by_sample: dict[str, dict[str, str]] = {}
    for row in rows:
        sample = row["sample"]
        if sample in by_sample:
            raise QCError(f"duplicate QC decision for sample: {sample!r}")
        by_sample[sample] = row
    expected = set(expected_samples)
    missing = sorted(expected - set(by_sample))
    unknown = sorted(set(by_sample) - expected)
    if missing or unknown:
        raise QCError(
            f"QC decisions do not match assessed samples: missing={missing}, unknown={unknown}"
        )
    ordered = [by_sample[sample] for sample in expected_samples]
    for row in ordered:
        sample = row["sample"]
        if row["decision"] not in {"include", "exclude"}:
            raise QCError(f"sample {sample!r} requires decision 'include' or 'exclude'")
        if not row["reviewer"]:
            raise QCError(f"sample {sample!r} requires a reviewer")
        if row["decision"] == "exclude" and not row["reason"]:
            raise QCError(f"excluded sample {sample!r} requires a reason")
        if row["decision"] == "include" and sample in flagged_samples and not row["reason"]:
            raise QCError(f"flagged included sample {sample!r} requires a review reason")
    return resolved, ordered


def _verified_plan_from_assessment(
    assessment: dict[str, object]
) -> tuple[Path, dict[str, object]]:
    completion_record = assessment.get("primary_completion")
    if not isinstance(completion_record, dict):
        raise QCError("QC assessment has no primary completion record")
    completion_path = Path(str(completion_record.get("path", ""))).resolve()
    _verified_record(completion_path, completion_record, "primary completion receipt")
    _, completion, _ = _read_completion(completion_path)
    plan_record = completion.get("run_plan")
    if not isinstance(plan_record, dict):
        raise QCError("primary completion receipt has no run-plan record")
    plan_path = Path(str(plan_record.get("path", ""))).resolve()
    _verified_record(plan_path, plan_record, "run plan")
    _, plan = _read_json(plan_path, "run plan")
    return plan_path, plan


def _verified_control(plan: dict[str, object], role: str) -> Path:
    controls = plan.get("control_files")
    if not isinstance(controls, dict) or not isinstance(controls.get(role), dict):
        raise QCError(f"run plan has no {role} control")
    record = controls[role]
    path = Path(str(record.get("path", ""))).resolve()
    _verified_record(path, record, f"run-plan {role} control")
    return path


def _read_table(path: Path) -> list[dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            return [
                {str(key).strip(): (value or "").strip() for key, value in row.items()}
                for row in csv.DictReader(handle)
            ]
    except (OSError, csv.Error, UnicodeError) as exc:
        raise QCError(f"could not read control table {path}: {exc}") from exc


def _contrast_balance(
    plan: dict[str, object], accepted: set[str]
) -> tuple[int, list[dict[str, object]]]:
    policy = plan.get("intake_policy")
    if not isinstance(policy, dict) or not isinstance(policy.get("min_replicates"), int):
        raise QCError("run plan has no bound minimum-replicate policy")
    minimum = policy["min_replicates"]
    design_path = _verified_control(plan, "design")
    contrasts_path = _verified_control(plan, "contrasts")
    design_rows = _read_table(design_path)
    design: dict[str, dict[str, str]] = {}
    for row in design_rows:
        sample = row.get("sample", "")
        if not sample or sample in design:
            raise QCError("design control must have one non-empty row per sample")
        design[sample] = row
    balances: list[dict[str, object]] = []
    for contrast in _read_table(contrasts_path):
        contrast_id = contrast.get("contrast_id", "")
        variable = contrast.get("variable", "")
        reference = contrast.get("reference", "")
        target = contrast.get("target", "")
        counts = {
            reference: sum(
                sample in accepted and row.get(variable) == reference
                for sample, row in design.items()
            ),
            target: sum(
                sample in accepted and row.get(variable) == target
                for sample, row in design.items()
            ),
        }
        balances.append(
            {
                "contrast_id": contrast_id,
                "variable": variable,
                "reference": reference,
                "target": target,
                "accepted_reference_replicates": counts[reference],
                "accepted_target_replicates": counts[target],
                "minimum_replicates": minimum,
            }
        )
        if counts[reference] < minimum or counts[target] < minimum:
            raise QCError(
                f"QC decisions leave contrast {contrast_id!r} below {minimum} replicates: "
                f"{reference}={counts[reference]}, {target}={counts[target]}"
            )
    return minimum, balances


def finalize_qc(
    *, assessment_path: Path, decisions_path: Path, outdir: Path
) -> dict[str, object]:
    """Seal explicit human decisions into an accepted-sample manifest."""

    if outdir.exists():
        raise QCError(f"refusing to overwrite existing QC acceptance directory: {outdir}")
    assessment_resolved, assessment = _read_json(assessment_path, "QC assessment")
    if assessment.get("schema_version") != 1 or assessment.get("status") != "review_required":
        raise QCError("QC assessment is not eligible for finalization")
    policy_record = assessment.get("policy")
    if not isinstance(policy_record, dict):
        raise QCError("QC assessment has no policy record")
    policy_path = Path(str(policy_record.get("path", ""))).resolve()
    _verified_record(policy_path, policy_record, "QC policy")
    evidence = assessment.get("evidence")
    if not isinstance(evidence, dict):
        raise QCError("QC assessment has no evidence inventory")
    for metric, record in evidence.items():
        if not isinstance(record, dict) or not record.get("bound_by_completion"):
            continue
        evidence_path = Path(str(record.get("path", ""))).resolve()
        _verified_record(evidence_path, record, f"QC evidence for {metric}")
    raw_samples = assessment.get("samples")
    if not isinstance(raw_samples, list) or not all(isinstance(row, dict) for row in raw_samples):
        raise QCError("QC assessment has an invalid sample list")
    samples = [str(row.get("sample", "")) for row in raw_samples]
    flagged = {str(sample) for sample in assessment.get("flagged_samples", [])}
    decisions_resolved, decisions = _read_decisions(decisions_path, samples, flagged)
    accepted = [row["sample"] for row in decisions if row["decision"] == "include"]
    excluded = [row for row in decisions if row["decision"] == "exclude"]
    plan_path, plan = _verified_plan_from_assessment(assessment)
    minimum, balances = _contrast_balance(plan, set(accepted))
    payload: dict[str, object] = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "stage": "qc_accepted_for_differential_analysis",
        "status": "accepted",
        "qc_assessment": {
            "path": str(assessment_resolved),
            "size_bytes": assessment_resolved.stat().st_size,
            "sha256": sha256_file(assessment_resolved),
        },
        "decisions": {
            "path": str(decisions_resolved),
            "size_bytes": decisions_resolved.stat().st_size,
            "sha256": sha256_file(decisions_resolved),
        },
        "run_plan": {"path": str(plan_path), "sha256": sha256_file(plan_path)},
        "minimum_replicates": minimum,
        "contrast_balance": balances,
        "n_assessed": len(samples),
        "n_accepted": len(accepted),
        "n_excluded": len(excluded),
        "accepted_samples": accepted,
        "excluded_samples": excluded,
        "review_decisions": decisions,
        "automatic_exclusions": 0,
        "human_review_complete": True,
    }
    outdir.mkdir(parents=True, mode=0o700)
    try:
        receipt = outdir / "qc_acceptance.json"
        with receipt.open("x", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.chmod(receipt, 0o600)
        manifest = outdir / "accepted_samples.tsv"
        with manifest.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["sample"])
            writer.writerows([sample] for sample in accepted)
        os.chmod(manifest, 0o600)
    except Exception:
        shutil.rmtree(outdir)
        raise
    return payload
