"""Build checksum-bound nf-core/differentialabundance handoff plans."""

from __future__ import annotations

import csv
import json
import math
import os
import re
import shlex
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

from .bundle import BundleError, verify_bundle
from .plan import (
    CONTAINER_ENGINES,
    EXECUTORS,
    NETWORK_MODES,
    PROXY_VARIABLES,
    sha256_file,
)
from .workflow_lock import WorkflowLockError, load_workflow_lock


MatrixProgress = Callable[[str, int, int], None]
STUDY_NAME = re.compile(r"^[A-Za-z0-9]+$")


class DifferentialHandoffError(ValueError):
    """Raised when an accepted differential handoff cannot be trusted."""


def _read_json(path: Path, label: str) -> tuple[Path, dict[str, object]]:
    try:
        resolved = path.resolve(strict=True)
        payload = json.loads(resolved.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise DifferentialHandoffError(f"could not read {label}: {exc}") from exc
    if not isinstance(payload, dict):
        raise DifferentialHandoffError(f"{label} must contain a JSON object")
    return resolved, payload


def _verify_record(record: dict[str, object], label: str) -> Path:
    try:
        path = Path(str(record["path"])).resolve(strict=True)
    except (KeyError, OSError) as exc:
        raise DifferentialHandoffError(f"{label} is unavailable") from exc
    if not path.is_file() or path.stat().st_size != record.get("size_bytes"):
        raise DifferentialHandoffError(f"{label} size changed")
    if sha256_file(path) != record.get("sha256"):
        raise DifferentialHandoffError(f"{label} checksum changed")
    return path


def _control_record(path: Path) -> dict[str, object]:
    resolved = path.resolve(strict=True)
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256_file(resolved),
    }


def _accepted_chain(
    acceptance_path: Path,
) -> tuple[
    Path,
    dict[str, object],
    Path,
    dict[str, object],
    dict[str, dict[str, object]],
]:
    acceptance_resolved, acceptance = _read_json(acceptance_path, "QC acceptance receipt")
    if acceptance.get("schema_version") != 1:
        raise DifferentialHandoffError("QC acceptance receipt must use schema version 1")
    if (
        acceptance.get("stage") != "qc_accepted_for_differential_analysis"
        or acceptance.get("status") != "accepted"
        or acceptance.get("human_review_complete") is not True
        or acceptance.get("automatic_exclusions") != 0
    ):
        raise DifferentialHandoffError("QC acceptance receipt is not eligible for handoff")
    accepted = acceptance.get("accepted_samples")
    decisions = acceptance.get("review_decisions")
    if (
        not isinstance(accepted, list)
        or not accepted
        or not all(isinstance(sample, str) and sample for sample in accepted)
        or len(accepted) != len(set(accepted))
        or not isinstance(decisions, list)
    ):
        raise DifferentialHandoffError("QC acceptance receipt has an invalid sample set")
    included_from_decisions = [
        row.get("sample")
        for row in decisions
        if isinstance(row, dict) and row.get("decision") == "include"
    ]
    if included_from_decisions != accepted:
        raise DifferentialHandoffError("accepted samples conflict with reviewed decisions")
    assessment_record = acceptance.get("qc_assessment")
    if not isinstance(assessment_record, dict):
        raise DifferentialHandoffError("QC acceptance receipt has no assessment record")
    assessment_path = _verify_record(assessment_record, "QC assessment")
    _, assessment = _read_json(assessment_path, "QC assessment")
    completion_record = assessment.get("primary_completion")
    if not isinstance(completion_record, dict):
        raise DifferentialHandoffError("QC assessment has no primary completion record")
    completion_path = _verify_record(completion_record, "primary completion receipt")
    _, completion = _read_json(completion_path, "primary completion receipt")
    if completion.get("schema_version") != 2 or completion.get("stage") != "pipeline_complete_qc_pending":
        raise DifferentialHandoffError("primary completion receipt is not supported")
    artifacts_raw = completion.get("artifacts")
    if not isinstance(artifacts_raw, list):
        raise DifferentialHandoffError("primary completion receipt has no artifacts")
    artifacts: dict[str, dict[str, object]] = {}
    for record in artifacts_raw:
        if not isinstance(record, dict) or not isinstance(record.get("role"), str):
            raise DifferentialHandoffError("primary completion receipt has an invalid artifact")
        role = record["role"]
        if role in artifacts:
            raise DifferentialHandoffError(f"duplicate primary artifact role: {role}")
        artifacts[role] = record
    for role in ("gene_counts", "gene_lengths"):
        if role not in artifacts:
            raise DifferentialHandoffError(
                f"primary completion receipt lacks {role}; create a current completion receipt"
            )
    return acceptance_resolved, acceptance, completion_path, completion, artifacts


def _source_plan(completion: dict[str, object]) -> tuple[Path, dict[str, object]]:
    record = completion.get("run_plan")
    if not isinstance(record, dict):
        raise DifferentialHandoffError("primary completion receipt has no source run plan")
    path = _verify_record(record, "source RNA-seq run plan")
    _, plan = _read_json(path, "source RNA-seq run plan")
    if plan.get("schema_version") != 1 or plan.get("stage") != "planned_not_executed":
        raise DifferentialHandoffError("source RNA-seq run plan is not supported")
    reference = plan.get("reference")
    if (
        not isinstance(reference, dict)
        or reference.get("mode") != "igenomes"
        or not isinstance(reference.get("genome"), str)
    ):
        raise DifferentialHandoffError("source run plan has no reviewed iGenomes reference")
    return path, plan


def _source_control(plan: dict[str, object], role: str) -> Path:
    controls = plan.get("control_files")
    if not isinstance(controls, dict) or not isinstance(controls.get(role), dict):
        raise DifferentialHandoffError(f"source run plan has no {role} control")
    return _verify_record(controls[role], f"source {role} control")


def _subset_matrix(
    *,
    source: Path,
    output: Path,
    samples: list[str],
    label: str,
    progress: MatrixProgress | None,
) -> list[str]:
    total = source.stat().st_size
    completed = 0
    feature_ids: list[str] = []
    seen_features: set[str] = set()
    try:
        with source.open("r", encoding="utf-8", newline="") as src, output.open(
            "x", encoding="utf-8", newline=""
        ) as dst:
            header_line = src.readline()
            completed += len(header_line.encode("utf-8"))
            header = header_line.rstrip("\r\n").split("\t")
            if not header or header[0] != "gene_id" or len(header) != len(set(header)):
                raise DifferentialHandoffError(
                    f"{label} must start with gene_id and have unique columns"
                )
            missing = sorted(set(samples) - set(header))
            if missing:
                raise DifferentialHandoffError(
                    f"{label} is missing accepted samples: {', '.join(missing[:5])}"
                )
            indices = [0, *[header.index(sample) for sample in samples]]
            writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
            writer.writerow([header[index] for index in indices])
            for line_number, line in enumerate(src, start=2):
                completed += len(line.encode("utf-8"))
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != len(header):
                    raise DifferentialHandoffError(
                        f"{label} row {line_number} has {len(fields)} columns; expected {len(header)}"
                    )
                feature = fields[0].strip()
                if not feature or feature in seen_features:
                    raise DifferentialHandoffError(
                        f"{label} has an empty or duplicate gene_id at row {line_number}"
                    )
                seen_features.add(feature)
                for index in indices[1:]:
                    try:
                        value = float(fields[index])
                    except ValueError as exc:
                        raise DifferentialHandoffError(
                            f"{label} has a non-numeric value at row {line_number}"
                        ) from exc
                    if not math.isfinite(value) or value < 0:
                        raise DifferentialHandoffError(
                            f"{label} has a negative or non-finite value at row {line_number}"
                        )
                writer.writerow([fields[index] for index in indices])
                feature_ids.append(feature)
                if progress is not None:
                    progress(label, min(completed, total), total)
    except DifferentialHandoffError:
        raise
    except (OSError, UnicodeError, csv.Error) as exc:
        raise DifferentialHandoffError(f"could not subset {label}: {exc}") from exc
    if not feature_ids:
        raise DifferentialHandoffError(f"{label} contains no feature rows")
    if progress is not None:
        progress(label, total, total)
    os.chmod(output, 0o600)
    return feature_ids


def _read_csv(path: Path, label: str) -> tuple[list[str], list[dict[str, str]]]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                raise DifferentialHandoffError(f"{label} has no header")
            headers = [header.strip() for header in reader.fieldnames]
            rows = [
                {str(key).strip(): (value or "").strip() for key, value in row.items()}
                for row in reader
                if any((value or "").strip() for value in row.values())
            ]
    except DifferentialHandoffError:
        raise
    except (OSError, UnicodeError, csv.Error) as exc:
        raise DifferentialHandoffError(f"could not read {label}: {exc}") from exc
    return headers, rows


def _write_observations(source: Path, output: Path, samples: list[str]) -> list[str]:
    headers, rows = _read_csv(source, "design")
    if not headers or headers[0] != "sample" or len(headers) != len(set(headers)):
        raise DifferentialHandoffError("design must start with sample and have unique columns")
    by_sample: dict[str, dict[str, str]] = {}
    for row in rows:
        sample = row.get("sample", "")
        if not sample or sample in by_sample:
            raise DifferentialHandoffError("design must have one non-empty row per sample")
        by_sample[sample] = row
    missing = sorted(set(samples) - set(by_sample))
    if missing:
        raise DifferentialHandoffError(
            f"design is missing accepted samples: {', '.join(missing[:5])}"
        )
    with output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({header: by_sample[sample].get(header, "") for header in headers} for sample in samples)
    os.chmod(output, 0o600)
    return headers


def _write_contrasts(source: Path, output: Path, observation_headers: list[str]) -> int:
    headers, rows = _read_csv(source, "contrasts")
    required = {"contrast_id", "variable", "reference", "target"}
    if not required.issubset(headers) or not rows:
        raise DifferentialHandoffError("contrasts lack required columns or data rows")
    seen: set[str] = set()
    output_rows: list[dict[str, str]] = []
    for row in rows:
        contrast_id = row.get("contrast_id", "")
        variable = row.get("variable", "")
        blocking = row.get("blocking", "")
        if not contrast_id or contrast_id in seen:
            raise DifferentialHandoffError("contrast IDs must be unique and non-empty")
        seen.add(contrast_id)
        if variable not in observation_headers:
            raise DifferentialHandoffError(
                f"contrast {contrast_id!r} uses unknown variable {variable!r}"
            )
        blocking_factors = [value.strip() for value in blocking.split(";") if value.strip()]
        unknown = sorted(set(blocking_factors) - set(observation_headers))
        if unknown:
            raise DifferentialHandoffError(
                f"contrast {contrast_id!r} uses unknown blocking factors: {', '.join(unknown)}"
            )
        output_rows.append(
            {
                "id": contrast_id,
                "variable": variable,
                "reference": row.get("reference", ""),
                "target": row.get("target", ""),
                "blocking": ";".join(blocking_factors),
            }
        )
    output_headers = ["id", "variable", "reference", "target", "blocking"]
    with output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_headers, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    os.chmod(output, 0o600)
    return len(output_rows)


def _build_argv(
    *,
    workflow_name: str,
    revision: str | None,
    container_engine: str,
    infrastructure_config: Path | None,
    workdir: Path,
    execution_root: Path,
    pipeline_outdir: Path,
    observations: Path,
    counts: Path,
    lengths: Path,
    contrasts: Path,
    genome: str,
    study_name: str,
    seed: int,
) -> list[str]:
    argv = ["nextflow", "run", workflow_name]
    if revision is not None:
        argv.extend(["-r", revision])
    argv.extend(
        [
            "-profile",
            f"rnaseq,{container_engine}",
            "-work-dir",
            str(workdir),
        ]
    )
    if infrastructure_config is not None:
        argv.extend(["-c", str(infrastructure_config)])
    argv.extend(
        [
            "-resume",
            "-with-report",
            str(execution_root / "report.html"),
            "-with-trace",
            str(execution_root / "trace.tsv"),
            "-with-timeline",
            str(execution_root / "timeline.html"),
            "-with-dag",
            str(execution_root / "dag.html"),
            "--input",
            str(observations),
            "--matrix",
            str(counts),
            "--feature_length_matrix",
            str(lengths),
            "--contrasts",
            str(contrasts),
            "--genome",
            genome,
            "--study_name",
            study_name,
            "--study_type",
            "rnaseq",
            "--observations_id_col",
            "sample",
            "--features_id_col",
            "gene_id",
            "--seed",
            str(seed),
            "--outdir",
            str(pipeline_outdir),
        ]
    )
    return argv


def create_differential_plan(
    *,
    qc_acceptance: Path,
    handoff_dir: Path,
    output: Path,
    outdir: Path,
    workdir: Path,
    study_name: str,
    network_mode: str,
    container_engine: str,
    executor: str,
    seed: int,
    infrastructure_config: Path | None = None,
    offline_manifest: Path | None = None,
    progress: MatrixProgress | None = None,
) -> dict[str, object]:
    """Subset accepted inputs and write a non-executing differential plan."""

    if handoff_dir.exists():
        raise DifferentialHandoffError(f"refusing to overwrite handoff directory: {handoff_dir}")
    if output.exists():
        raise DifferentialHandoffError(f"refusing to overwrite differential plan: {output}")
    if STUDY_NAME.fullmatch(study_name) is None:
        raise DifferentialHandoffError(
            "study name must contain only letters and numbers"
        )
    if network_mode not in NETWORK_MODES:
        raise DifferentialHandoffError(f"unsupported network mode: {network_mode!r}")
    if container_engine not in CONTAINER_ENGINES:
        raise DifferentialHandoffError(f"unsupported container engine: {container_engine!r}")
    if executor not in EXECUTORS:
        raise DifferentialHandoffError(f"unsupported executor: {executor!r}")
    if executor == "slurm" and infrastructure_config is None:
        raise DifferentialHandoffError("--infrastructure-config is required for slurm")
    if executor == "local" and infrastructure_config is not None:
        raise DifferentialHandoffError("--infrastructure-config is only accepted with slurm")
    if network_mode == "offline" and offline_manifest is None:
        raise DifferentialHandoffError("--offline-manifest is required in offline mode")
    if network_mode != "offline" and offline_manifest is not None:
        raise DifferentialHandoffError("--offline-manifest is only accepted in offline mode")
    if network_mode == "offline" and container_engine == "docker":
        raise DifferentialHandoffError("offline bundles support apptainer/singularity only")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise DifferentialHandoffError("seed must be a non-negative integer")

    (
        acceptance_path,
        acceptance,
        completion_path,
        completion,
        artifacts,
    ) = _accepted_chain(qc_acceptance)
    source_plan_path, source_plan = _source_plan(completion)
    workflow_lock_path = _source_control(source_plan, "workflow_lock")
    design_path = _source_control(source_plan, "design")
    contrasts_path = _source_control(source_plan, "contrasts")
    try:
        workflow_lock = load_workflow_lock(workflow_lock_path)
    except WorkflowLockError as exc:
        raise DifferentialHandoffError(str(exc)) from exc
    counts_source = _verify_record(artifacts["gene_counts"], "gene-count matrix")
    lengths_source = _verify_record(artifacts["gene_lengths"], "gene-length matrix")
    infrastructure_resolved = (
        infrastructure_config.resolve(strict=True) if infrastructure_config is not None else None
    )
    if infrastructure_resolved is not None:
        source_record = _control_record(infrastructure_resolved)
    else:
        source_record = None
    workflow_name = workflow_lock.differential.name
    workflow_revision: str | None = workflow_lock.differential.revision
    required_environment: dict[str, str] = {}
    offline_resolved: Path | None = None
    if offline_manifest is not None:
        offline_resolved = offline_manifest.resolve(strict=True)
        try:
            bundle = verify_bundle(offline_resolved, expected_workflow_lock=workflow_lock_path)
        except (BundleError, OSError) as exc:
            raise DifferentialHandoffError(f"offline bundle verification failed: {exc}") from exc
        components = bundle["components"]
        workflow_name = str((offline_resolved.parent / components["differential_workflow"]).resolve())
        workflow_revision = None
        required_environment = {
            "NXF_OFFLINE": "true",
            "NXF_SINGULARITY_CACHEDIR": str(
                (offline_resolved.parent / components["container_root"]).resolve()
            ),
            "NXF_PLUGINS_DIR": str(
                (offline_resolved.parent / components["plugin_root"]).resolve()
            ),
        }
    samples = list(acceptance["accepted_samples"])
    handoff_dir.mkdir(parents=True, mode=0o700)
    observations = handoff_dir / "observations.accepted.tsv"
    counts = handoff_dir / "gene_counts.accepted.tsv"
    lengths = handoff_dir / "gene_lengths.accepted.tsv"
    converted_contrasts = handoff_dir / "contrasts.accepted.tsv"
    try:
        observation_headers = _write_observations(design_path, observations, samples)
        count_features = _subset_matrix(
            source=counts_source,
            output=counts,
            samples=samples,
            label="gene_counts",
            progress=progress,
        )
        length_features = _subset_matrix(
            source=lengths_source,
            output=lengths,
            samples=samples,
            label="gene_lengths",
            progress=progress,
        )
        if count_features != length_features:
            raise DifferentialHandoffError(
                "gene-count and gene-length matrices have different feature IDs or row order"
            )
        n_contrasts = _write_contrasts(
            contrasts_path, converted_contrasts, observation_headers
        )
        outdir_resolved = outdir.resolve()
        workdir_resolved = workdir.resolve()
        execution_root = outdir_resolved / "execution"
        pipeline_outdir = outdir_resolved / "differentialabundance"
        argv = _build_argv(
            workflow_name=workflow_name,
            revision=workflow_revision,
            container_engine=container_engine,
            infrastructure_config=infrastructure_resolved,
            workdir=workdir_resolved,
            execution_root=execution_root,
            pipeline_outdir=pipeline_outdir,
            observations=observations.resolve(),
            counts=counts.resolve(),
            lengths=lengths.resolve(),
            contrasts=converted_contrasts.resolve(),
            genome=str(source_plan["reference"]["genome"]),
            study_name=study_name,
            seed=seed,
        )
        controls = {
            "qc_acceptance": _control_record(acceptance_path),
            "workflow_lock": _control_record(workflow_lock_path),
            "observations": _control_record(observations),
            "gene_counts": _control_record(counts),
            "gene_lengths": _control_record(lengths),
            "contrasts": _control_record(converted_contrasts),
        }
        if source_record is not None:
            controls["infrastructure_config"] = source_record
        if offline_resolved is not None:
            controls["offline_manifest"] = _control_record(offline_resolved)
        plan: dict[str, object] = {
            "schema_version": 1,
            "created_at": datetime.now(timezone.utc).isoformat(),
            "stage": "planned_not_executed",
            "analysis": "differential_abundance",
            "workflow": {
                "name": workflow_lock.differential.name,
                "revision": workflow_lock.differential.revision,
            },
            "runtime": {
                "nextflow_version": workflow_lock.runtime.nextflow_version,
                "nf_core_tools_version": workflow_lock.runtime.nf_core_tools_version,
            },
            "reference": source_plan["reference"],
            "source_primary": {
                "run_plan": _control_record(source_plan_path),
                "completion_receipt": _control_record(completion_path),
            },
            "control_files": controls,
            "accepted_sample_count": len(samples),
            "accepted_samples": samples,
            "feature_count": len(count_features),
            "contrast_count": n_contrasts,
            "seed": seed,
            "execution": {
                "executor": executor,
                "container_engine": container_engine,
                "network_mode": network_mode,
                "outdir": str(outdir_resolved),
                "workdir": str(workdir_resolved),
                "proxy_environment_present": [
                    name for name in PROXY_VARIABLES if os.environ.get(name)
                ],
                "proxy_values_recorded": False,
                "required_environment": required_environment,
            },
            "command_argv": argv,
            "command_preview": shlex.join(argv),
            "contains_client_results": True,
            "contains_secrets": False,
        }
        output.parent.mkdir(parents=True, exist_ok=True)
        descriptor = os.open(output, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(plan, handle, indent=2, sort_keys=True)
            handle.write("\n")
    except Exception:
        shutil.rmtree(handoff_dir)
        output.unlink(missing_ok=True)
        raise
    return plan
