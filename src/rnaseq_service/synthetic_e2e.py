"""Exercise every service gate with clearly labeled synthetic artifacts."""

from __future__ import annotations

import csv
import json
import math
import os
import shutil
import zipfile
from pathlib import Path
from typing import Callable

from .differential import create_differential_plan, create_differential_receipt
from .inputs import create_input_manifest
from .plan import create_rnaseq_plan, sha256_file
from .preflight import run_preflight
from .primary import create_primary_receipt
from .qc import evaluate_qc, finalize_qc
from .synthetic import GENES, PAIR_COUNTS, create_synthetic_study


ContractProgress = Callable[[str, int, int], None]
STAGES = (
    "synthetic_study",
    "intake_preflight",
    "input_manifest",
    "primary_plan",
    "primary_completion",
    "qc_acceptance",
    "differential_plan",
    "differential_completion",
)
PNG_1X1 = bytes.fromhex(
    "89504e470d0a1a0a0000000d49484452000000010000000108060000001f15c489"
    "0000000d49444154789c6360000000020001e221bc330000000049454e44ae426082"
)


class SyntheticContractError(ValueError):
    """Raised when the synthetic end-to-end contract is not satisfied."""


def _advance(progress: ContractProgress | None, stage: str) -> None:
    if progress is not None:
        progress(stage, STAGES.index(stage) + 1, len(STAGES))


def _write_tsv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def _matrix_rows() -> list[list[object]]:
    samples = list(PAIR_COUNTS)
    return [
        [gene, *[PAIR_COUNTS[sample][index] for sample in samples]]
        for index, gene in enumerate(GENES)
    ]


def _materialize_primary_outputs(plan: dict[str, object]) -> None:
    run_root = Path(str(plan["execution"]["outdir"]))
    samples = list(PAIR_COUNTS)
    execution = run_root / "execution"
    execution.mkdir(parents=True)
    for name in ("report.html", "timeline.html", "dag.html"):
        (execution / name).write_text(f"<html>synthetic {name}</html>\n", encoding="utf-8")
    (execution / "trace.tsv").write_text(
        "name\tstatus\tduration\trealtime\t%cpu\tpeak_rss\n"
        "FASTQC\tCOMPLETED\t1s\t1s\t80%\t20 MB\n"
        "STAR_ALIGN\tCACHED\t2s\t2s\t90%\t30 MB\n"
        "SALMON\tCOMPLETED\t1s\t1s\t75%\t20 MB\n",
        encoding="utf-8",
    )
    rnaseq = run_root / "rnaseq"
    info = rnaseq / "pipeline_info"
    info.mkdir(parents=True)
    (info / "nf_core_rnaseq_software_mqc_versions.yml").write_text(
        "synthetic-contract: 1\n", encoding="utf-8"
    )
    samplesheet = Path(str(plan["control_files"]["samplesheet"]["path"]))
    (info / "params.json").write_text(
        json.dumps(
            {
                "input": str(samplesheet),
                "outdir": str(rnaseq.resolve()),
                "genome": "GRCh38",
                "synthetic_contract_artifact": True,
            }
        ),
        encoding="utf-8",
    )
    (info / "samplesheet.valid.csv").write_text(
        samplesheet.read_text(encoding="utf-8"), encoding="utf-8"
    )
    quantification = rnaseq / "star_salmon"
    matrix_header = ["gene_id", *samples]
    rows = _matrix_rows()
    _write_tsv(quantification / "salmon.merged.gene_counts.tsv", matrix_header, rows)
    _write_tsv(
        quantification / "salmon.merged.gene_lengths.tsv",
        matrix_header,
        [[gene, *([500] * len(samples))] for gene in GENES],
    )
    _write_tsv(quantification / "salmon.merged.gene_tpm.tsv", matrix_header, rows)
    multiqc = rnaseq / "multiqc" / "star_salmon"
    data = multiqc / "multiqc_report_data"
    data.mkdir(parents=True)
    (multiqc / "multiqc_report.html").write_text(
        "<html>synthetic MultiQC contract artifact</html>\n", encoding="utf-8"
    )
    (data / "multiqc_data.json").write_text(
        json.dumps({"synthetic_contract_artifact": True}) + "\n", encoding="utf-8"
    )
    _write_tsv(
        data / "multiqc_general_stats.txt",
        ["Sample", "mapped_percent"],
        [[sample, 95] for sample in samples],
    )


def _write_qc_controls(root: Path) -> tuple[Path, Path]:
    policy = root / "private" / "qc" / "policy.toml"
    policy.parent.mkdir(parents=True)
    policy.write_text(
        "schema_version = 1\n"
        "[study]\n"
        "robust_z_threshold = 4.5\n"
        "min_samples_for_outliers = 4\n"
        "[[metrics]]\n"
        'id = "mapping_percent"\n'
        'kind = "numeric"\n'
        'source_file = "multiqc_general_stats.txt"\n'
        'sample_column = "Sample"\n'
        'value_column = "mapped_percent"\n'
        "required = true\n"
        "lower = 70.0\n"
        "upper = 100.0\n"
        "outlier = true\n",
        encoding="utf-8",
    )
    decisions = root / "private" / "qc" / "reviewed_decisions.tsv"
    _write_tsv(
        decisions,
        ["sample", "decision", "reason", "reviewer"],
        [[sample, "include", "", "synthetic-contract-validator"] for sample in PAIR_COUNTS],
    )
    return policy, decisions


def _log2_fold_changes() -> dict[str, float]:
    controls = [sample for sample in PAIR_COUNTS if sample.startswith("control_")]
    treated = [sample for sample in PAIR_COUNTS if sample.startswith("treated_")]
    changes: dict[str, float] = {}
    for index, gene in enumerate(GENES):
        reference_mean = sum(PAIR_COUNTS[sample][index] for sample in controls) / len(controls)
        target_mean = sum(PAIR_COUNTS[sample][index] for sample in treated) / len(treated)
        changes[gene] = math.log2(target_mean / reference_mean)
    return changes


def _materialize_differential_outputs(plan: dict[str, object]) -> None:
    run_root = Path(str(plan["execution"]["outdir"]))
    pipeline = run_root / "differentialabundance"
    execution = run_root / "execution"
    execution.mkdir(parents=True)
    for name in ("report.html", "timeline.html", "dag.html"):
        (execution / name).write_text(f"<html>synthetic {name}</html>\n", encoding="utf-8")
    (execution / "trace.tsv").write_text(
        "name\tstatus\tduration\trealtime\t%cpu\tpeak_rss\n"
        "DESEQ2_NORM\tCOMPLETED\t1s\t1s\t80%\t20 MB\n"
        "DESEQ2_DIFFERENTIAL\tCOMPLETED\t1s\t1s\t80%\t20 MB\n",
        encoding="utf-8",
    )
    controls = plan["control_files"]
    info = pipeline / "pipeline_info"
    info.mkdir(parents=True)
    (info / "nf_core_differentialabundance_software_versions.yml").write_text(
        "synthetic-contract: 1\n", encoding="utf-8"
    )
    (info / "params.json").write_text(
        json.dumps(
            {
                "input": controls["observations"]["path"],
                "matrix": controls["gene_counts"]["path"],
                "feature_length_matrix": controls["gene_lengths"]["path"],
                "contrasts": controls["contrasts"]["path"],
                "outdir": str(pipeline.resolve()),
                "study_name": "synthetice2e",
                "study_type": "rnaseq",
                "observations_id_col": "sample",
                "features_id_col": "gene_id",
                "genome": "GRCh38",
                "differential_method": "deseq2",
                "deseq2_vs_method": "vst",
                "functional_method": "none",
                "skip_reports": False,
                "synthetic_contract_artifact": True,
            }
        ),
        encoding="utf-8",
    )
    (info / "samplesheet.valid.csv").write_text(
        Path(str(controls["observations"]["path"])).read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    report = pipeline / "report" / "synthetice2e"
    report.mkdir(parents=True)
    html = report / "synthetice2e_differentialabundance_report.html"
    html.write_text("<html>synthetic differential report</html>\n", encoding="utf-8")
    with zipfile.ZipFile(report / "synthetice2e.zip", "x") as archive:
        archive.writestr("README.txt", "Synthetic contract artifact; no scientific execution.\n")

    samples = list(PAIR_COUNTS)
    matrix_header = ["gene_id", *samples]
    rows = _matrix_rows()
    processed = pipeline / "tables" / "processed_abundance" / "synthetice2e"
    _write_tsv(processed / "all.normalised_counts.tsv", matrix_header, rows)
    _write_tsv(
        processed / "all.vst.tsv",
        matrix_header,
        [[row[0], *[math.log2(float(value) + 1) for value in row[1:]]] for row in rows],
    )
    changes = _log2_fold_changes()
    result_header = ["gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    result_rows = [
        [gene, 10, changes[gene], 0.25, changes[gene] / 0.25, 0.001, 0.01]
        for gene in GENES
    ]
    differential = pipeline / "tables" / "differential" / "synthetice2e"
    _write_tsv(
        differential / "treated_vs_control.deseq2.results.tsv",
        result_header,
        result_rows,
    )
    _write_tsv(
        differential / "treated_vs_control.deseq2.results_filtered.tsv",
        result_header,
        [row for row in result_rows if abs(float(row[2])) >= 2],
    )
    _write_tsv(
        differential / "treated_vs_control_deseq2.annotated.tsv",
        result_header,
        result_rows,
    )
    volcano = (
        pipeline
        / "plots"
        / "differential"
        / "synthetice2e"
        / "treated_vs_control"
        / "png"
    )
    volcano.mkdir(parents=True)
    (volcano / "volcano.png").write_bytes(PNG_1X1)
    qc = pipeline / "plots" / "qc" / "synthetice2e"
    qc.mkdir(parents=True)
    (qc / "treated_vs_control.deseq2.dispersion.png").write_bytes(PNG_1X1)
    exploratory = pipeline / "plots" / "exploratory" / "synthetice2e" / "condition" / "png"
    exploratory.mkdir(parents=True)
    (exploratory / "pca2d.png").write_bytes(PNG_1X1)
    (exploratory / "mad_correlation.png").write_bytes(PNG_1X1)


def _verify_direction(completion: dict[str, object]) -> dict[str, object]:
    record = next(
        item for item in completion["artifacts"] if item["role"] == "deseq2_results"
    )
    with Path(str(record["path"])).open("r", encoding="utf-8", newline="") as handle:
        values = {
            row["gene_id"]: float(row["log2FoldChange"])
            for row in csv.DictReader(handle, delimiter="\t")
        }
    checks = {
        "gene_treated_marker_positive": values["gene_treated_marker"] > 2,
        "gene_control_marker_negative": values["gene_control_marker"] < -2,
        "gene_housekeeping_near_zero": abs(values["gene_housekeeping"]) < 0.1,
        "gene_neutral_near_zero": abs(values["gene_neutral"]) < 0.1,
    }
    if not all(checks.values()):
        raise SyntheticContractError(f"synthetic direction validation failed: {checks}")
    return {"passed": True, "checks": checks, "log2_fold_changes": values}


def run_contract_validation(
    *,
    output_dir: Path,
    workflow_lock: Path,
    progress: ContractProgress | None = None,
) -> dict[str, object]:
    """Run the complete wrapper contract without executing scientific workflows."""

    if output_dir.exists():
        raise SyntheticContractError(f"refusing to overwrite contract validation: {output_dir}")
    root = output_dir.resolve()
    try:
        root.mkdir(parents=True)
        study = root / "study"
        create_synthetic_study(output_dir=study)
        _advance(progress, "synthetic_study")
        samplesheet = study / "intake" / "samplesheet.csv"
        design = study / "intake" / "design.csv"
        contrasts = study / "intake" / "contrasts.csv"
        preflight = run_preflight(samplesheet, design, contrasts, check_files=True)
        if not preflight.valid:
            raise SyntheticContractError("synthetic intake preflight failed")
        _advance(progress, "intake_preflight")

        private = root / "private"
        input_manifest = private / "inputs" / "fastq-manifest.json"
        create_input_manifest(samplesheet=samplesheet, output=input_manifest)
        _advance(progress, "input_manifest")
        primary_plan_path = private / "plans" / "primary.json"
        primary_plan = create_rnaseq_plan(
            samplesheet=samplesheet,
            design=design,
            contrasts=contrasts,
            workflow_lock=workflow_lock,
            output=primary_plan_path,
            outdir=root / "results" / "primary",
            workdir=root / "work" / "primary",
            network_mode="direct",
            container_engine="apptainer",
            executor="local",
            input_manifest=input_manifest,
            genome="GRCh38",
        )
        _advance(progress, "primary_plan")
        _materialize_primary_outputs(primary_plan)
        primary_completion_path = private / "receipts" / "primary-completion.json"
        create_primary_receipt(run_plan=primary_plan_path, output=primary_completion_path)
        _advance(progress, "primary_completion")

        policy, decisions = _write_qc_controls(root)
        assessment_dir = private / "qc" / "assessment"
        assessment = evaluate_qc(
            completion_receipt=primary_completion_path,
            policy_path=policy,
            outdir=assessment_dir,
        )
        if assessment["status"] != "review_required" or assessment["flagged_samples"]:
            raise SyntheticContractError("synthetic QC did not reach clean human review")
        acceptance_dir = private / "qc" / "accepted"
        acceptance = finalize_qc(
            assessment_path=assessment_dir / "qc_assessment.json",
            decisions_path=decisions,
            outdir=acceptance_dir,
        )
        _advance(progress, "qc_acceptance")

        differential_plan_path = private / "plans" / "differential.json"
        differential_plan = create_differential_plan(
            qc_acceptance=acceptance_dir / "qc_acceptance.json",
            handoff_dir=private / "differential" / "inputs",
            output=differential_plan_path,
            outdir=root / "results" / "differential",
            workdir=root / "work" / "differential",
            study_name="synthetice2e",
            network_mode="direct",
            container_engine="apptainer",
            executor="local",
            seed=20260909,
        )
        _advance(progress, "differential_plan")
        _materialize_differential_outputs(differential_plan)
        differential_completion_path = private / "receipts" / "differential-completion.json"
        differential_completion = create_differential_receipt(
            run_plan=differential_plan_path,
            output=differential_completion_path,
        )
        direction = _verify_direction(differential_completion)
        _advance(progress, "differential_completion")

        stage_paths = {
            "input_manifest": input_manifest,
            "primary_plan": primary_plan_path,
            "primary_completion": primary_completion_path,
            "qc_assessment": assessment_dir / "qc_assessment.json",
            "qc_acceptance": acceptance_dir / "qc_acceptance.json",
            "differential_plan": differential_plan_path,
            "differential_completion": differential_completion_path,
        }
        receipt: dict[str, object] = {
            "schema_version": 1,
            "stage": "synthetic_contract_validation_complete",
            "scientific_execution_performed": False,
            "synthetic_upstream_artifacts": True,
            "validates": "service_orchestration_and_artifact_contracts_only",
            "n_samples": len(PAIR_COUNTS),
            "n_genes": len(GENES),
            "n_contrasts": 1,
            "automatic_exclusions": acceptance["automatic_exclusions"],
            "direction_validation": direction,
            "artifacts": {
                name: {
                    "path": str(path),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
                for name, path in stage_paths.items()
            },
            "contains_client_data": False,
            "contains_secrets": False,
        }
        receipt_path = root / "synthetic_contract_receipt.json"
        descriptor = os.open(receipt_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(receipt, handle, indent=2, sort_keys=True)
            handle.write("\n")
        return receipt
    except Exception:
        shutil.rmtree(root, ignore_errors=True)
        raise

