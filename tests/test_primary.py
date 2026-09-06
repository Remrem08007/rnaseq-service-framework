from __future__ import annotations

import csv
import json
import stat
from pathlib import Path

import pytest

from rnaseq_service.plan import create_rnaseq_plan
from rnaseq_service.primary import PrimaryDeliveryError, create_primary_receipt


ROOT = Path(__file__).parents[1]


def write_csv(path: Path, headers: list[str], rows: list[list[str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        writer.writerows(rows)
    return path


def completed_fixture(tmp_path: Path) -> tuple[Path, Path]:
    reads = tmp_path / "reads"
    reads.mkdir()
    sample_rows = []
    design_rows = []
    samples = []
    for condition in ("control", "treated"):
        for replicate in (1, 2):
            sample = f"{condition}_{replicate}"
            samples.append(sample)
            r1 = reads / f"{sample}_R1.fastq.gz"
            r2 = reads / f"{sample}_R2.fastq.gz"
            r1.write_bytes(b"synthetic-r1")
            r2.write_bytes(b"synthetic-r2")
            sample_rows.append([sample, str(r1), str(r2), "auto"])
            design_rows.append([sample, condition])
    samplesheet = write_csv(
        tmp_path / "intake" / "samplesheet.csv",
        ["sample", "fastq_1", "fastq_2", "strandedness"],
        sample_rows,
    )
    design = write_csv(
        tmp_path / "intake" / "design.csv",
        ["sample", "condition"],
        design_rows,
    )
    contrasts = write_csv(
        tmp_path / "intake" / "contrasts.csv",
        ["contrast_id", "variable", "reference", "target"],
        [["treated_vs_control", "condition", "control", "treated"]],
    )
    plan_path = tmp_path / "private" / "run-plan.json"
    run_root = tmp_path / "results"
    create_rnaseq_plan(
        samplesheet=samplesheet,
        design=design,
        contrasts=contrasts,
        workflow_lock=ROOT / "config" / "workflows.toml",
        output=plan_path,
        outdir=run_root,
        workdir=tmp_path / "work",
        network_mode="direct",
        container_engine="apptainer",
        executor="local",
    )

    execution = run_root / "execution"
    execution.mkdir(parents=True)
    for name in ("report.html", "timeline.html", "dag.html"):
        (execution / name).write_text(f"<html>{name}</html>\n", encoding="utf-8")
    (execution / "trace.tsv").write_text(
        "name\tstatus\tduration\trealtime\t%cpu\tpeak_rss\n"
        "FASTQC (control_1)\tCOMPLETED\t2s\t1s\t90%\t100 MB\n"
        "STAR_ALIGN (treated_1)\tCACHED\t2m\t1m\t300%\t2 GB\n",
        encoding="utf-8",
    )
    rnaseq = run_root / "rnaseq"
    pipeline_info = rnaseq / "pipeline_info"
    pipeline_info.mkdir(parents=True)
    (pipeline_info / "software_versions.yml").write_text("STAR: 2.7.11b\n", encoding="utf-8")
    (pipeline_info / "params.json").write_text(
        json.dumps({"input": str(samplesheet.resolve()), "outdir": str(rnaseq.resolve())}),
        encoding="utf-8",
    )
    (pipeline_info / "samplesheet.valid.csv").write_text(
        samplesheet.read_text(encoding="utf-8"), encoding="utf-8"
    )
    quantification = rnaseq / "star_salmon"
    quantification.mkdir()
    header = "gene_id\t" + "\t".join(samples) + "\n"
    row = "ENSG000001\t1\t2\t3\t4\n"
    (quantification / "salmon.merged.gene_counts.tsv").write_text(header + row, encoding="utf-8")
    (quantification / "salmon.merged.gene_tpm.tsv").write_text(header + row, encoding="utf-8")
    multiqc = rnaseq / "multiqc" / "star_salmon"
    data = multiqc / "multiqc_data"
    data.mkdir(parents=True)
    (multiqc / "multiqc_report.html").write_text("<html>MultiQC</html>\n", encoding="utf-8")
    (data / "multiqc_data.json").write_text("{}\n", encoding="utf-8")
    return plan_path, run_root


def test_completed_run_creates_private_immutable_receipt(tmp_path: Path) -> None:
    plan, _ = completed_fixture(tmp_path)
    output = tmp_path / "private" / "primary-completion.json"
    progress: list[tuple[str, int, int]] = []

    result = create_primary_receipt(
        run_plan=plan,
        output=output,
        progress=lambda role, current, total: progress.append((role, current, total)),
    )

    assert result["stage"] == "pipeline_complete_qc_pending"
    assert result["qc_status"] == "pending_review"
    assert result["sample_count"] == 4
    assert result["task_summary"]["task_count"] == 2
    assert result["multiqc"]["data_storage"]["file_count"] == 1
    assert {item["role"] for item in result["artifacts"]} >= {
        "gene_counts", "gene_tpm", "multiqc_report", "software_versions"
    }
    assert progress
    assert stat.S_IMODE(output.stat().st_mode) == 0o600
    with pytest.raises(PrimaryDeliveryError, match="refusing to overwrite"):
        create_primary_receipt(run_plan=plan, output=output)


def test_failed_trace_is_not_completion(tmp_path: Path) -> None:
    plan, run_root = completed_fixture(tmp_path)
    trace = run_root / "execution" / "trace.tsv"
    trace.write_text(
        trace.read_text(encoding="utf-8").replace("COMPLETED", "FAILED"),
        encoding="utf-8",
    )

    with pytest.raises(PrimaryDeliveryError, match="not complete and successful"):
        create_primary_receipt(run_plan=plan, output=tmp_path / "receipt.json")


def test_missing_multiqc_report_is_rejected(tmp_path: Path) -> None:
    plan, run_root = completed_fixture(tmp_path)
    (run_root / "rnaseq" / "multiqc" / "star_salmon" / "multiqc_report.html").unlink()

    with pytest.raises(PrimaryDeliveryError, match="exactly one MultiQC report"):
        create_primary_receipt(run_plan=plan, output=tmp_path / "receipt.json")


def test_matrix_must_include_every_planned_sample(tmp_path: Path) -> None:
    plan, run_root = completed_fixture(tmp_path)
    matrix = run_root / "rnaseq" / "star_salmon" / "salmon.merged.gene_counts.tsv"
    matrix.write_text(
        matrix.read_text(encoding="utf-8").replace("\ttreated_2", ""),
        encoding="utf-8",
    )

    with pytest.raises(PrimaryDeliveryError, match="missing 1 planned sample"):
        create_primary_receipt(run_plan=plan, output=tmp_path / "receipt.json")


def test_pipeline_params_must_match_plan(tmp_path: Path) -> None:
    plan, run_root = completed_fixture(tmp_path)
    params = run_root / "rnaseq" / "pipeline_info" / "params.json"
    params.write_text(json.dumps({"input": "/wrong", "outdir": str(run_root / "rnaseq")}), encoding="utf-8")

    with pytest.raises(PrimaryDeliveryError, match="'input' does not match"):
        create_primary_receipt(run_plan=plan, output=tmp_path / "receipt.json")


def test_changed_control_file_is_rejected(tmp_path: Path) -> None:
    plan, _ = completed_fixture(tmp_path)
    payload = json.loads(plan.read_text(encoding="utf-8"))
    design = Path(payload["control_files"]["design"]["path"])
    design.write_text(design.read_text(encoding="utf-8") + "\n", encoding="utf-8")

    with pytest.raises(PrimaryDeliveryError, match="size changed after planning"):
        create_primary_receipt(run_plan=plan, output=tmp_path / "receipt.json")
