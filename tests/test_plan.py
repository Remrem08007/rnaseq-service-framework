from __future__ import annotations

import csv
import json
import stat
from pathlib import Path

import pytest

from rnaseq_service.bundle import seal_bundle
from rnaseq_service.inputs import create_input_manifest
from rnaseq_service.plan import PlanError, create_rnaseq_plan
from rnaseq_service.workflow_lock import WorkflowLockError, load_workflow_lock


ROOT = Path(__file__).parents[1]


def write_csv(path: Path, headers: list[str], rows: list[list[str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        writer.writerows(rows)
    return path


def make_intake(tmp_path: Path) -> tuple[Path, Path, Path]:
    reads = tmp_path / "reads"
    reads.mkdir(parents=True)
    sample_rows: list[list[str]] = []
    design_rows: list[list[str]] = []
    for group in ("control", "treated"):
        for replicate in (1, 2):
            sample = f"{group}_{replicate}"
            r1 = reads / f"{sample}_R1.fastq.gz"
            r2 = reads / f"{sample}_R2.fastq.gz"
            r1.write_bytes(b"synthetic-r1\n")
            r2.write_bytes(b"synthetic-r2\n")
            sample_rows.append([sample, str(r1), str(r2), "auto"])
            design_rows.append([sample, group])
    samplesheet = write_csv(
        tmp_path / "samplesheet.csv",
        ["sample", "fastq_1", "fastq_2", "strandedness"],
        sample_rows,
    )
    design = write_csv(
        tmp_path / "design.csv",
        ["sample", "condition"],
        design_rows,
    )
    contrasts = write_csv(
        tmp_path / "contrasts.csv",
        ["contrast_id", "variable", "reference", "target"],
        [["treated_vs_control", "condition", "control", "treated"]],
    )
    return samplesheet, design, contrasts


def plan_kwargs(tmp_path: Path) -> dict[str, object]:
    samplesheet, design, contrasts = make_intake(tmp_path)
    return {
        "samplesheet": samplesheet,
        "design": design,
        "contrasts": contrasts,
        "workflow_lock": ROOT / "config" / "workflows.toml",
        "output": tmp_path / "private" / "run-plan.json",
        "outdir": tmp_path / "results",
        "workdir": tmp_path / "work",
        "network_mode": "direct",
        "container_engine": "apptainer",
        "executor": "local",
    }


def make_offline_manifest(tmp_path: Path) -> Path:
    root = tmp_path / "offline-bundle"
    rnaseq = root / "pipelines" / "rnaseq" / "workflow"
    differential = root / "pipelines" / "differential" / "workflow"
    containers = root / "containers"
    plugins = root / "plugins" / "nf-schema-2.5.1"
    for workflow in (rnaseq, differential):
        workflow.mkdir(parents=True)
        (workflow / "main.nf").write_text("nextflow.enable.dsl=2\n", encoding="utf-8")
        (workflow / "nextflow.config").write_text("plugins {}\n", encoding="utf-8")
    containers.mkdir(parents=True)
    (containers / "tool.sif").write_bytes(b"synthetic-container")
    plugins.mkdir(parents=True)
    (plugins / "plugin.jar").write_bytes(b"synthetic-plugin")
    seal_bundle(
        bundle_dir=root,
        workflow_lock=ROOT / "config" / "workflows.toml",
        rnaseq_workflow=rnaseq,
        differential_workflow=differential,
        container_root=containers,
        plugin_root=plugins.parent,
    )
    return root / "offline_bundle.manifest.json"


def test_repository_workflow_lock_is_exact() -> None:
    lock = load_workflow_lock(ROOT / "config" / "workflows.toml")

    assert lock.rnaseq.name == "nf-core/rnaseq"
    assert lock.rnaseq.revision == "3.26.0"
    assert lock.differential.name == "nf-core/differentialabundance"
    assert lock.differential.revision == "2.0.0"
    assert lock.runtime.nextflow_version == "26.04.6"
    assert lock.runtime.nf_core_tools_version == "4.1.0"


def test_floating_revision_is_rejected(tmp_path: Path) -> None:
    lock = tmp_path / "workflows.toml"
    lock.write_text(
        """[workflows]
rnaseq_name = "nf-core/rnaseq"
rnaseq_revision = "latest"
differential_name = "nf-core/differentialabundance"
differential_revision = "2.0.0"
[runtime]
nextflow_version = "26.04.6"
nf_core_tools_version = "4.1.0"
[policy]
versions_are_pinned = true
allow_development_revisions = false
last_reviewed = "2026-09-06"
""",
        encoding="utf-8",
    )

    with pytest.raises(WorkflowLockError, match="not an exact semantic version"):
        load_workflow_lock(lock)


def test_plan_records_exact_argv_hashes_and_no_results(tmp_path: Path) -> None:
    kwargs = plan_kwargs(tmp_path)

    plan = create_rnaseq_plan(**kwargs)

    assert plan["stage"] == "planned_not_executed"
    assert plan["workflow"] == {"name": "nf-core/rnaseq", "revision": "3.26.0"}
    assert plan["command_argv"][:6] == [
        "nextflow",
        "run",
        "nf-core/rnaseq",
        "-r",
        "3.26.0",
        "-profile",
    ]
    assert plan["contains_client_results"] is False
    assert plan["contains_secrets"] is False
    controls = plan["control_files"]
    assert set(controls) == {"samplesheet", "design", "contrasts", "workflow_lock"}
    assert all(len(item["sha256"]) == 64 for item in controls.values())


def test_proxy_values_are_never_serialized(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    kwargs = plan_kwargs(tmp_path)
    kwargs["network_mode"] = "proxy"
    monkeypatch.setenv("HTTPS_PROXY", "http://person:secret@example.invalid:8080")

    plan = create_rnaseq_plan(**kwargs)
    rendered = json.dumps(plan)

    assert "HTTPS_PROXY" in plan["execution"]["proxy_environment_present"]
    assert plan["execution"]["proxy_values_recorded"] is False
    assert "person:secret" not in rendered
    assert "example.invalid" not in rendered


def test_offline_mode_requires_and_verifies_local_bundle(tmp_path: Path) -> None:
    kwargs = plan_kwargs(tmp_path)
    kwargs["network_mode"] = "offline"

    with pytest.raises(PlanError, match="offline-manifest is required"):
        create_rnaseq_plan(**kwargs)

    manifest = make_offline_manifest(tmp_path)
    kwargs["offline_manifest"] = manifest
    plan = create_rnaseq_plan(**kwargs)

    execution = plan["execution"]
    assert execution["offline_bundle_verified"] is True
    assert execution["required_environment"]["NXF_OFFLINE"] == "true"
    assert execution["required_environment"]["NXF_PLUGINS_DIR"].endswith("/plugins")
    assert "-r" not in plan["command_argv"]
    assert plan["command_argv"][2].endswith("/pipelines/rnaseq/workflow")
    assert "offline_manifest" in plan["control_files"]


def test_offline_mode_rejects_modified_bundle(tmp_path: Path) -> None:
    kwargs = plan_kwargs(tmp_path)
    manifest = make_offline_manifest(tmp_path)
    (manifest.parent / "containers" / "tool.sif").write_bytes(b"changed")
    kwargs.update(network_mode="offline", offline_manifest=manifest)

    with pytest.raises(PlanError, match="offline bundle verification failed"):
        create_rnaseq_plan(**kwargs)


def test_plan_file_is_private_and_never_overwritten(tmp_path: Path) -> None:
    kwargs = plan_kwargs(tmp_path)
    output = kwargs["output"]

    create_rnaseq_plan(**kwargs)

    assert stat.S_IMODE(output.stat().st_mode) == 0o600
    with pytest.raises(PlanError, match="refusing to overwrite"):
        create_rnaseq_plan(**kwargs)


def test_slurm_requires_explicit_infrastructure_config(tmp_path: Path) -> None:
    kwargs = plan_kwargs(tmp_path)
    kwargs["executor"] = "slurm"

    with pytest.raises(PlanError, match="infrastructure-config is required"):
        create_rnaseq_plan(**kwargs)


def test_shell_metacharacters_remain_one_argv_value(tmp_path: Path) -> None:
    unusual = tmp_path / "client;do-not-run"
    unusual.mkdir()
    kwargs = plan_kwargs(unusual)

    plan = create_rnaseq_plan(**kwargs)

    samplesheet = str(kwargs["samplesheet"].resolve())
    index = plan["command_argv"].index("--input")
    assert plan["command_argv"][index + 1] == samplesheet
    assert "'" in plan["command_preview"]
    assert not (tmp_path / "do-not-run").exists()


def test_plan_binds_verified_input_manifest(tmp_path: Path) -> None:
    kwargs = plan_kwargs(tmp_path)
    manifest = tmp_path / "private" / "inputs.json"
    create_input_manifest(samplesheet=kwargs["samplesheet"], output=manifest)
    kwargs["input_manifest"] = manifest

    plan = create_rnaseq_plan(**kwargs)

    assert "input_manifest" in plan["control_files"]
    assert plan["input_dataset"]["n_fastq_files"] == 8
    assert plan["input_dataset"]["metadata_verified_at_planning"] is True
