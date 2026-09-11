from __future__ import annotations

import json
import stat
import subprocess
from pathlib import Path

import pytest

import rnaseq_service.upstream_data as upstream_data
from rnaseq_service.bundle import seal_bundle
from rnaseq_service.upstream_smoke import (
    UpstreamSmokeError,
    create_upstream_smoke_launcher,
    create_upstream_smoke_plan,
    create_upstream_smoke_receipt,
    inspect_upstream_smoke_completion,
)
from rnaseq_service.upstream_smoke_cli import main as smoke_main
from rnaseq_service.upstream_data import stage_upstream_test_data
from test_upstream_data import FakeOpener


ROOT = Path(__file__).parents[1]


def _file(path: Path, content: str = "evidence\n") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _executable(path: Path, content: str) -> None:
    _file(path, content)
    path.chmod(0o700)


def _complete_fixture(plan: dict[str, object]) -> None:
    run_root = Path(str(plan["execution"]["run_root"]))
    _file(
        run_root / "runtime" / "versions.tsv",
        "tool\tversion\n"
        f"nextflow\t{plan['runtime']['nextflow_version']}\n"
        f"container_engine\t{plan['execution']['container_engine']}\n"
        "container_runtime\tapptainer version 1.4.5\n",
    )
    requested = str(plan["execution"]["network_mode"])
    offline = bool(plan["execution"]["offline_bundle_verified"])
    selected = "offline" if offline else ("direct" if requested == "auto" else requested)
    _file(
        run_root / "runtime" / "network_selection.tsv",
        "requested_mode\tselected_mode\toffline_bundle_available\n"
        f"{requested}\t{selected}\t{str(offline).lower()}\n",
    )
    trace = (
        "task_id\tname\tstatus\trealtime\t%cpu\tpeak_rss\tduration\n"
        "1\tTASK\tCOMPLETED\t1s\t100%\t1 MB\t1s\n"
    )
    for stage_id in ("rnaseq", "differential"):
        stage_root = run_root / stage_id
        execution = stage_root / "execution"
        results = stage_root / "results"
        _file(execution / "trace.tsv", trace)
        for name in ("report.html", "timeline.html", "dag.html"):
            _file(execution / name)
        _file(
            results / "pipeline_info" / "params.json",
            json.dumps({"outdir": str(results.resolve())}),
        )
        version_name = (
            "nf_core_rnaseq_software_mqc_versions.yml"
            if stage_id == "rnaseq"
            else "nf_core_differentialabundance_software_versions.yml"
        )
        _file(results / "pipeline_info" / version_name)
    rnaseq = run_root / "rnaseq" / "results"
    _file(rnaseq / "multiqc" / "star_salmon" / "multiqc_report.html")
    _file(rnaseq / "star_salmon" / "salmon.merged.gene_counts.tsv")
    _file(rnaseq / "star_salmon" / "salmon.merged.gene_tpm.tsv")
    differential = run_root / "differential" / "results"
    prefix = "treatment_mCherry_hND6__SRP254919"
    _file(
        differential
        / "report"
        / "rnaseq_deseq2_gsea"
        / "SRP254919_differentialabundance_report.html"
    )
    _file(
        differential
        / "tables"
        / "processed_abundance"
        / "rnaseq_deseq2_gsea"
        / "all.normalised_counts.tsv"
    )
    _file(
        differential
        / "tables"
        / "processed_abundance"
        / "rnaseq_deseq2_gsea"
        / "all.vst.tsv"
    )
    _file(
        differential
        / "tables"
        / "differential"
        / "rnaseq_deseq2_gsea"
        / f"{prefix}.deseq2.results.tsv"
    )
    _file(
        differential
        / "tables"
        / "differential"
        / "rnaseq_deseq2_gsea"
        / f"{prefix}.deseq2.results_filtered.tsv"
    )
    _file(
        differential
        / "plots"
        / "differential"
        / "rnaseq_deseq2_gsea"
        / prefix
        / "png"
        / "volcano.png"
    )
    _file(
        differential
        / "report"
        / "gsea"
        / "rnaseq_deseq2_gsea"
        / prefix
        / "hallmarks"
        / f"{prefix}.hallmarks.Gsea.rpt"
    )


def _offline_bundle(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    root = tmp_path / "offline"
    for workflow in ("rnaseq", "differential"):
        _file(root / "pipelines" / workflow / "main.nf")
        _file(root / "pipelines" / workflow / "nextflow.config")
    _file(root / "containers" / "tool.sif")
    _file(root / "plugins" / "nf-validation.jar")
    monkeypatch.setattr(
        upstream_data, "select_network_mode", lambda *args, **kwargs: "direct"
    )
    monkeypatch.setattr(upstream_data, "_opener", lambda mode: FakeOpener())
    data_root = root / "upstream-test-data"
    stage_upstream_test_data(output_root=data_root)
    seal_bundle(
        bundle_dir=root,
        workflow_lock=ROOT / "config" / "workflows.toml",
        rnaseq_workflow=root / "pipelines" / "rnaseq",
        differential_workflow=root / "pipelines" / "differential",
        container_root=root / "containers",
        plugin_root=root / "plugins",
        upstream_test_data_root=data_root,
    )
    return root / "offline_bundle.manifest.json"


def test_plan_pins_both_official_public_test_profiles(tmp_path: Path) -> None:
    output = tmp_path / "private" / "upstream-smoke.json"
    plan = create_upstream_smoke_plan(
        workflow_lock=ROOT / "config" / "workflows.toml",
        output=output,
        run_root=tmp_path / "runs",
        network_mode="direct",
    )

    assert plan["stage"] == "upstream_smoke_planned_not_executed"
    assert plan["uses_public_upstream_test_data"] is True
    assert plan["scientific_execution_expected"] is True
    assert plan["contains_client_data"] is False
    assert plan["contains_secrets"] is False
    assert [stage["id"] for stage in plan["stages"]] == ["rnaseq", "differential"]
    rnaseq, differential = plan["stages"]
    assert rnaseq["workflow"] == {"name": "nf-core/rnaseq", "revision": "3.26.0"}
    assert differential["workflow"] == {
        "name": "nf-core/differentialabundance",
        "revision": "2.0.0",
    }
    assert "test,apptainer" in rnaseq["command_argv"]
    assert "test_rnaseq_deseq2_gsea,apptainer" in differential["command_argv"]
    assert all("-resume" in stage["command_argv"] for stage in plan["stages"])
    assert stat.S_IMODE(output.stat().st_mode) == 0o600


def test_plan_requires_offline_manifest_and_rejects_overwrite(tmp_path: Path) -> None:
    kwargs = {
        "workflow_lock": ROOT / "config" / "workflows.toml",
        "output": tmp_path / "plan.json",
        "run_root": tmp_path / "runs",
        "network_mode": "offline",
    }
    with pytest.raises(UpstreamSmokeError, match="offline-manifest is required"):
        create_upstream_smoke_plan(**kwargs)
    kwargs["network_mode"] = "direct"
    create_upstream_smoke_plan(**kwargs)
    with pytest.raises(UpstreamSmokeError, match="refusing to overwrite"):
        create_upstream_smoke_plan(**kwargs)


def test_auto_plan_prefers_verified_offline_bundle(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest = _offline_bundle(tmp_path, monkeypatch)
    plan = create_upstream_smoke_plan(
        workflow_lock=ROOT / "config" / "workflows.toml",
        output=tmp_path / "plan.json",
        run_root=tmp_path / "runs",
        network_mode="auto",
        offline_manifest=manifest,
    )
    execution = plan["execution"]
    assert execution["selected_network_mode"] == "offline"
    assert execution["offline_bundle_verified"] is True
    assert execution["required_environment"]["NXF_OFFLINE"] == "true"
    for stage in plan["stages"]:
        assert "-r" not in stage["command_argv"]
        assert stage["command_argv"][4].startswith(str(manifest.parent.resolve()))
        assert all("https://" not in value for value in stage["command_argv"])

    launcher = tmp_path / "smoke.sbatch"
    create_upstream_smoke_launcher(
        run_plan=tmp_path / "plan.json",
        output=launcher,
        account="def-project",
        partition="compute",
        wall_time="04:00:00",
        cpus=8,
        memory_gb=32,
        modules=["nextflow", "apptainer"],
    )
    rendered = launcher.read_text(encoding="utf-8")
    assert "export NXF_OFFLINE=true" in rendered
    assert "curl --silent" not in rendered
    assert "selected_network=offline" in rendered


def test_offline_smoke_rejects_software_only_bundle(tmp_path: Path) -> None:
    root = tmp_path / "software-only"
    for workflow in ("rnaseq", "differential"):
        _file(root / "pipelines" / workflow / "main.nf")
        _file(root / "pipelines" / workflow / "nextflow.config")
    _file(root / "containers" / "tool.sif")
    _file(root / "plugins" / "plugin.jar")
    seal_bundle(
        bundle_dir=root,
        workflow_lock=ROOT / "config" / "workflows.toml",
        rnaseq_workflow=root / "pipelines/rnaseq",
        differential_workflow=root / "pipelines/differential",
        container_root=root / "containers",
        plugin_root=root / "plugins",
    )
    with pytest.raises(UpstreamSmokeError, match="no verified upstream test-data"):
        create_upstream_smoke_plan(
            workflow_lock=ROOT / "config" / "workflows.toml",
            output=tmp_path / "plan.json",
            run_root=tmp_path / "runs",
            network_mode="offline",
            offline_manifest=root / "offline_bundle.manifest.json",
        )


def test_prepare_reverifies_offline_test_data(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest = _offline_bundle(tmp_path, monkeypatch)
    plan = tmp_path / "plan.json"
    create_upstream_smoke_plan(
        workflow_lock=ROOT / "config" / "workflows.toml",
        output=plan,
        run_root=tmp_path / "runs",
        network_mode="offline",
        offline_manifest=manifest,
    )
    (manifest.parent / "upstream-test-data/differential/mouse_network.tsv").write_text(
        "tampered\n", encoding="utf-8"
    )
    with pytest.raises(UpstreamSmokeError, match="offline bundle verification failed"):
        create_upstream_smoke_launcher(
            run_plan=plan,
            output=tmp_path / "smoke.sbatch",
            account="def-project",
            partition="compute",
            wall_time="04:00:00",
            cpus=8,
            memory_gb=32,
            modules=["nextflow", "apptainer"],
        )
def test_upstream_smoke_cli_plans_without_running(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    output = tmp_path / "plan.json"
    assert smoke_main(
        [
            "plan",
            "--workflow-lock",
            str(ROOT / "config" / "workflows.toml"),
            "--output",
            str(output),
            "--run-root",
            str(tmp_path / "runs"),
            "--network-mode",
            "direct",
        ]
    ) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary == {
        "offline_bundle_verified": False,
        "requested_network_mode": "direct",
        "scientific_execution_expected": True,
        "selected_network_mode": "direct",
        "stages": ["rnaseq", "differential"],
        "status": "upstream_smoke_planned_not_executed",
        "uses_public_upstream_test_data": True,
    }


def test_prepare_renders_non_submitting_heartbeat_launcher(tmp_path: Path) -> None:
    plan_path = tmp_path / "private" / "plan.json"
    create_upstream_smoke_plan(
        workflow_lock=ROOT / "config" / "workflows.toml",
        output=plan_path,
        run_root=tmp_path / "runs",
        network_mode="direct",
    )
    launcher = tmp_path / "private" / "launch" / "smoke.sbatch"

    result = create_upstream_smoke_launcher(
        run_plan=plan_path,
        output=launcher,
        account="def-project",
        partition="compute",
        wall_time="04:00:00",
        cpus=8,
        memory_gb=32,
        modules=["StdEnv/2023", "nextflow/25.10.2", "apptainer/1.3.5"],
    )

    text = launcher.read_text(encoding="utf-8")
    assert result["prepared"] is True
    assert result["submitted"] is False
    assert result["heartbeat_seconds"] == 60
    assert "#SBATCH --account=def-project" in text
    assert "stage=${ordinal}/2" in text
    assert "Nextflow version mismatch" in text
    assert "runtime/versions.tsv" in text
    assert "test,apptainer" in text
    assert "test_rnaseq_deseq2_gsea,apptainer" in text
    assert text.count("run_stage ") == 2
    assert "sbatch" not in text
    assert stat.S_IMODE(launcher.stat().st_mode) == 0o700
    assert subprocess.run(["bash", "-n", str(launcher)], check=False).returncode == 0


def test_offline_launcher_executes_without_calling_curl(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest = _offline_bundle(tmp_path, monkeypatch)
    plan_path = tmp_path / "plan.json"
    plan = create_upstream_smoke_plan(
        workflow_lock=ROOT / "config" / "workflows.toml",
        output=plan_path,
        run_root=tmp_path / "runs",
        network_mode="auto",
        offline_manifest=manifest,
    )
    launcher = tmp_path / "smoke.sbatch"
    create_upstream_smoke_launcher(
        run_plan=plan_path,
        output=launcher,
        account="def-project",
        partition="compute",
        wall_time="04:00:00",
        cpus=2,
        memory_gb=4,
        modules=["nextflow", "apptainer"],
        purge_modules=False,
    )
    binaries = tmp_path / "bin"
    marker = tmp_path / "curl-called"
    _executable(binaries / "module", "#!/usr/bin/env bash\nexit 0\n")
    _executable(
        binaries / "nextflow",
        "#!/usr/bin/env bash\n"
        "if [[ ${1:-} == -version ]]; then echo 'nextflow version 25.10.2'; fi\n"
        "exit 0\n",
    )
    _executable(binaries / "apptainer", "#!/usr/bin/env bash\necho 'apptainer version 1.4.5'\n")
    _executable(binaries / "curl", f"#!/usr/bin/env bash\ntouch {marker}\nexit 0\n")
    result = subprocess.run(
        ["bash", str(launcher)],
        env={"PATH": f"{binaries}:/usr/bin:/bin"},
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert not marker.exists()
    evidence = (
        Path(str(plan["execution"]["run_root"]))
        / "runtime"
        / "network_selection.tsv"
    ).read_text(encoding="utf-8")
    assert "auto\toffline\ttrue" in evidence


def test_auto_launcher_falls_back_to_working_proxy(
    tmp_path: Path,
) -> None:
    plan_path = tmp_path / "plan.json"
    plan = create_upstream_smoke_plan(
        workflow_lock=ROOT / "config" / "workflows.toml",
        output=plan_path,
        run_root=tmp_path / "runs",
        network_mode="auto",
    )
    launcher = tmp_path / "smoke.sbatch"
    create_upstream_smoke_launcher(
        run_plan=plan_path,
        output=launcher,
        account="def-project",
        partition="compute",
        wall_time="04:00:00",
        cpus=2,
        memory_gb=4,
        modules=["nextflow", "apptainer"],
        purge_modules=False,
    )
    binaries = tmp_path / "bin"
    _executable(binaries / "module", "#!/usr/bin/env bash\nexit 0\n")
    _executable(
        binaries / "nextflow",
        "#!/usr/bin/env bash\n"
        "if [[ ${1:-} == -version ]]; then echo 'nextflow version 25.10.2'; fi\n"
        "exit 0\n",
    )
    _executable(binaries / "apptainer", "#!/usr/bin/env bash\necho 'apptainer version 1.4.5'\n")
    _executable(
        binaries / "curl",
        "#!/usr/bin/env bash\n"
        "if [[ -n ${HTTPS_PROXY:-}${https_proxy:-} ]]; then exit 0; fi\n"
        "exit 7\n",
    )
    result = subprocess.run(
        ["bash", str(launcher)],
        env={
            "PATH": f"{binaries}:/usr/bin:/bin",
            "HTTPS_PROXY": "http://proxy.invalid:8080",
        },
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    evidence = (
        Path(str(plan["execution"]["run_root"]))
        / "runtime"
        / "network_selection.tsv"
    ).read_text(encoding="utf-8")
    assert "auto\tproxy\tfalse" in evidence


def test_prepare_rejects_changed_plan_control_and_overwrite(tmp_path: Path) -> None:
    plan_path = tmp_path / "plan.json"
    lock = tmp_path / "workflows.toml"
    lock.write_text(
        (ROOT / "config" / "workflows.toml").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    create_upstream_smoke_plan(
        workflow_lock=lock,
        output=plan_path,
        run_root=tmp_path / "runs",
        network_mode="direct",
    )
    kwargs = {
        "run_plan": plan_path,
        "output": tmp_path / "smoke.sbatch",
        "account": "def-project",
        "partition": "compute",
        "wall_time": "04:00:00",
        "cpus": 8,
        "memory_gb": 32,
        "modules": ["nextflow", "apptainer"],
    }
    create_upstream_smoke_launcher(**kwargs)
    with pytest.raises(UpstreamSmokeError, match="refusing to overwrite"):
        create_upstream_smoke_launcher(**kwargs)

    original = lock.read_text(encoding="utf-8")
    try:
        lock.write_text(original + "\n", encoding="utf-8")
        with pytest.raises(UpstreamSmokeError, match="workflow lock .* changed"):
            create_upstream_smoke_launcher(**{**kwargs, "output": tmp_path / "other.sbatch"})
    finally:
        lock.write_text(original, encoding="utf-8")


def test_prepare_cli(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    plan_path = tmp_path / "plan.json"
    create_upstream_smoke_plan(
        workflow_lock=ROOT / "config" / "workflows.toml",
        output=plan_path,
        run_root=tmp_path / "runs",
        network_mode="proxy",
    )
    launcher = tmp_path / "smoke.sbatch"
    assert smoke_main(
        [
            "prepare",
            "--run-plan",
            str(plan_path),
            "--output",
            str(launcher),
            "--account",
            "def-project",
            "--partition",
            "compute",
            "--time",
            "04:00:00",
            "--module",
            "nextflow",
            "--module",
            "apptainer",
        ]
    ) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["prepared"] is True
    assert summary["submitted"] is False


def test_completion_requires_both_successful_scientific_stages(tmp_path: Path) -> None:
    plan_path = tmp_path / "plan.json"
    plan = create_upstream_smoke_plan(
        workflow_lock=ROOT / "config" / "workflows.toml",
        output=plan_path,
        run_root=tmp_path / "runs",
        network_mode="direct",
    )
    _complete_fixture(plan)
    progress: list[tuple[str, int, int]] = []

    result = inspect_upstream_smoke_completion(
        plan_path,
        progress=lambda role, current, total: progress.append((role, current, total)),
    )

    assert result["stage"] == "upstream_smoke_complete"
    assert result["scientific_execution_performed"] is True
    assert result["uses_public_upstream_test_data"] is True
    assert result["contains_client_data"] is False
    assert [stage["id"] for stage in result["stages"]] == ["rnaseq", "differential"]
    assert all(stage["task_summary"]["task_count"] == 1 for stage in result["stages"])
    assert {artifact["role"] for artifact in result["stages"][0]["artifacts"]} >= {
        "multiqc_report",
        "gene_counts",
        "gene_tpm",
    }
    assert {artifact["role"] for artifact in result["stages"][1]["artifacts"]} >= {
        "deseq2_results",
        "volcano_plot",
        "gsea_report",
    }
    assert progress


def test_completion_receipt_is_immutable_and_private(tmp_path: Path) -> None:
    plan_path = tmp_path / "plan.json"
    plan = create_upstream_smoke_plan(
        workflow_lock=ROOT / "config" / "workflows.toml",
        output=plan_path,
        run_root=tmp_path / "runs",
        network_mode="proxy",
    )
    _complete_fixture(plan)
    receipt = tmp_path / "private" / "upstream-smoke.complete.json"

    create_upstream_smoke_receipt(run_plan=plan_path, output=receipt)

    assert stat.S_IMODE(receipt.stat().st_mode) == 0o600
    with pytest.raises(UpstreamSmokeError, match="refusing to overwrite"):
        create_upstream_smoke_receipt(run_plan=plan_path, output=receipt)


def test_completion_rejects_nonterminal_trace_and_missing_science(tmp_path: Path) -> None:
    plan_path = tmp_path / "plan.json"
    plan = create_upstream_smoke_plan(
        workflow_lock=ROOT / "config" / "workflows.toml",
        output=plan_path,
        run_root=tmp_path / "runs",
        network_mode="direct",
    )
    _complete_fixture(plan)
    trace = Path(str(plan["execution"]["run_root"])) / "rnaseq" / "execution" / "trace.tsv"
    trace.write_text(
        "task_id\tname\tstatus\trealtime\t%cpu\tpeak_rss\tduration\n"
        "1\tTASK\tRUNNING\t1s\t100%\t1 MB\t1s\n",
        encoding="utf-8",
    )
    with pytest.raises(UpstreamSmokeError, match="not complete and successful"):
        inspect_upstream_smoke_completion(plan_path)

    trace.write_text(
        "task_id\tname\tstatus\trealtime\t%cpu\tpeak_rss\tduration\n"
        "1\tTASK\tCOMPLETED\t1s\t100%\t1 MB\t1s\n",
        encoding="utf-8",
    )
    gsea = next(
        (Path(str(plan["execution"]["run_root"])) / "differential" / "results").glob(
            "report/gsea/**/*.Gsea.rpt"
        )
    )
    gsea.unlink()
    with pytest.raises(UpstreamSmokeError, match="GSEA report"):
        inspect_upstream_smoke_completion(plan_path)


def test_completion_rejects_runtime_version_mismatch(tmp_path: Path) -> None:
    plan_path = tmp_path / "plan.json"
    plan = create_upstream_smoke_plan(
        workflow_lock=ROOT / "config" / "workflows.toml",
        output=plan_path,
        run_root=tmp_path / "runs",
        network_mode="direct",
    )
    _complete_fixture(plan)
    runtime = Path(str(plan["execution"]["run_root"])) / "runtime" / "versions.tsv"
    runtime.write_text(
        "tool\tversion\n"
        "nextflow\t26.04.6\n"
        "container_engine\tapptainer\n"
        "container_runtime\tapptainer version 1.4.5\n",
        encoding="utf-8",
    )

    with pytest.raises(UpstreamSmokeError, match="Nextflow version"):
        inspect_upstream_smoke_completion(plan_path)


def test_complete_cli_reports_real_upstream_scope(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    plan_path = tmp_path / "plan.json"
    plan = create_upstream_smoke_plan(
        workflow_lock=ROOT / "config" / "workflows.toml",
        output=plan_path,
        run_root=tmp_path / "runs",
        network_mode="direct",
    )
    _complete_fixture(plan)
    output = tmp_path / "receipt.json"

    assert smoke_main(
        ["complete", "--run-plan", str(plan_path), "--output", str(output), "--quiet"]
    ) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["status"] == "upstream_smoke_complete"
    assert summary["scientific_execution_performed"] is True
    assert [stage["id"] for stage in summary["stages"]] == ["rnaseq", "differential"]
