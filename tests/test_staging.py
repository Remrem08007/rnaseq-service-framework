from __future__ import annotations

import json
import stat
from pathlib import Path

import pytest

from rnaseq_service.staging import (
    StagingError,
    _discover_workflow_dir,
    _redact,
    check_staging_environment,
    create_staging_plan,
)


ROOT = Path(__file__).parents[1]


def kwargs(tmp_path: Path) -> dict[str, object]:
    return {
        "workflow_lock": ROOT / "config" / "workflows.toml",
        "bundle_dir": tmp_path / "offline-bundle",
        "output": tmp_path / "private" / "staging-plan.json",
        "network_mode": "direct",
        "parallel_downloads": 6,
    }


def test_staging_plan_pins_both_workflows_and_shared_cache(tmp_path: Path) -> None:
    values = kwargs(tmp_path)

    plan = create_staging_plan(**values)

    assert plan["stage"] == "planned_not_executed"
    assert [item["workflow"]["revision"] for item in plan["commands"]] == [
        "3.26.0",
        "2.0.0",
    ]
    for item in plan["commands"]:
        argv = item["argv"]
        assert argv[:3] == ["nf-core", "pipelines", "download"]
        assert "--container-cache-utilisation" in argv
        assert argv[argv.index("--container-cache-utilisation") + 1] == "amend"
        assert argv[argv.index("--parallel-downloads") + 1] == "6"
        assert argv[argv.index("--download-configuration") + 1] == "yes"
    assert plan["runtime_environment"]["NXF_SINGULARITY_CACHEDIR"].endswith(
        "/offline-bundle/containers"
    )
    assert plan["runtime_environment"]["NXF_PLUGINS_DIR"].endswith(
        "/offline-bundle/plugins"
    )
    assert stat.S_IMODE(Path(values["output"]).stat().st_mode) == 0o600


def test_proxy_plan_records_names_but_not_values(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    values = kwargs(tmp_path)
    values["network_mode"] = "proxy"
    monkeypatch.setenv("https_proxy", "http://user:secret@proxy.invalid:8080")

    plan = create_staging_plan(**values)
    rendered = json.dumps(plan)

    assert "https_proxy" in plan["runtime_environment"]["proxy_environment_present"]
    assert plan["runtime_environment"]["proxy_values_recorded"] is False
    assert "user:secret" not in rendered
    assert "proxy.invalid" not in rendered


def test_proxy_plan_requires_runtime_configuration(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    for name in ("HTTP_PROXY", "HTTPS_PROXY", "NO_PROXY", "http_proxy", "https_proxy", "no_proxy"):
        monkeypatch.delenv(name, raising=False)
    values = kwargs(tmp_path)
    values["network_mode"] = "proxy"

    with pytest.raises(StagingError, match="requires a standard proxy"):
        create_staging_plan(**values)


def test_staging_plan_refuses_overwrite(tmp_path: Path) -> None:
    values = kwargs(tmp_path)
    create_staging_plan(**values)

    with pytest.raises(StagingError, match="refusing to overwrite"):
        create_staging_plan(**values)


def test_redaction_removes_every_proxy_value() -> None:
    output = "connecting through http://user:pw@proxy.invalid:8080\n"

    rendered = _redact(output, ["http://user:pw@proxy.invalid:8080"])

    assert rendered == "connecting through <redacted-proxy>\n"


def test_workflow_discovery_ignores_nested_workflow_fixtures(tmp_path: Path) -> None:
    output_root = tmp_path / "rnaseq-3.26.0"
    workflow = output_root / "3_26_0"
    nested = workflow / "tests" / "fixture"
    for directory in (workflow, nested):
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "main.nf").write_text("workflow {}\n", encoding="utf-8")
        (directory / "nextflow.config").write_text("params {}\n", encoding="utf-8")

    assert _discover_workflow_dir(output_root) == workflow.resolve()


def test_workflow_discovery_rejects_ambiguous_top_level_roots(tmp_path: Path) -> None:
    output_root = tmp_path / "download"
    for name in ("first", "second"):
        directory = output_root / name
        directory.mkdir(parents=True)
        (directory / "main.nf").write_text("workflow {}\n", encoding="utf-8")
        (directory / "nextflow.config").write_text("params {}\n", encoding="utf-8")

    with pytest.raises(StagingError, match="top-level downloaded workflow.*found 2"):
        _discover_workflow_dir(output_root)


def test_workflow_discovery_rejects_missing_root(tmp_path: Path) -> None:
    output_root = tmp_path / "empty"
    output_root.mkdir()

    with pytest.raises(StagingError, match="no downloaded workflow"):
        _discover_workflow_dir(output_root)


def test_environment_check_verifies_versions_and_endpoints(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    values = kwargs(tmp_path)
    create_staging_plan(**values)

    class Response:
        status = 200

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

    monkeypatch.setattr("rnaseq_service.staging.shutil.which", lambda *args, **kwargs: "/bin/tool")
    monkeypatch.setattr(
        "rnaseq_service.staging.subprocess.check_output",
        lambda argv, **kwargs: "nf-core 4.1.0" if argv[0] == "nf-core" else "Nextflow 25.10.2",
    )
    monkeypatch.setattr(
        "rnaseq_service.staging.urllib.request.urlopen",
        lambda request, timeout: Response(),
    )

    result = check_staging_environment(Path(values["output"]))

    assert result["ready"] is True
    assert all(check["ok"] for check in result["checks"])
    assert result["proxy_values_recorded"] is False
