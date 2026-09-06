from __future__ import annotations

import json
import stat
from pathlib import Path

import pytest

from rnaseq_service.bundle import BundleError, seal_bundle, verify_bundle


ROOT = Path(__file__).parents[1]


def make_bundle(tmp_path: Path) -> tuple[Path, dict[str, Path]]:
    root = tmp_path / "offline"
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
    return root, {
        "rnaseq_workflow": rnaseq,
        "differential_workflow": differential,
        "container_root": containers,
        "plugin_root": plugins.parent,
    }


def seal(tmp_path: Path) -> tuple[Path, dict[str, object]]:
    root, components = make_bundle(tmp_path)
    payload = seal_bundle(
        bundle_dir=root,
        workflow_lock=ROOT / "config" / "workflows.toml",
        **components,
    )
    return root / "offline_bundle.manifest.json", payload


def test_seal_and_verify_exact_bundle(tmp_path: Path) -> None:
    manifest, payload = seal(tmp_path)

    result = verify_bundle(
        manifest,
        expected_workflow_lock=ROOT / "config" / "workflows.toml",
    )

    assert result["verified"] is True
    assert result["artifact_count"] == 6
    assert payload["workflows"]["rnaseq"]["revision"] == "3.26.0"
    assert stat.S_IMODE(manifest.stat().st_mode) == 0o600


def test_modified_artifact_is_rejected(tmp_path: Path) -> None:
    manifest, payload = seal(tmp_path)
    member = manifest.parent / payload["artifacts"][0]["path"]
    member.write_text("modified\n", encoding="utf-8")

    with pytest.raises(BundleError, match="size changed|checksum changed"):
        verify_bundle(manifest)


def test_missing_and_unexpected_artifacts_are_rejected(tmp_path: Path) -> None:
    manifest, payload = seal(tmp_path)
    member = manifest.parent / payload["artifacts"][0]["path"]
    member.unlink()

    with pytest.raises(BundleError, match="membership changed"):
        verify_bundle(manifest)

    member.write_text("restored", encoding="utf-8")
    (manifest.parent / "untracked.txt").write_text("surprise", encoding="utf-8")
    with pytest.raises(BundleError, match="membership changed"):
        verify_bundle(manifest)


def test_workflow_lock_mismatch_is_rejected(tmp_path: Path) -> None:
    manifest, _ = seal(tmp_path)
    other_lock = tmp_path / "other.toml"
    text = (ROOT / "config" / "workflows.toml").read_text(encoding="utf-8")
    other_lock.write_text(text.replace("3.26.0", "3.25.1"), encoding="utf-8")

    with pytest.raises(BundleError, match="do not match"):
        verify_bundle(manifest, expected_workflow_lock=other_lock)


def test_manifest_refuses_overwrite(tmp_path: Path) -> None:
    manifest, _ = seal(tmp_path)
    root = manifest.parent
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    components = payload["components"]

    with pytest.raises(BundleError, match="refusing to overwrite"):
        seal_bundle(
            bundle_dir=root,
            workflow_lock=ROOT / "config" / "workflows.toml",
            rnaseq_workflow=root / components["rnaseq_workflow"],
            differential_workflow=root / components["differential_workflow"],
            container_root=root / components["container_root"],
            plugin_root=root / components["plugin_root"],
        )


def test_symlink_is_rejected(tmp_path: Path) -> None:
    root, components = make_bundle(tmp_path)
    (root / "linked.sif").symlink_to(components["container_root"] / "tool.sif")

    with pytest.raises(BundleError, match="symbolic links"):
        seal_bundle(
            bundle_dir=root,
            workflow_lock=ROOT / "config" / "workflows.toml",
            **components,
        )
