"""Seal and verify checksum-bound offline workflow bundles."""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Callable

from .upstream_data import UpstreamDataError, verify_upstream_test_data
from .workflow_lock import WorkflowLockError, load_workflow_lock


MANIFEST_NAME = "offline_bundle.manifest.json"
ProgressCallback = Callable[[int, int, str], None]


class BundleError(ValueError):
    """Raised when an offline bundle is unsafe, incomplete, or modified."""


def _inside(root: Path, candidate: Path, label: str) -> Path:
    try:
        resolved = candidate.resolve(strict=True)
    except FileNotFoundError as exc:
        raise BundleError(f"{label} does not exist: {candidate}") from exc
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise BundleError(f"{label} must be inside the bundle: {candidate}") from exc
    return resolved


def _safe_relative(value: object) -> PurePosixPath:
    if not isinstance(value, str) or not value:
        raise BundleError("bundle file paths must be non-empty strings")
    path = PurePosixPath(value)
    if path.is_absolute() or ".." in path.parts or str(path) in {".", ""}:
        raise BundleError(f"unsafe bundle-relative path: {value!r}")
    return path


def _sha256(
    path: Path,
    *,
    offset: int,
    total: int,
    progress: ProgressCallback | None,
) -> str:
    digest = hashlib.sha256()
    consumed = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
            consumed += len(chunk)
            if progress is not None:
                progress(offset + consumed, total, path.name)
    if progress is not None and path.stat().st_size == 0:
        progress(offset, total, path.name)
    return digest.hexdigest()


def _bundle_files(root: Path, manifest: Path) -> list[Path]:
    files: list[Path] = []
    for path in sorted(root.rglob("*")):
        if path == manifest:
            continue
        if path.is_symlink():
            raise BundleError(f"symbolic links are not allowed in bundles: {path}")
        if path.is_file():
            files.append(path)
    if not files:
        raise BundleError("bundle contains no artifacts")
    return files


def _validate_components(
    root: Path,
    *,
    rnaseq_workflow: Path,
    differential_workflow: Path,
    container_root: Path,
    plugin_root: Path,
    upstream_test_data_root: Path | None = None,
) -> dict[str, str]:
    rnaseq = _inside(root, rnaseq_workflow, "RNA-seq workflow directory")
    differential = _inside(
        root, differential_workflow, "differential workflow directory"
    )
    containers = _inside(root, container_root, "container directory")
    plugins = _inside(root, plugin_root, "Nextflow plugin directory")
    for label, workflow in (("RNA-seq", rnaseq), ("differential", differential)):
        if not workflow.is_dir():
            raise BundleError(f"{label} workflow component is not a directory")
        for required in ("main.nf", "nextflow.config"):
            if not (workflow / required).is_file():
                raise BundleError(f"{label} workflow is missing {required}")
    if not containers.is_dir():
        raise BundleError("container component is not a directory")
    images = [
        path
        for path in containers.rglob("*")
        if path.is_file() and path.suffix.lower() in {".img", ".sif"}
    ]
    if not images:
        raise BundleError("container directory has no .img or .sif images")
    if not plugins.is_dir() or not any(path.is_file() for path in plugins.rglob("*")):
        raise BundleError("Nextflow plugin directory contains no plugin artifacts")
    components = {
        "rnaseq_workflow": rnaseq.relative_to(root).as_posix(),
        "differential_workflow": differential.relative_to(root).as_posix(),
        "container_root": containers.relative_to(root).as_posix(),
        "plugin_root": plugins.relative_to(root).as_posix(),
    }
    if upstream_test_data_root is not None:
        test_data = _inside(root, upstream_test_data_root, "upstream test-data directory")
        if not test_data.is_dir():
            raise BundleError("upstream test-data component is not a directory")
        try:
            verify_upstream_test_data(test_data)
        except (UpstreamDataError, OSError, UnicodeError) as exc:
            raise BundleError(f"upstream test-data verification failed: {exc}") from exc
        components["upstream_test_data_root"] = test_data.relative_to(root).as_posix()
    return components


def seal_bundle(
    *,
    bundle_dir: Path,
    workflow_lock: Path,
    rnaseq_workflow: Path,
    differential_workflow: Path,
    container_root: Path,
    plugin_root: Path,
    upstream_test_data_root: Path | None = None,
    manifest: Path | None = None,
    progress: ProgressCallback | None = None,
) -> dict[str, object]:
    """Hash all staged artifacts and exclusively create a bundle manifest."""

    root = bundle_dir.resolve(strict=True)
    if not root.is_dir():
        raise BundleError(f"bundle path is not a directory: {bundle_dir}")
    manifest_path = (manifest or (root / MANIFEST_NAME)).resolve()
    if manifest_path.parent != root:
        raise BundleError("bundle manifest must be written at the bundle root")
    if manifest_path.exists():
        raise BundleError(f"refusing to overwrite existing manifest: {manifest_path}")
    components = _validate_components(
        root,
        rnaseq_workflow=rnaseq_workflow,
        differential_workflow=differential_workflow,
        container_root=container_root,
        plugin_root=plugin_root,
        upstream_test_data_root=upstream_test_data_root,
    )
    try:
        lock = load_workflow_lock(workflow_lock)
    except WorkflowLockError as exc:
        raise BundleError(str(exc)) from exc

    files = _bundle_files(root, manifest_path)
    total_bytes = sum(path.stat().st_size for path in files)
    offset = 0
    records: list[dict[str, object]] = []
    for path in files:
        size = path.stat().st_size
        records.append(
            {
                "path": path.relative_to(root).as_posix(),
                "size_bytes": size,
                "sha256": _sha256(
                    path,
                    offset=offset,
                    total=total_bytes,
                    progress=progress,
                ),
            }
        )
        offset += size

    payload: dict[str, object] = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "bundle_status": "sealed",
        "workflows": {
            "rnaseq": asdict(lock.rnaseq),
            "differential": asdict(lock.differential),
        },
        "runtime": asdict(lock.runtime),
        "components": components,
        "artifacts": records,
        "artifact_count": len(records),
        "total_bytes": total_bytes,
        "contains_client_data": False,
        "contains_secrets": False,
    }
    try:
        descriptor = os.open(
            manifest_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600
        )
    except FileExistsError as exc:
        raise BundleError(f"refusing to overwrite existing manifest: {manifest_path}") from exc
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
    except Exception:
        manifest_path.unlink(missing_ok=True)
        raise
    return payload


def verify_bundle(
    manifest: Path,
    *,
    expected_workflow_lock: Path | None = None,
    progress: ProgressCallback | None = None,
) -> dict[str, object]:
    """Verify paths, exact file membership, sizes, hashes, and workflow pins."""

    manifest_path = manifest.resolve(strict=True)
    root = manifest_path.parent
    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise BundleError(f"could not read bundle manifest: {exc}") from exc
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise BundleError("unsupported bundle manifest schema")
    if payload.get("bundle_status") != "sealed":
        raise BundleError("bundle manifest is not sealed")
    if payload.get("contains_client_data") is not False:
        raise BundleError("bundle is not marked free of client data")
    if payload.get("contains_secrets") is not False:
        raise BundleError("bundle is not marked free of secrets")

    artifacts = payload.get("artifacts")
    if not isinstance(artifacts, list) or not artifacts:
        raise BundleError("bundle manifest has no artifacts")
    records: dict[str, dict[str, object]] = {}
    for record in artifacts:
        if not isinstance(record, dict):
            raise BundleError("invalid artifact record")
        relative = _safe_relative(record.get("path")).as_posix()
        if relative == manifest_path.name or relative in records:
            raise BundleError(f"duplicate or reserved artifact path: {relative}")
        records[relative] = record

    actual_files = {
        path.relative_to(root).as_posix()
        for path in _bundle_files(root, manifest_path)
    }
    expected_files = set(records)
    if actual_files != expected_files:
        missing = sorted(expected_files - actual_files)
        unexpected = sorted(actual_files - expected_files)
        raise BundleError(
            f"bundle membership changed; missing={missing}, unexpected={unexpected}"
        )

    total_bytes = sum(int(records[name].get("size_bytes", -1)) for name in records)
    if total_bytes < 0 or payload.get("total_bytes") != total_bytes:
        raise BundleError("bundle total byte count is invalid")
    offset = 0
    for relative in sorted(records):
        record = records[relative]
        path = root / relative
        size = path.stat().st_size
        if size != record.get("size_bytes"):
            raise BundleError(f"artifact size changed: {relative}")
        observed = _sha256(
            path,
            offset=offset,
            total=total_bytes,
            progress=progress,
        )
        if observed != record.get("sha256"):
            raise BundleError(f"artifact checksum changed: {relative}")
        offset += size

    components = payload.get("components")
    if not isinstance(components, dict):
        raise BundleError("bundle components are missing")
    try:
        test_data_component = components.get("upstream_test_data_root")
        checked_components = _validate_components(
            root,
            rnaseq_workflow=root / str(components["rnaseq_workflow"]),
            differential_workflow=root / str(components["differential_workflow"]),
            container_root=root / str(components["container_root"]),
            plugin_root=root / str(components["plugin_root"]),
            upstream_test_data_root=(
                root / test_data_component if test_data_component is not None else None
            ),
        )
    except KeyError as exc:
        raise BundleError(f"bundle component is missing: {exc.args[0]}") from exc

    if expected_workflow_lock is not None:
        try:
            lock = load_workflow_lock(expected_workflow_lock)
        except WorkflowLockError as exc:
            raise BundleError(str(exc)) from exc
        expected = {
            "rnaseq": asdict(lock.rnaseq),
            "differential": asdict(lock.differential),
        }
        if payload.get("workflows") != expected:
            raise BundleError("bundle workflow revisions do not match the workflow lock")
        if payload.get("runtime") != asdict(lock.runtime):
            raise BundleError("bundle runtime versions do not match the workflow lock")

    return {
        "verified": True,
        "bundle_status": "sealed",
        "artifact_count": len(records),
        "total_bytes": total_bytes,
        "components": checked_components,
        "workflows": payload.get("workflows"),
        "runtime": payload.get("runtime"),
    }
