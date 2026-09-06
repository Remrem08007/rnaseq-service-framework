"""Create and verify checksum-bound FASTQ input manifests."""

from __future__ import annotations

import csv
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

from .preflight import PreflightReport, validate_samplesheet


InputProgress = Callable[[str, int, int], None]


class InputManifestError(ValueError):
    """Raised when a FASTQ manifest cannot be created or verified."""


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _fastq_rows(samplesheet: Path) -> list[dict[str, object]]:
    report = PreflightReport()
    validate_samplesheet(samplesheet, report, check_files=True)
    if not report.valid:
        raise InputManifestError("samplesheet validation failed: " + "; ".join(report.errors))
    records: list[dict[str, object]] = []
    seen: set[Path] = set()
    try:
        with samplesheet.open("r", encoding="utf-8-sig", newline="") as handle:
            for line_number, row in enumerate(csv.DictReader(handle), start=2):
                for read in ("fastq_1", "fastq_2"):
                    value = (row.get(read) or "").strip()
                    if not value:
                        continue
                    source = Path(os.path.expandvars(os.path.expanduser(value)))
                    if not source.is_absolute():
                        source = samplesheet.parent / source
                    resolved = source.resolve(strict=True)
                    if resolved in seen:
                        raise InputManifestError(
                            f"resolved FASTQ appears more than once: {resolved}"
                        )
                    seen.add(resolved)
                    records.append(
                        {
                            "sample": (row.get("sample") or "").strip(),
                            "read": read,
                            "samplesheet_line": line_number,
                            "path": resolved,
                        }
                    )
    except (OSError, csv.Error) as exc:
        raise InputManifestError(f"could not read samplesheet: {exc}") from exc
    return records


def _hash_stable_file(
    path: Path,
    *,
    completed_before: int,
    total_bytes: int,
    progress: InputProgress | None,
) -> tuple[dict[str, int], str]:
    before = path.stat()
    digest = hashlib.sha256()
    current = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
            current += len(chunk)
            if progress is not None:
                progress(path.name, completed_before + current, total_bytes)
    after = path.stat()
    identity_before = {"size_bytes": before.st_size, "mtime_ns": before.st_mtime_ns}
    identity_after = {"size_bytes": after.st_size, "mtime_ns": after.st_mtime_ns}
    if identity_before != identity_after:
        raise InputManifestError(f"FASTQ changed while it was being hashed: {path}")
    if progress is not None and before.st_size == 0:
        progress(path.name, completed_before, total_bytes)
    return identity_after, digest.hexdigest()


def create_input_manifest(
    *,
    samplesheet: Path,
    output: Path,
    progress: InputProgress | None = None,
) -> dict[str, object]:
    """Hash every planned FASTQ once and exclusively write its manifest."""

    if output.exists():
        raise InputManifestError(f"refusing to overwrite existing manifest: {output}")
    samplesheet_resolved = samplesheet.resolve(strict=True)
    rows = _fastq_rows(samplesheet_resolved)
    total_bytes = sum(Path(row["path"]).stat().st_size for row in rows)
    completed = 0
    files = []
    for row in rows:
        path = Path(row["path"])
        identity, checksum = _hash_stable_file(
            path,
            completed_before=completed,
            total_bytes=total_bytes,
            progress=progress,
        )
        completed += identity["size_bytes"]
        files.append(
            {
                **{key: value for key, value in row.items() if key != "path"},
                "path": str(path),
                **identity,
                "sha256": checksum,
            }
        )
    payload: dict[str, object] = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "samplesheet": {
            "path": str(samplesheet_resolved),
            "size_bytes": samplesheet_resolved.stat().st_size,
            "sha256": _sha256_file(samplesheet_resolved),
        },
        "n_fastq_files": len(files),
        "total_bytes": total_bytes,
        "files": files,
        "contains_client_data": True,
        "contains_secrets": False,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        descriptor = os.open(output, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    except FileExistsError as exc:
        raise InputManifestError(f"refusing to overwrite existing manifest: {output}") from exc
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
    except Exception:
        output.unlink(missing_ok=True)
        raise
    return payload


def verify_input_manifest(
    manifest: Path,
    *,
    expected_samplesheet: Path | None = None,
    verify_hashes: bool = False,
) -> dict[str, object]:
    """Verify manifest structure and current FASTQ identity metadata."""

    try:
        payload = json.loads(manifest.resolve(strict=True).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise InputManifestError(f"could not read input manifest: {exc}") from exc
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise InputManifestError("unsupported input manifest schema")
    samplesheet_record = payload.get("samplesheet")
    if not isinstance(samplesheet_record, dict):
        raise InputManifestError("input manifest has no samplesheet identity")
    samplesheet = Path(str(samplesheet_record.get("path", ""))).resolve(strict=True)
    if expected_samplesheet is not None and samplesheet != expected_samplesheet.resolve(strict=True):
        raise InputManifestError("input manifest belongs to a different samplesheet")
    if samplesheet.stat().st_size != samplesheet_record.get("size_bytes"):
        raise InputManifestError("input-manifest samplesheet size has changed")
    if _sha256_file(samplesheet) != samplesheet_record.get("sha256"):
        raise InputManifestError("input-manifest samplesheet checksum has changed")
    files = payload.get("files")
    if not isinstance(files, list) or not files:
        raise InputManifestError("input manifest contains no FASTQ files")
    if payload.get("n_fastq_files") != len(files):
        raise InputManifestError("input manifest FASTQ count is inconsistent")
    seen: set[Path] = set()
    total = 0
    for index, record in enumerate(files, start=1):
        if not isinstance(record, dict):
            raise InputManifestError(f"invalid FASTQ record {index}")
        path = Path(str(record.get("path", ""))).resolve(strict=True)
        if path in seen:
            raise InputManifestError(f"duplicate FASTQ path in manifest: {path}")
        seen.add(path)
        metadata = path.stat()
        if metadata.st_size != record.get("size_bytes"):
            raise InputManifestError(f"FASTQ size changed after manifest creation: {path}")
        if metadata.st_mtime_ns != record.get("mtime_ns"):
            raise InputManifestError(f"FASTQ modification time changed after manifest creation: {path}")
        checksum = record.get("sha256")
        if not isinstance(checksum, str) or len(checksum) != 64:
            raise InputManifestError(f"invalid FASTQ checksum in manifest: {path}")
        if verify_hashes and _sha256_file(path) != checksum:
            raise InputManifestError(f"FASTQ checksum changed after manifest creation: {path}")
        total += metadata.st_size
    if payload.get("total_bytes") != total:
        raise InputManifestError("input manifest total byte count is inconsistent")
    return payload
