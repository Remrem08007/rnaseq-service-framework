"""Read and validate exact upstream workflow revisions."""

from __future__ import annotations

import re
import tomllib
from dataclasses import dataclass
from pathlib import Path


EXACT_REVISION = re.compile(r"^\d+\.\d+\.\d+(?:[-+][A-Za-z0-9.-]+)?$")


class WorkflowLockError(ValueError):
    """Raised when the upstream workflow lock is unsafe or incomplete."""


@dataclass(frozen=True)
class WorkflowSpec:
    name: str
    revision: str


@dataclass(frozen=True)
class RuntimeSpec:
    nextflow_version: str
    nf_core_tools_version: str


@dataclass(frozen=True)
class WorkflowLock:
    rnaseq: WorkflowSpec
    differential: WorkflowSpec
    runtime: RuntimeSpec
    last_reviewed: str


def _require_string(table: dict[str, object], key: str, path: Path) -> str:
    value = table.get(key)
    if not isinstance(value, str) or not value.strip():
        raise WorkflowLockError(f"{path}: {key!r} must be a non-empty string")
    return value.strip()


def _validate_spec(name: str, revision: str, path: Path) -> WorkflowSpec:
    if not name.startswith("nf-core/") or name.count("/") != 1:
        raise WorkflowLockError(
            f"{path}: workflow name {name!r} must use the nf-core/<pipeline> form"
        )
    if not EXACT_REVISION.fullmatch(revision):
        raise WorkflowLockError(
            f"{path}: revision {revision!r} is not an exact semantic version"
        )
    return WorkflowSpec(name=name, revision=revision)


def load_workflow_lock(path: Path) -> WorkflowLock:
    """Load a TOML lock and reject floating/development workflow revisions."""

    try:
        with path.open("rb") as handle:
            payload = tomllib.load(handle)
    except FileNotFoundError as exc:
        raise WorkflowLockError(f"{path}: workflow lock does not exist") from exc
    except (OSError, tomllib.TOMLDecodeError) as exc:
        raise WorkflowLockError(f"{path}: could not read workflow lock: {exc}") from exc

    workflows = payload.get("workflows")
    runtime = payload.get("runtime")
    policy = payload.get("policy")
    if (
        not isinstance(workflows, dict)
        or not isinstance(runtime, dict)
        or not isinstance(policy, dict)
    ):
        raise WorkflowLockError(
            f"{path}: [workflows], [runtime], and [policy] tables are required"
        )
    if policy.get("versions_are_pinned") is not True:
        raise WorkflowLockError(f"{path}: policy.versions_are_pinned must be true")
    if policy.get("allow_development_revisions") is not False:
        raise WorkflowLockError(
            f"{path}: policy.allow_development_revisions must be false"
        )

    rnaseq = _validate_spec(
        _require_string(workflows, "rnaseq_name", path),
        _require_string(workflows, "rnaseq_revision", path),
        path,
    )
    differential = _validate_spec(
        _require_string(workflows, "differential_name", path),
        _require_string(workflows, "differential_revision", path),
        path,
    )
    nextflow_version = _require_string(runtime, "nextflow_version", path)
    nf_core_tools_version = _require_string(runtime, "nf_core_tools_version", path)
    for label, version in (
        ("nextflow_version", nextflow_version),
        ("nf_core_tools_version", nf_core_tools_version),
    ):
        if not EXACT_REVISION.fullmatch(version):
            raise WorkflowLockError(
                f"{path}: runtime {label} {version!r} is not an exact semantic version"
            )
    return WorkflowLock(
        rnaseq=rnaseq,
        differential=differential,
        runtime=RuntimeSpec(
            nextflow_version=nextflow_version,
            nf_core_tools_version=nf_core_tools_version,
        ),
        last_reviewed=_require_string(policy, "last_reviewed", path),
    )
