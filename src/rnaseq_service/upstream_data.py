"""Stage and verify public inputs used by the pinned nf-core smoke profiles."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import urllib.error
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Callable

MANIFEST_NAME = "upstream_test_data.manifest.json"
NETWORK_MODES = {"auto", "direct", "proxy"}
PROXY_VARIABLES = (
    "HTTPS_PROXY",
    "https_proxy",
    "HTTP_PROXY",
    "http_proxy",
    "NO_PROXY",
    "no_proxy",
)
ProgressCallback = Callable[[str, int, int], None]

RNASEQ_PROFILE_COMMIT = "626c8fab639062eade4b10747e919341cbf9b41a"
RNASEQ_DATA_COMMIT = "e07c1b158d1c4c9ea7978959d31e651098bec581"
RNASEQ_KRAKEN_COMMIT = "eb0cbf73c3f103f8aeda9878ba200e92b4d045d8"
MODULES_DATA_COMMIT = "d3b81d0289e70db899a7c5b8ab01730a2da7b435"
DIFFERENTIAL_DATA_COMMIT = "96f0de878b5a5758e0ab5d0a461a6a00c5bf24f3"


class UpstreamDataError(ValueError):
    """Raised when public smoke data cannot be staged or verified safely."""


@dataclass(frozen=True)
class Asset:
    role: str
    relative_path: str
    url: str


def _raw(commit: str, path: str) -> str:
    return f"https://raw.githubusercontent.com/nf-core/test-datasets/{commit}/{path}"


_RNASEQ_FASTQ_NAMES = (
    "SRR6357070_1.fastq.gz",
    "SRR6357070_2.fastq.gz",
    "SRR6357071_1.fastq.gz",
    "SRR6357071_2.fastq.gz",
    "SRR6357072_1.fastq.gz",
    "SRR6357072_2.fastq.gz",
    "SRR6357073_1.fastq.gz",
    "SRR6357074_1.fastq.gz",
    "SRR6357075_1.fastq.gz",
    "SRR6357076_1.fastq.gz",
    "SRR6357076_2.fastq.gz",
)

ASSETS = (
    Asset(
        "rnaseq_source_samplesheet",
        "rnaseq/source/samplesheet_test.csv",
        _raw(RNASEQ_PROFILE_COMMIT, "samplesheet/v3.10/samplesheet_test.csv"),
    ),
    Asset(
        "rnaseq_fasta",
        "rnaseq/reference/genome.fasta",
        _raw(RNASEQ_PROFILE_COMMIT, "reference/genome.fasta"),
    ),
    Asset(
        "rnaseq_gtf",
        "rnaseq/reference/genes_with_empty_tid.gtf.gz",
        _raw(RNASEQ_PROFILE_COMMIT, "reference/genes_with_empty_tid.gtf.gz"),
    ),
    Asset(
        "rnaseq_gff",
        "rnaseq/reference/genes.gff.gz",
        _raw(RNASEQ_PROFILE_COMMIT, "reference/genes.gff.gz"),
    ),
    Asset(
        "rnaseq_transcript_fasta",
        "rnaseq/reference/transcriptome.fasta",
        _raw(RNASEQ_PROFILE_COMMIT, "reference/transcriptome.fasta"),
    ),
    Asset(
        "rnaseq_additional_fasta",
        "rnaseq/reference/gfp.fa.gz",
        _raw(RNASEQ_PROFILE_COMMIT, "reference/gfp.fa.gz"),
    ),
    Asset(
        "rnaseq_source_bbsplit_list",
        "rnaseq/source/bbsplit_fasta_list.txt",
        _raw(RNASEQ_PROFILE_COMMIT, "reference/bbsplit_fasta_list.txt"),
    ),
    Asset(
        "rnaseq_hisat2_index",
        "rnaseq/reference/hisat2.tar.gz",
        _raw(RNASEQ_PROFILE_COMMIT, "reference/hisat2.tar.gz"),
    ),
    Asset(
        "rnaseq_salmon_index",
        "rnaseq/reference/salmon.tar.gz",
        _raw(RNASEQ_PROFILE_COMMIT, "reference/salmon.tar.gz"),
    ),
    Asset(
        "rnaseq_kraken_db",
        "rnaseq/reference/kraken2.tar.gz",
        _raw(
            RNASEQ_KRAKEN_COMMIT,
            "data/genomics/sarscov2/genome/db/kraken2.tar.gz",
        ),
    ),
    *(
        Asset(
            f"rnaseq_fastq:{name}",
            f"rnaseq/fastq/{name}",
            _raw(RNASEQ_DATA_COMMIT, f"testdata/GSE110004/{name}"),
        )
        for name in _RNASEQ_FASTQ_NAMES
    ),
    Asset(
        "rnaseq_bbsplit_sarscov2",
        "rnaseq/reference/GCA_009858895.3_ASM985889v3_genomic.200409.fna",
        _raw(RNASEQ_DATA_COMMIT, "reference/GCA_009858895.3_ASM985889v3_genomic.200409.fna"),
    ),
    Asset(
        "rnaseq_bbsplit_human",
        "rnaseq/reference/chr22_23800000-23980000.fa",
        _raw(RNASEQ_DATA_COMMIT, "reference/chr22_23800000-23980000.fa"),
    ),
    Asset(
        "differential_input",
        "differential/SRP254919.samplesheet.csv",
        _raw(
            MODULES_DATA_COMMIT,
            "data/genomics/mus_musculus/rnaseq_expression/SRP254919.samplesheet.csv",
        ),
    ),
    Asset(
        "differential_matrix",
        "differential/SRP254919.salmon.merged.gene_counts.top1000cov.tsv",
        _raw(
            MODULES_DATA_COMMIT,
            "data/genomics/mus_musculus/rnaseq_expression/"
            "SRP254919.salmon.merged.gene_counts.top1000cov.tsv",
        ),
    ),
    Asset(
        "differential_feature_length_matrix",
        "differential/SRP254919.spoofed_lengths.tsv",
        _raw(
            MODULES_DATA_COMMIT,
            "data/genomics/mus_musculus/rnaseq_expression/SRP254919.spoofed_lengths.tsv",
        ),
    ),
    Asset(
        "differential_contrasts",
        "differential/SRP254919.contrasts.yaml",
        _raw(DIFFERENTIAL_DATA_COMMIT, "testdata/SRP254919.contrasts.yaml"),
    ),
    Asset(
        "differential_gtf",
        "differential/Mus_musculus.GRCm38.81.gtf.gz",
        "https://ftp.ensembl.org/pub/release-81/gtf/mus_musculus/Mus_musculus.GRCm38.81.gtf.gz",
    ),
    Asset(
        "differential_gene_sets_files",
        "differential/mh.all.v2022.1.Mm.symbols.gmt",
        _raw(
            MODULES_DATA_COMMIT,
            "data/genomics/mus_musculus/gene_set_analysis/"
            "mh.all.v2022.1.Mm.symbols.gmt",
        ),
    ),
    Asset(
        "differential_decoupler_network",
        "differential/mouse_network.tsv",
        _raw(DIFFERENTIAL_DATA_COMMIT, "modules_testdata/progeny/mouse_network.tsv"),
    ),
)


def _safe_relative(value: object) -> PurePosixPath:
    if not isinstance(value, str) or not value:
        raise UpstreamDataError("test-data paths must be non-empty strings")
    path = PurePosixPath(value)
    if path.is_absolute() or ".." in path.parts or str(path) in {"", "."}:
        raise UpstreamDataError(f"unsafe test-data path: {value!r}")
    return path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _opener(network_mode: str) -> urllib.request.OpenerDirector:
    if network_mode == "direct":
        return urllib.request.build_opener(urllib.request.ProxyHandler({}))
    return urllib.request.build_opener()


def _proxy_present() -> bool:
    return any(os.environ.get(name) for name in PROXY_VARIABLES)


def _redact(text: str) -> str:
    rendered = text
    secrets = [os.environ[name] for name in PROXY_VARIABLES if os.environ.get(name)]
    for value in sorted(set(secrets), key=len, reverse=True):
        rendered = rendered.replace(value, "<redacted-proxy>")
    return rendered


def select_network_mode(requested: str, *, timeout_seconds: int = 10) -> str:
    """Select direct first, then a configured proxy, without downloading data."""

    if requested not in NETWORK_MODES:
        raise UpstreamDataError("network mode must be auto, direct, or proxy")
    if timeout_seconds < 1:
        raise UpstreamDataError("network timeout must be at least 1 second")
    if requested == "proxy":
        if not _proxy_present():
            raise UpstreamDataError("proxy mode requires a standard proxy environment variable")
        return "proxy"
    if requested == "direct":
        return "direct"
    probe = urllib.request.Request(
        ASSETS[0].url,
        method="HEAD",
        headers={"User-Agent": "rnaseq-service-framework"},
    )
    try:
        with _opener("direct").open(probe, timeout=timeout_seconds):
            return "direct"
    except (urllib.error.URLError, OSError):
        if _proxy_present():
            try:
                with _opener("proxy").open(probe, timeout=timeout_seconds):
                    return "proxy"
            except (urllib.error.URLError, OSError):
                pass
    raise UpstreamDataError(
        "no direct HTTPS access or working configured proxy; stage on a connected host"
    )


def _download(
    asset: Asset,
    target: Path,
    *,
    opener: urllib.request.OpenerDirector,
    timeout_seconds: int,
    progress: ProgressCallback | None,
) -> dict[str, object]:
    if target.exists():
        raise UpstreamDataError(f"refusing to overwrite staged asset: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.name + ".part")
    if temporary.exists():
        raise UpstreamDataError(f"partial download already exists: {temporary}")
    request = urllib.request.Request(
        asset.url, headers={"User-Agent": "rnaseq-service-framework"}
    )
    digest = hashlib.sha256()
    downloaded = 0
    try:
        with opener.open(request, timeout=timeout_seconds) as response, temporary.open(
            "xb"
        ) as handle:
            total = int(response.headers.get("Content-Length", "0"))
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)
                digest.update(chunk)
                downloaded += len(chunk)
                if progress is not None:
                    progress(asset.role, downloaded, total)
        if downloaded == 0:
            raise UpstreamDataError(f"downloaded asset is empty: {asset.url}")
        if total > 0 and downloaded != total:
            raise UpstreamDataError(
                f"incomplete download for {asset.role}: expected {total} bytes, got {downloaded}"
            )
        temporary.replace(target)
    except (urllib.error.URLError, OSError) as exc:
        temporary.unlink(missing_ok=True)
        raise UpstreamDataError(
            f"download failed for {asset.role}: {_redact(str(exc))}"
        ) from exc
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    if progress is not None:
        progress(asset.role, downloaded, downloaded)
    return {
        "role": asset.role,
        "path": asset.relative_path,
        "source_url": asset.url,
        "size_bytes": downloaded,
        "sha256": digest.hexdigest(),
    }


def _write_rnaseq_controls(root: Path) -> tuple[Path, Path]:
    source = root / "rnaseq/source/samplesheet_test.csv"
    output = root / "rnaseq/samplesheet.csv"
    fastq_root = root / "rnaseq/fastq"
    with source.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != ["sample", "fastq_1", "fastq_2", "strandedness"]:
            raise UpstreamDataError("upstream RNA-seq samplesheet columns changed")
        rows = list(reader)
    with output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=reader.fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            for field in ("fastq_1", "fastq_2"):
                if row[field]:
                    local = fastq_root / Path(row[field]).name
                    if not local.is_file():
                        raise UpstreamDataError(f"samplesheet FASTQ was not staged: {local}")
                    row[field] = local.relative_to(root).as_posix()
            writer.writerow(row)

    bbsplit = root / "rnaseq/bbsplit_fasta_list.csv"
    reference = root / "rnaseq/reference"
    with bbsplit.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "sarscov2",
                (
                    reference / "GCA_009858895.3_ASM985889v3_genomic.200409.fna"
                ).relative_to(root).as_posix(),
            ]
        )
        writer.writerow(
            [
                "human",
                (reference / "chr22_23800000-23980000.fa")
                .relative_to(root)
                .as_posix(),
            ]
        )
    return output, bbsplit


def _parameter_paths(root: Path, samplesheet: Path, bbsplit: Path) -> dict[str, dict[str, str]]:
    rnaseq_reference = Path("rnaseq/reference")
    differential = Path("differential")
    return {
        "rnaseq": {
            "input": samplesheet.relative_to(root).as_posix(),
            "fasta": str(rnaseq_reference / "genome.fasta"),
            "gtf": str(rnaseq_reference / "genes_with_empty_tid.gtf.gz"),
            "gff": str(rnaseq_reference / "genes.gff.gz"),
            "transcript_fasta": str(rnaseq_reference / "transcriptome.fasta"),
            "additional_fasta": str(rnaseq_reference / "gfp.fa.gz"),
            "bbsplit_fasta_list": bbsplit.relative_to(root).as_posix(),
            "hisat2_index": str(rnaseq_reference / "hisat2.tar.gz"),
            "salmon_index": str(rnaseq_reference / "salmon.tar.gz"),
            "kraken_db": str(rnaseq_reference / "kraken2.tar.gz"),
        },
        "differential": {
            "input": str(differential / "SRP254919.samplesheet.csv"),
            "matrix": str(differential / "SRP254919.salmon.merged.gene_counts.top1000cov.tsv"),
            "feature_length_matrix": str(differential / "SRP254919.spoofed_lengths.tsv"),
            "contrasts": str(differential / "SRP254919.contrasts.yaml"),
            "gtf": str(differential / "Mus_musculus.GRCm38.81.gtf.gz"),
            "gene_sets_files": str(differential / "mh.all.v2022.1.Mm.symbols.gmt"),
            "decoupler_network": str(differential / "mouse_network.tsv"),
        },
    }


def _verify_local_controls(root: Path) -> None:
    source_path = root / "rnaseq/source/samplesheet_test.csv"
    local_path = root / "rnaseq/samplesheet.csv"
    try:
        with source_path.open("r", encoding="utf-8", newline="") as handle:
            source_reader = csv.DictReader(handle)
            source_rows = list(source_reader)
        with local_path.open("r", encoding="utf-8", newline="") as handle:
            local_reader = csv.DictReader(handle)
            local_rows = list(local_reader)
    except (OSError, UnicodeError, csv.Error) as exc:
        raise UpstreamDataError(f"could not read RNA-seq smoke controls: {exc}") from exc
    fields = ["sample", "fastq_1", "fastq_2", "strandedness"]
    if source_reader.fieldnames != fields or local_reader.fieldnames != fields:
        raise UpstreamDataError("RNA-seq smoke control columns changed")
    if len(source_rows) != len(local_rows) or not source_rows:
        raise UpstreamDataError("RNA-seq local samplesheet rows changed")
    for source, local in zip(source_rows, local_rows, strict=True):
        if (
            source["sample"] != local["sample"]
            or source["strandedness"] != local["strandedness"]
        ):
            raise UpstreamDataError("RNA-seq local samplesheet metadata changed")
        for field in ("fastq_1", "fastq_2"):
            expected = (
                f"rnaseq/fastq/{Path(source[field]).name}" if source[field] else ""
            )
            if local[field] != expected:
                raise UpstreamDataError("RNA-seq local samplesheet paths changed")
    bbsplit_path = root / "rnaseq/bbsplit_fasta_list.csv"
    try:
        with bbsplit_path.open("r", encoding="utf-8", newline="") as handle:
            bbsplit_rows = list(csv.reader(handle))
    except (OSError, UnicodeError, csv.Error) as exc:
        raise UpstreamDataError(f"could not read local BBsplit control: {exc}") from exc
    if bbsplit_rows != [
        [
            "sarscov2",
            "rnaseq/reference/GCA_009858895.3_ASM985889v3_genomic.200409.fna",
        ],
        ["human", "rnaseq/reference/chr22_23800000-23980000.fa"],
    ]:
        raise UpstreamDataError("RNA-seq local BBsplit paths changed")


def stage_upstream_test_data(
    *,
    output_root: Path,
    network_mode: str = "auto",
    timeout_seconds: int = 60,
    progress: ProgressCallback | None = None,
) -> dict[str, object]:
    """Download all pinned public inputs and write local-only smoke controls."""

    if output_root.exists():
        raise UpstreamDataError(f"refusing existing test-data root: {output_root}")
    if timeout_seconds < 1:
        raise UpstreamDataError("network timeout must be at least 1 second")
    selected = select_network_mode(network_mode, timeout_seconds=min(timeout_seconds, 10))
    root = output_root.resolve()
    root.mkdir(parents=True, exist_ok=False)
    records: list[dict[str, object]] = []
    try:
        opener = _opener(selected)
        for asset in ASSETS:
            records.append(
                _download(
                    asset,
                    root / asset.relative_path,
                    opener=opener,
                    timeout_seconds=timeout_seconds,
                    progress=progress,
                )
            )
        samplesheet, bbsplit = _write_rnaseq_controls(root)
        for role, path in (
            ("rnaseq_local_samplesheet", samplesheet),
            ("rnaseq_local_bbsplit_list", bbsplit),
        ):
            records.append(
                {
                    "role": role,
                    "path": path.relative_to(root).as_posix(),
                    "source_url": None,
                    "size_bytes": path.stat().st_size,
                    "sha256": _sha256(path),
                }
            )
        parameters = _parameter_paths(root, samplesheet, bbsplit)
        payload: dict[str, object] = {
            "schema_version": 1,
            "created_at": datetime.now(timezone.utc).isoformat(),
            "status": "complete",
            "requested_network_mode": network_mode,
            "selected_network_mode": selected,
            "source_revisions": {
                "rnaseq_profile": RNASEQ_PROFILE_COMMIT,
                "rnaseq_data": RNASEQ_DATA_COMMIT,
                "rnaseq_kraken": RNASEQ_KRAKEN_COMMIT,
                "modules_data": MODULES_DATA_COMMIT,
                "differential_data": DIFFERENTIAL_DATA_COMMIT,
            },
            "parameters": parameters,
            "artifacts": records,
            "artifact_count": len(records),
            "total_bytes": sum(int(record["size_bytes"]) for record in records),
            "contains_client_data": False,
            "contains_secrets": False,
        }
        manifest = root / MANIFEST_NAME
        manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        os.chmod(manifest, 0o600)
        return verify_upstream_test_data(root)
    except Exception:
        # Keep partial downloads for diagnosis; a new destination is required for retry.
        raise


def verify_upstream_test_data(root: Path) -> dict[str, object]:
    """Verify exact membership, hashes, and local-only parameter paths."""

    resolved = root.resolve(strict=True)
    manifest = resolved / MANIFEST_NAME
    try:
        payload = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise UpstreamDataError(f"could not read upstream test-data manifest: {exc}") from exc
    if (
        not isinstance(payload, dict)
        or payload.get("schema_version") != 1
        or payload.get("status") != "complete"
        or payload.get("contains_client_data") is not False
        or payload.get("contains_secrets") is not False
    ):
        raise UpstreamDataError("unsupported or unsafe upstream test-data manifest")
    artifacts = payload.get("artifacts")
    if not isinstance(artifacts, list) or not artifacts:
        raise UpstreamDataError("upstream test-data manifest has no artifacts")
    records: dict[str, dict[str, object]] = {}
    for record in artifacts:
        if not isinstance(record, dict):
            raise UpstreamDataError("invalid upstream test-data artifact record")
        relative = _safe_relative(record.get("path")).as_posix()
        if relative in records or relative == MANIFEST_NAME:
            raise UpstreamDataError(f"duplicate or reserved test-data path: {relative}")
        records[relative] = record
    paths = list(resolved.rglob("*"))
    links = [path for path in paths if path.is_symlink()]
    if links:
        raise UpstreamDataError(f"symbolic links are not allowed: {links[0]}")
    actual = {
        path.relative_to(resolved).as_posix()
        for path in paths
        if path.is_file() and path != manifest
    }
    if actual != set(records):
        raise UpstreamDataError(
            f"upstream test-data membership changed; missing={sorted(set(records) - actual)}, "
            f"unexpected={sorted(actual - set(records))}"
        )
    expected_sources = {asset.relative_path: (asset.role, asset.url) for asset in ASSETS}
    expected_sources.update(
        {
            "rnaseq/samplesheet.csv": ("rnaseq_local_samplesheet", None),
            "rnaseq/bbsplit_fasta_list.csv": ("rnaseq_local_bbsplit_list", None),
        }
    )
    if set(records) != set(expected_sources):
        raise UpstreamDataError("upstream test-data artifact catalog changed")
    total_bytes = 0
    for relative, record in records.items():
        path = resolved / relative
        if path.stat().st_size != record.get("size_bytes"):
            raise UpstreamDataError(f"upstream test-data size changed: {relative}")
        if _sha256(path) != record.get("sha256"):
            raise UpstreamDataError(f"upstream test-data checksum changed: {relative}")
        expected_role, expected_url = expected_sources[relative]
        if record.get("role") != expected_role or record.get("source_url") != expected_url:
            raise UpstreamDataError(f"upstream test-data provenance changed: {relative}")
        total_bytes += path.stat().st_size
    if payload.get("artifact_count") != len(records) or payload.get("total_bytes") != total_bytes:
        raise UpstreamDataError("upstream test-data totals are invalid")
    expected_revisions = {
        "rnaseq_profile": RNASEQ_PROFILE_COMMIT,
        "rnaseq_data": RNASEQ_DATA_COMMIT,
        "rnaseq_kraken": RNASEQ_KRAKEN_COMMIT,
        "modules_data": MODULES_DATA_COMMIT,
        "differential_data": DIFFERENTIAL_DATA_COMMIT,
    }
    if payload.get("source_revisions") != expected_revisions:
        raise UpstreamDataError("upstream test-data source revisions changed")
    for relative in ("rnaseq/samplesheet.csv", "rnaseq/bbsplit_fasta_list.csv"):
        try:
            text = (resolved / relative).read_text(encoding="utf-8")
        except (OSError, UnicodeError) as exc:
            raise UpstreamDataError(f"could not read local control {relative}: {exc}") from exc
        if "://" in text or str(resolved) in text:
            raise UpstreamDataError(f"non-portable path remains in local control: {relative}")
    _verify_local_controls(resolved)
    parameters = payload.get("parameters")
    if not isinstance(parameters, dict) or set(parameters) != {"rnaseq", "differential"}:
        raise UpstreamDataError("upstream test-data parameters are missing")
    checked: dict[str, dict[str, str]] = {}
    for stage, values in parameters.items():
        if not isinstance(values, dict) or not values:
            raise UpstreamDataError(f"upstream test-data parameters are invalid: {stage}")
        checked[stage] = {}
        for name, value in values.items():
            try:
                relative = _safe_relative(value)
                path = (resolved / relative).resolve(strict=True)
                path.relative_to(resolved)
            except (OSError, ValueError) as exc:
                raise UpstreamDataError(
                    f"upstream test-data parameter {stage}.{name} is not local to the bundle"
                ) from exc
            if not path.is_file():
                raise UpstreamDataError(f"upstream test-data parameter is not a file: {path}")
            checked[stage][str(name)] = str(path)
    return {
        "verified": True,
        "manifest": str(manifest),
        "artifact_count": len(records),
        "total_bytes": total_bytes,
        "parameters": checked,
        "source_revisions": payload.get("source_revisions"),
        "contains_client_data": False,
        "contains_secrets": False,
    }
