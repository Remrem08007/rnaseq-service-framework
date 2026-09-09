from __future__ import annotations

import io
import json
from pathlib import Path

import pytest

import rnaseq_service.upstream_data as upstream_data
from rnaseq_service.upstream_data import (
    ASSETS,
    UpstreamDataError,
    stage_upstream_test_data,
    verify_upstream_test_data,
)
from rnaseq_service.upstream_data_cli import main as data_main


SAMPLESHEET = """sample,fastq_1,fastq_2,strandedness
WT_REP1,https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357070_1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357070_2.fastq.gz,auto
WT_REP1,https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357071_1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357071_2.fastq.gz,auto
WT_REP2,https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357072_1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357072_2.fastq.gz,reverse
RAP1_UNINDUCED_REP1,https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357073_1.fastq.gz,,reverse
RAP1_UNINDUCED_REP2,https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357074_1.fastq.gz,,reverse
RAP1_UNINDUCED_REP2,https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357075_1.fastq.gz,,reverse
RAP1_IAA_30M_REP1,https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357076_1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357076_2.fastq.gz,reverse
"""


class Response(io.BytesIO):
    def __init__(self, content: bytes) -> None:
        super().__init__(content)
        self.headers = {"Content-Length": str(len(content))}

    def __enter__(self) -> "Response":
        return self

    def __exit__(self, *args: object) -> None:
        self.close()


class FakeOpener:
    def open(self, request: object, timeout: int) -> Response:
        url = str(getattr(request, "full_url"))
        if url.endswith("samplesheet_test.csv"):
            return Response(SAMPLESHEET.encode())
        return Response((url + "\n").encode())


def _stage(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    monkeypatch.setattr(
        upstream_data, "select_network_mode", lambda *args, **kwargs: "direct"
    )
    monkeypatch.setattr(upstream_data, "_opener", lambda mode: FakeOpener())
    root = tmp_path / "upstream-test-data"
    result = stage_upstream_test_data(output_root=root)
    assert result["verified"] is True
    return root


def test_asset_catalog_is_pinned_and_complete() -> None:
    assert len(ASSETS) == 30
    assert len({asset.relative_path for asset in ASSETS}) == 30
    github_urls = [asset.url for asset in ASSETS if "github" in asset.url]
    assert all("/rnaseq/" not in url for url in github_urls)
    assert all("/modules/" not in url for url in github_urls)
    assert all("/differentialabundance/" not in url for url in github_urls)


def test_stage_rewrites_nested_urls_and_verifies(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = _stage(tmp_path, monkeypatch)
    samplesheet = (root / "rnaseq/samplesheet.csv").read_text(encoding="utf-8")
    bbsplit = (root / "rnaseq/bbsplit_fasta_list.csv").read_text(encoding="utf-8")
    assert "http" not in samplesheet
    assert "http" not in bbsplit
    assert "rnaseq/fastq/SRR6357070_1.fastq.gz" in samplesheet
    assert "rnaseq/reference/chr22_23800000-23980000.fa" in bbsplit
    assert str(root.resolve()) not in samplesheet
    assert str(root.resolve()) not in bbsplit
    result = verify_upstream_test_data(root)
    assert result["artifact_count"] == 32
    assert set(result["parameters"]) == {"rnaseq", "differential"}


def test_verify_rejects_modified_or_unexpected_data(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = _stage(tmp_path, monkeypatch)
    target = root / "differential/mouse_network.tsv"
    original = target.read_text(encoding="utf-8")
    target.write_text("changed\n", encoding="utf-8")
    with pytest.raises(UpstreamDataError, match="size changed|checksum changed"):
        verify_upstream_test_data(root)

    target.write_text(original, encoding="utf-8")
    (root / "unexpected.txt").write_text("x", encoding="utf-8")
    with pytest.raises(UpstreamDataError, match="membership changed"):
        verify_upstream_test_data(root)


def test_verify_rejects_changed_source_provenance(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = _stage(tmp_path, monkeypatch)
    manifest = root / "upstream_test_data.manifest.json"
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    payload["artifacts"][0]["source_url"] = "https://example.invalid/replacement"
    manifest.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(UpstreamDataError, match="provenance changed"):
        verify_upstream_test_data(root)


def test_cli_verify(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    root = _stage(tmp_path, monkeypatch)
    assert data_main(["verify", "--root", str(root)]) == 0
    assert '"verified": true' in capsys.readouterr().out


def test_stage_refuses_existing_root(tmp_path: Path) -> None:
    root = tmp_path / "existing"
    root.mkdir()
    with pytest.raises(UpstreamDataError, match="refusing existing"):
        stage_upstream_test_data(output_root=root)
