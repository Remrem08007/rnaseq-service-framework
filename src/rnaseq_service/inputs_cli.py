"""CLI for checksum-bound RNA-seq FASTQ manifests."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .inputs import InputManifestError, create_input_manifest


class InputProgress:
    def __init__(self) -> None:
        self.bucket = -1

    def __call__(self, filename: str, current: int, total: int) -> None:
        percent = 100 if total == 0 else min(100, int(current * 100 / total))
        bucket = percent // 2
        if bucket == self.bucket:
            return
        self.bucket = bucket
        width = 25
        filled = min(width, int(percent * width / 100))
        bar = "#" * filled + "-" * (width - filled)
        print(
            f"\rHashing FASTQs [{bar}] {percent:3d}%  {filename[:30]:<30}",
            end="\n" if percent == 100 else "",
            file=sys.stderr,
            flush=True,
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rnaseq-service-inputs",
        description="Hash FASTQs and write an immutable input manifest.",
    )
    parser.add_argument("--samplesheet", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--quiet", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = create_input_manifest(
            samplesheet=args.samplesheet,
            output=args.output,
            progress=None if args.quiet else InputProgress(),
        )
    except (InputManifestError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
