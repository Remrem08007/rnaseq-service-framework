"""CLI for deterministic synthetic RNA-seq validation studies."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .synthetic import SyntheticStudyError, create_synthetic_study


class SampleProgress:
    def __init__(self) -> None:
        self.last = 0

    def __call__(self, sample: str, current: int, total: int) -> None:
        if current == self.last:
            return
        self.last = current
        percent = int(current * 100 / total)
        width = 20
        filled = int(current * width / total)
        bar = "#" * filled + "-" * (width - filled)
        print(
            f"\rGenerating [{bar}] {percent:3d}% ({current}/{total}) {sample}",
            end="\n" if current == total else "",
            file=sys.stderr,
            flush=True,
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rnaseq-service-synthetic",
        description="Create a deterministic non-client paired-end validation study.",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260909)
    parser.add_argument("--quiet", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = create_synthetic_study(
            output_dir=args.output_dir,
            seed=args.seed,
            progress=None if args.quiet else SampleProgress(),
        )
    except (SyntheticStudyError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "dataset": result["dataset"],
                "n_samples": len(result["samples"]),
                "n_fastq_files": len(result["files"]),
                "truth": str((args.output_dir.resolve() / "truth.json")),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

