"""CLI for configurable, non-automatic RNA-seq sample QC assessment."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .qc import QCError, evaluate_qc


class MetricProgress:
    def __init__(self) -> None:
        self.last: tuple[str, int] | None = None

    def __call__(self, metric: str, current: int, total: int) -> None:
        state = (metric, current)
        if state == self.last:
            return
        self.last = state
        width = 20
        filled = width if total == 0 else int(width * current / total)
        bar = "#" * filled + "-" * (width - filled)
        print(
            f"\rEvaluating QC [{bar}] {current:>2}/{total:<2} {metric}",
            end="\n" if current == total else "",
            file=sys.stderr,
            flush=True,
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rnaseq-service-qc",
        description="Assess checksum-bound MultiQC evidence without automatically excluding samples.",
    )
    parser.add_argument("--completion-receipt", type=Path, required=True)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--quiet", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = evaluate_qc(
            completion_receipt=args.completion_receipt,
            policy_path=args.policy,
            outdir=args.outdir,
            progress=None if args.quiet else MetricProgress(),
        )
    except (QCError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "status": result["status"],
                "sample_count": result["sample_count"],
                "flagged_sample_count": result["flagged_sample_count"],
                "automatic_exclusions": result["automatic_exclusions"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
