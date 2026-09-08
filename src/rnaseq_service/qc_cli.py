"""CLI for configurable, non-automatic RNA-seq sample QC assessment."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .qc import QCError, evaluate_qc, finalize_qc


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
    commands = parser.add_subparsers(dest="command", required=True)
    assess = commands.add_parser("assess", help="evaluate MultiQC evidence and create a review template")
    assess.add_argument("--completion-receipt", type=Path, required=True)
    assess.add_argument("--policy", type=Path, required=True)
    assess.add_argument("--outdir", type=Path, required=True)
    assess.add_argument("--quiet", action="store_true")
    finalize = commands.add_parser("finalize", help="seal reviewed include/exclude decisions")
    finalize.add_argument("--assessment", type=Path, required=True)
    finalize.add_argument("--decisions", type=Path, required=True)
    finalize.add_argument("--outdir", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "assess":
            result = evaluate_qc(
                completion_receipt=args.completion_receipt,
                policy_path=args.policy,
                outdir=args.outdir,
                progress=None if args.quiet else MetricProgress(),
            )
            summary = {
                "status": result["status"],
                "sample_count": result["sample_count"],
                "flagged_sample_count": result["flagged_sample_count"],
                "automatic_exclusions": result["automatic_exclusions"],
            }
        else:
            result = finalize_qc(
                assessment_path=args.assessment,
                decisions_path=args.decisions,
                outdir=args.outdir,
            )
            summary = {
                "status": result["status"],
                "n_assessed": result["n_assessed"],
                "n_accepted": result["n_accepted"],
                "n_excluded": result["n_excluded"],
                "automatic_exclusions": result["automatic_exclusions"],
            }
    except (QCError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(
        json.dumps(
            summary,
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
