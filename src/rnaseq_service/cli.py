"""Command-line interface for service intake validation."""

from __future__ import annotations

import argparse
from pathlib import Path

from .preflight import run_preflight


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rnaseq-service-preflight",
        description="Validate RNA-seq sample, design, and contrast tables before launch.",
    )
    parser.add_argument("--samplesheet", type=Path, required=True)
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument("--contrasts", type=Path, required=True)
    parser.add_argument(
        "--check-files",
        action="store_true",
        help="Require every referenced FASTQ to exist on the current filesystem.",
    )
    parser.add_argument(
        "--min-replicates",
        type=int,
        default=2,
        help="Minimum samples required in each requested contrast level (default: 2).",
    )
    parser.add_argument(
        "--json-out",
        type=Path,
        help="Write the complete machine-readable report to this path.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.min_replicates < 1:
        raise SystemExit("--min-replicates must be at least 1")

    report = run_preflight(
        args.samplesheet,
        args.design,
        args.contrasts,
        check_files=args.check_files,
        min_replicates=args.min_replicates,
    )
    rendered = report.to_json()
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0 if report.valid else 2


if __name__ == "__main__":
    raise SystemExit(main())
