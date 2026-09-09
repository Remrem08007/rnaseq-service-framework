"""CLI for deterministic synthetic RNA-seq validation studies."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .synthetic import SyntheticStudyError, create_synthetic_study
from .synthetic_e2e import SyntheticContractError, run_contract_validation


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


class StageProgress:
    def __call__(self, stage: str, current: int, total: int) -> None:
        width = 20
        filled = int(current * width / total)
        bar = "#" * filled + "-" * (width - filled)
        print(
            f"\rValidating [{bar}] {current}/{total} {stage}",
            end="\n" if current == total else "",
            file=sys.stderr,
            flush=True,
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rnaseq-service-synthetic",
        description="Create and validate deterministic non-client RNA-seq fixtures.",
    )
    commands = parser.add_subparsers(dest="command", required=True)
    create = commands.add_parser("create", help="create paired FASTQs and known truth")
    create.add_argument("--output-dir", type=Path, required=True)
    create.add_argument("--seed", type=int, default=20260909)
    create.add_argument("--quiet", action="store_true")
    validate = commands.add_parser(
        "validate-contract", help="exercise all wrapper gates with synthetic artifacts"
    )
    validate.add_argument("--output-dir", type=Path, required=True)
    validate.add_argument("--workflow-lock", type=Path, required=True)
    validate.add_argument("--quiet", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "create":
            result = create_synthetic_study(
                output_dir=args.output_dir,
                seed=args.seed,
                progress=None if args.quiet else SampleProgress(),
            )
            summary = {
                "dataset": result["dataset"],
                "n_samples": len(result["samples"]),
                "n_fastq_files": len(result["files"]),
                "truth": str((args.output_dir.resolve() / "truth.json")),
            }
        else:
            result = run_contract_validation(
                output_dir=args.output_dir,
                workflow_lock=args.workflow_lock,
                progress=None if args.quiet else StageProgress(),
            )
            summary = {
                "status": result["stage"],
                "scientific_execution_performed": result["scientific_execution_performed"],
                "n_samples": result["n_samples"],
                "n_genes": result["n_genes"],
                "direction_validation_passed": result["direction_validation"]["passed"],
                "receipt": str(
                    args.output_dir.resolve() / "synthetic_contract_receipt.json"
                ),
            }
    except (SyntheticStudyError, SyntheticContractError, OSError, ValueError) as exc:
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
