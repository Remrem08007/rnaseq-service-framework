"""CLI for sealing completed nf-core/rnaseq primary-processing runs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .primary import (
    PrimaryDeliveryError,
    create_primary_receipt,
    diagnose_primary_run,
    prepare_safe_restart,
)


class ArtifactProgress:
    def __init__(self) -> None:
        self.role: str | None = None
        self.bucket = -1

    def __call__(self, role: str, current: int, total: int) -> None:
        percent = 100 if total == 0 else min(100, int(current * 100 / total))
        bucket = percent // 5
        if role == self.role and bucket == self.bucket:
            return
        self.role = role
        self.bucket = bucket
        width = 20
        filled = min(width, int(percent * width / 100))
        bar = "#" * filled + "-" * (width - filled)
        print(
            f"\rHashing {role:<22} [{bar}] {percent:3d}%",
            end="\n" if percent == 100 else "",
            file=sys.stderr,
            flush=True,
        )


class LogProgress:
    def __init__(self) -> None:
        self.last: dict[str, int] = {}

    def __call__(self, source: str, line: int, done: bool) -> None:
        if not done and line != 1 and line % 1000:
            return
        if done and self.last.get(source) == line:
            print(file=sys.stderr, flush=True)
            return
        print(
            f"\rScanning {source}: {line} lines",
            end="\n" if done else "",
            file=sys.stderr,
            flush=True,
        )
        self.last[source] = line


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rnaseq-service-primary",
        description="Validate primary outputs and write an immutable completion receipt.",
    )
    commands = parser.add_subparsers(dest="command", required=True)
    complete = commands.add_parser("complete", help="validate and seal completed outputs")
    complete.add_argument("--run-plan", type=Path, required=True)
    complete.add_argument("--output", type=Path, required=True)
    complete.add_argument("--quiet", action="store_true")
    diagnose = commands.add_parser("diagnose", help="classify an incomplete run")
    diagnose.add_argument("--run-plan", type=Path, required=True)
    diagnose.add_argument("--submission-receipt", type=Path)
    diagnose.add_argument("--nextflow-log", type=Path)
    diagnose.add_argument("--quiet", action="store_true")
    restart = commands.add_parser("restart", help="prepare a reviewed -resume launcher")
    restart.add_argument("--run-plan", type=Path, required=True)
    restart.add_argument("--settings", type=Path, required=True)
    restart.add_argument("--previous-receipt", type=Path, required=True)
    restart.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "complete":
            result = create_primary_receipt(
                run_plan=args.run_plan,
                output=args.output,
                progress=None if args.quiet else ArtifactProgress(),
            )
        elif args.command == "diagnose":
            result = diagnose_primary_run(
                run_plan=args.run_plan,
                submission_receipt=args.submission_receipt,
                nextflow_log=args.nextflow_log,
                progress=None if args.quiet else LogProgress(),
            )
        else:
            result = prepare_safe_restart(
                run_plan=args.run_plan,
                settings_path=args.settings,
                previous_receipt=args.previous_receipt,
                output=args.output,
            )
    except (PrimaryDeliveryError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
