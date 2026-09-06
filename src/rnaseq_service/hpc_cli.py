"""CLI for HPC configuration, launch, status, and resource inspection."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

from .hpc_config import HPCConfigError, write_nextflow_config
from .hpc_run import HPCRunError, prepare_launcher, query_job, submit_launcher
from .resources import ResourceSummaryError, write_resource_summary


class RowProgress:
    def __init__(self) -> None:
        self.last_bucket = -1

    def __call__(self, row: int, total: int) -> None:
        percent = 100 if total == 0 else min(100, int(row * 100 / total))
        bucket = percent // 5
        if bucket == self.last_bucket:
            return
        self.last_bucket = bucket
        width = 20
        filled = min(width, int(percent * width / 100))
        bar = "#" * filled + "-" * (width - filled)
        print(
            f"\rSummarizing [{bar}] {percent:3d}% ({row}/{total} tasks)",
            end="\n" if percent == 100 else "",
            file=sys.stderr,
            flush=True,
        )


class StorageProgress:
    def __init__(self) -> None:
        self.last_counts: dict[str, int] = {}

    def __call__(self, label: str, count: int) -> None:
        previous = self.last_counts.get(label)
        if previous == count:
            print(
                f"\rMeasuring {label}: {count} files\n",
                end="",
                file=sys.stderr,
                flush=True,
            )
            return
        if count == 1 or count % 1000 == 0:
            print(
                f"\rMeasuring {label}: {count} files",
                end="",
                file=sys.stderr,
                flush=True,
            )
        self.last_counts[label] = count


def _storage_argument(value: str) -> tuple[str, Path]:
    label, separator, raw_path = value.partition("=")
    if (
        not separator
        or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", label) is None
        or not raw_path
    ):
        raise argparse.ArgumentTypeError("storage must use LABEL=PATH")
    return label, Path(raw_path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="rnaseq-service-hpc")
    subparsers = parser.add_subparsers(dest="command", required=True)
    render = subparsers.add_parser("render-config", help="render a SLURM/Apptainer config")
    render.add_argument("--settings", type=Path, required=True)
    render.add_argument("--output", type=Path, required=True)
    prepare = subparsers.add_parser("prepare", help="render an immutable sbatch launcher")
    prepare.add_argument("--run-plan", type=Path, required=True)
    prepare.add_argument("--settings", type=Path, required=True)
    prepare.add_argument("--output", type=Path, required=True)
    submit = subparsers.add_parser("submit", help="submit one launcher and reserve a receipt")
    submit.add_argument("--launcher", type=Path, required=True)
    submit.add_argument("--receipt", type=Path, required=True)
    status = subparsers.add_parser("status", help="query a submitted controller job")
    status.add_argument("--receipt", type=Path, required=True)
    summary = subparsers.add_parser("summarize", help="summarize a Nextflow trace TSV")
    summary.add_argument("--trace", type=Path, required=True)
    summary.add_argument("--json-out", type=Path, required=True)
    summary.add_argument("--tsv-out", type=Path, required=True)
    summary.add_argument(
        "--storage",
        action="append",
        default=[],
        type=_storage_argument,
        metavar="LABEL=PATH",
        help="measure a work or results directory; repeat for multiple paths",
    )
    summary.add_argument("--quiet", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "render-config":
            result = write_nextflow_config(args.settings, args.output)
        elif args.command == "prepare":
            result = prepare_launcher(
                run_plan=args.run_plan,
                settings_path=args.settings,
                output=args.output,
            )
        elif args.command == "submit":
            result = submit_launcher(launcher=args.launcher, receipt=args.receipt)
        elif args.command == "status":
            result = query_job(args.receipt)
        else:
            storage_paths = dict(args.storage)
            if len(storage_paths) != len(args.storage):
                raise ResourceSummaryError("storage labels must be unique")
            result = write_resource_summary(
                trace=args.trace,
                json_output=args.json_out,
                tsv_output=args.tsv_out,
                progress=None if args.quiet else RowProgress(),
                storage_paths=storage_paths,
                storage_progress=None if args.quiet else StorageProgress(),
            )
    except (HPCConfigError, HPCRunError, ResourceSummaryError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
