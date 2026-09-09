"""CLI for accepted-sample differential-analysis handoff planning."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .differential import (
    DifferentialHandoffError,
    create_differential_plan,
    create_differential_receipt,
)


class ProgressBar:
    def __init__(self, verb: str) -> None:
        self.verb = verb
        self.last: tuple[str, int] | None = None

    def __call__(self, label: str, current: int, total: int) -> None:
        percent = 100 if total == 0 else min(100, int(current * 100 / total))
        bucket = percent // 5
        state = (label, bucket)
        if state == self.last:
            return
        self.last = state
        width = 20
        filled = int(percent * width / 100)
        bar = "#" * filled + "-" * (width - filled)
        print(
            f"\r{self.verb} {label:<22} [{bar}] {percent:3d}%",
            end="\n" if percent == 100 else "",
            file=sys.stderr,
            flush=True,
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rnaseq-service-differential",
        description="Plan and seal pinned differential-abundance analysis.",
    )
    commands = parser.add_subparsers(dest="command", required=True)
    plan = commands.add_parser("plan", help="create an accepted-sample handoff and run plan")
    plan.add_argument("--qc-acceptance", type=Path, required=True)
    plan.add_argument("--handoff-dir", type=Path, required=True)
    plan.add_argument("--output", type=Path, required=True)
    plan.add_argument("--outdir", type=Path, required=True)
    plan.add_argument("--workdir", type=Path, required=True)
    plan.add_argument("--study-name", required=True)
    plan.add_argument(
        "--network-mode", choices=("auto", "direct", "proxy", "offline"), required=True
    )
    plan.add_argument(
        "--container-engine", choices=("apptainer", "docker", "singularity"), default="apptainer"
    )
    plan.add_argument("--executor", choices=("local", "slurm"), default="local")
    plan.add_argument("--infrastructure-config", type=Path)
    plan.add_argument("--offline-manifest", type=Path)
    plan.add_argument("--seed", type=int, default=20260908)
    plan.add_argument("--quiet", action="store_true")
    complete = commands.add_parser("complete", help="validate outputs and seal a completion receipt")
    complete.add_argument("--run-plan", type=Path, required=True)
    complete.add_argument("--output", type=Path, required=True)
    complete.add_argument("--quiet", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "plan":
            result = create_differential_plan(
                qc_acceptance=args.qc_acceptance,
                handoff_dir=args.handoff_dir,
                output=args.output,
                outdir=args.outdir,
                workdir=args.workdir,
                study_name=args.study_name,
                network_mode=args.network_mode,
                container_engine=args.container_engine,
                executor=args.executor,
                seed=args.seed,
                infrastructure_config=args.infrastructure_config,
                offline_manifest=args.offline_manifest,
                progress=None if args.quiet else ProgressBar("Subsetting"),
            )
        else:
            result = create_differential_receipt(
                run_plan=args.run_plan,
                output=args.output,
                progress=None if args.quiet else ProgressBar("Hashing"),
            )
    except (DifferentialHandoffError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "status": result["stage"],
                "accepted_sample_count": result["accepted_sample_count"],
                "feature_count": result["feature_count"],
                "contrast_count": result["contrast_count"],
                "workflow": result["workflow"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
