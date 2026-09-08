"""CLI for accepted-sample differential-analysis handoff planning."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .differential import DifferentialHandoffError, create_differential_plan


class MatrixProgressBar:
    def __init__(self) -> None:
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
            f"\rSubsetting {label:<12} [{bar}] {percent:3d}%",
            end="\n" if percent == 100 else "",
            file=sys.stderr,
            flush=True,
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rnaseq-service-differential",
        description="Create a pinned differential-analysis plan from sealed QC acceptance.",
    )
    parser.add_argument("--qc-acceptance", type=Path, required=True)
    parser.add_argument("--handoff-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--study-name", required=True)
    parser.add_argument("--network-mode", choices=("direct", "proxy", "offline"), required=True)
    parser.add_argument(
        "--container-engine", choices=("apptainer", "docker", "singularity"), default="apptainer"
    )
    parser.add_argument("--executor", choices=("local", "slurm"), default="local")
    parser.add_argument("--infrastructure-config", type=Path)
    parser.add_argument("--offline-manifest", type=Path)
    parser.add_argument("--seed", type=int, default=20260908)
    parser.add_argument("--quiet", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        plan = create_differential_plan(
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
            progress=None if args.quiet else MatrixProgressBar(),
        )
    except (DifferentialHandoffError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "status": plan["stage"],
                "accepted_sample_count": plan["accepted_sample_count"],
                "feature_count": plan["feature_count"],
                "contrast_count": plan["contrast_count"],
                "workflow": plan["workflow"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
