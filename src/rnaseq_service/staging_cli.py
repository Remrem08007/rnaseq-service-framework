"""CLI for proxy-aware staging of pinned workflows and containers."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .staging import StagingError, create_staging_plan, run_staging_plan


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="rnaseq-service-stage")
    subparsers = parser.add_subparsers(dest="command", required=True)
    plan = subparsers.add_parser("plan", help="write a non-executing staging plan")
    plan.add_argument("--workflow-lock", type=Path, default=Path("config/workflows.toml"))
    plan.add_argument("--bundle-dir", type=Path, required=True)
    plan.add_argument("--output", type=Path, required=True)
    plan.add_argument("--network-mode", choices=("direct", "proxy"), required=True)
    plan.add_argument("--parallel-downloads", type=int, default=4)
    run = subparsers.add_parser("run", help="execute one immutable staging plan")
    run.add_argument("--plan", type=Path, required=True)
    run.add_argument("--heartbeat-seconds", type=int, default=30)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "plan":
            result = create_staging_plan(
                workflow_lock=args.workflow_lock,
                bundle_dir=args.bundle_dir,
                output=args.output,
                network_mode=args.network_mode,
                parallel_downloads=args.parallel_downloads,
            )
        else:
            result = run_staging_plan(
                args.plan,
                heartbeat_seconds=args.heartbeat_seconds,
            )
    except (StagingError, OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
