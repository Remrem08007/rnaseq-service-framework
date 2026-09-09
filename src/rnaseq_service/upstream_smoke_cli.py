"""CLI for public-data upstream workflow smoke validation."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .upstream_smoke import UpstreamSmokeError, create_upstream_smoke_plan


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="rnaseq-service-upstream-smoke")
    commands = parser.add_subparsers(dest="command", required=True)
    plan = commands.add_parser("plan", help="plan pinned public-data upstream test profiles")
    plan.add_argument("--workflow-lock", type=Path, required=True)
    plan.add_argument("--output", type=Path, required=True)
    plan.add_argument("--run-root", type=Path, required=True)
    plan.add_argument("--network-mode", choices=("direct", "proxy"), required=True)
    plan.add_argument(
        "--container-engine",
        choices=("apptainer", "docker", "singularity"),
        default="apptainer",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = create_upstream_smoke_plan(
            workflow_lock=args.workflow_lock,
            output=args.output,
            run_root=args.run_root,
            network_mode=args.network_mode,
            container_engine=args.container_engine,
        )
    except (UpstreamSmokeError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "status": result["stage"],
                "stages": [stage["id"] for stage in result["stages"]],
                "uses_public_upstream_test_data": result["uses_public_upstream_test_data"],
                "scientific_execution_expected": result["scientific_execution_expected"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

