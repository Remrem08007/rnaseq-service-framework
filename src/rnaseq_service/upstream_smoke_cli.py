"""CLI for public-data upstream workflow smoke validation."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .upstream_smoke import (
    UpstreamSmokeError,
    create_upstream_smoke_launcher,
    create_upstream_smoke_plan,
    create_upstream_smoke_receipt,
)


class ArtifactProgress:
    def __init__(self) -> None:
        self.label: str | None = None
        self.bucket = -1

    def __call__(self, label: str, current: int, total: int) -> None:
        percent = 100 if total == 0 else min(100, int(current * 100 / total))
        bucket = percent // 5
        if label == self.label and bucket == self.bucket:
            return
        self.label = label
        self.bucket = bucket
        width = 20
        filled = min(width, int(percent * width / 100))
        bar = "#" * filled + "-" * (width - filled)
        print(
            f"\rHashing {label:<34} [{bar}] {percent:3d}%",
            end="\n" if percent == 100 else "",
            file=sys.stderr,
            flush=True,
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="rnaseq-service-upstream-smoke")
    commands = parser.add_subparsers(dest="command", required=True)
    plan = commands.add_parser("plan", help="plan pinned public-data upstream test profiles")
    plan.add_argument("--workflow-lock", type=Path, required=True)
    plan.add_argument("--output", type=Path, required=True)
    plan.add_argument("--run-root", type=Path, required=True)
    plan.add_argument(
        "--network-mode",
        choices=("auto", "direct", "proxy", "offline"),
        default="auto",
    )
    plan.add_argument("--offline-manifest", type=Path)
    plan.add_argument(
        "--container-engine",
        choices=("apptainer", "docker", "singularity"),
        default="apptainer",
    )
    prepare = commands.add_parser("prepare", help="render a heartbeat SLURM launcher")
    prepare.add_argument("--run-plan", type=Path, required=True)
    prepare.add_argument("--output", type=Path, required=True)
    prepare.add_argument("--account", required=True)
    prepare.add_argument("--partition", required=True)
    prepare.add_argument("--time", required=True)
    prepare.add_argument("--cpus", type=int, default=8)
    prepare.add_argument("--memory-gb", type=int, default=32)
    prepare.add_argument("--module", action="append", dest="modules", required=True)
    prepare.add_argument("--no-module-purge", action="store_false", dest="purge_modules")
    complete = commands.add_parser(
        "complete", help="validate both upstream runs and seal a completion receipt"
    )
    complete.add_argument("--run-plan", type=Path, required=True)
    complete.add_argument("--output", type=Path, required=True)
    complete.add_argument("--quiet", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "plan":
            result = create_upstream_smoke_plan(
                workflow_lock=args.workflow_lock,
                output=args.output,
                run_root=args.run_root,
                network_mode=args.network_mode,
                container_engine=args.container_engine,
                offline_manifest=args.offline_manifest,
            )
            summary = {
                "status": result["stage"],
                "stages": [stage["id"] for stage in result["stages"]],
                "requested_network_mode": result["execution"]["network_mode"],
                "selected_network_mode": result["execution"][
                    "selected_network_mode"
                ],
                "offline_bundle_verified": result["execution"][
                    "offline_bundle_verified"
                ],
                "uses_public_upstream_test_data": result["uses_public_upstream_test_data"],
                "scientific_execution_expected": result["scientific_execution_expected"],
            }
        elif args.command == "prepare":
            result = create_upstream_smoke_launcher(
                run_plan=args.run_plan,
                output=args.output,
                account=args.account,
                partition=args.partition,
                wall_time=args.time,
                cpus=args.cpus,
                memory_gb=args.memory_gb,
                modules=args.modules,
                purge_modules=args.purge_modules,
            )
            summary = result
        else:
            result = create_upstream_smoke_receipt(
                run_plan=args.run_plan,
                output=args.output,
                progress=None if args.quiet else ArtifactProgress(),
            )
            summary = {
                "status": result["stage"],
                "scientific_execution_performed": result["scientific_execution_performed"],
                "uses_public_upstream_test_data": result["uses_public_upstream_test_data"],
                "stages": [
                    {
                        "id": stage["id"],
                        "tasks": stage["task_summary"]["task_count"],
                        "artifacts": len(stage["artifacts"]),
                    }
                    for stage in result["stages"]
                ],
            }
    except (UpstreamSmokeError, OSError) as exc:
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
