"""CLI for sealing and verifying pre-staged offline bundles."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .bundle import BundleError, seal_bundle, verify_bundle


class ProgressPrinter:
    """Small dependency-free byte progress indicator for large images."""

    def __init__(self, label: str) -> None:
        self.label = label
        self.last_bucket = -1

    def __call__(self, done: int, total: int, current: str) -> None:
        percent = 100 if total == 0 else min(100, int(done * 100 / total))
        bucket = percent // 5
        if bucket == self.last_bucket:
            return
        self.last_bucket = bucket
        width = 20
        filled = min(width, int(percent * width / 100))
        bar = "#" * filled + "-" * (width - filled)
        print(
            f"\r{self.label} [{bar}] {percent:3d}% {current}",
            end="\n" if percent == 100 else "",
            file=sys.stderr,
            flush=True,
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="rnaseq-service-bundle")
    subparsers = parser.add_subparsers(dest="command", required=True)

    seal = subparsers.add_parser("seal", help="hash artifacts and create a manifest")
    seal.add_argument("--bundle-dir", type=Path, required=True)
    seal.add_argument("--workflow-lock", type=Path, default=Path("config/workflows.toml"))
    seal.add_argument("--rnaseq-workflow", type=Path, required=True)
    seal.add_argument("--differential-workflow", type=Path, required=True)
    seal.add_argument("--container-root", type=Path, required=True)
    seal.add_argument("--plugin-root", type=Path, required=True)
    seal.add_argument("--quiet", action="store_true")

    verify = subparsers.add_parser("verify", help="verify every sealed artifact")
    verify.add_argument("--manifest", type=Path, required=True)
    verify.add_argument("--workflow-lock", type=Path, default=Path("config/workflows.toml"))
    verify.add_argument("--quiet", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "seal":
            result = seal_bundle(
                bundle_dir=args.bundle_dir,
                workflow_lock=args.workflow_lock,
                rnaseq_workflow=args.rnaseq_workflow,
                differential_workflow=args.differential_workflow,
                container_root=args.container_root,
                plugin_root=args.plugin_root,
                progress=None if args.quiet else ProgressPrinter("Sealing  "),
            )
        else:
            result = verify_bundle(
                args.manifest,
                expected_workflow_lock=args.workflow_lock,
                progress=None if args.quiet else ProgressPrinter("Verifying"),
            )
    except (BundleError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
