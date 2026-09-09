"""CLI for staging and verifying official nf-core smoke-test inputs."""

from __future__ import annotations

import argparse
import json
import sys
import urllib.error
from pathlib import Path

from .upstream_data import (
    UpstreamDataError,
    stage_upstream_test_data,
    verify_upstream_test_data,
)


class DownloadProgress:
    def __init__(self) -> None:
        self.label: str | None = None
        self.bucket = -1

    def __call__(self, label: str, current: int, total: int) -> None:
        percent = None if total <= 0 else min(100, int(current * 100 / total))
        bucket = current // (8 * 1024 * 1024) if percent is None else percent // 5
        if label == self.label and bucket == self.bucket:
            return
        self.label = label
        self.bucket = bucket
        if percent is None:
            status = f"{current / (1024 * 1024):.1f} MiB"
        else:
            width = 20
            filled = min(width, int(percent * width / 100))
            status = f"[{'#' * filled}{'-' * (width - filled)}] {percent:3d}%"
        print(
            f"\rDownloading {label:<42} {status}",
            end="\n" if percent == 100 else "",
            file=sys.stderr,
            flush=True,
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="rnaseq-service-upstream-data")
    commands = parser.add_subparsers(dest="command", required=True)
    stage = commands.add_parser(
        "stage", help="download all pinned public smoke-test inputs"
    )
    stage.add_argument("--output-root", type=Path, required=True)
    stage.add_argument(
        "--network-mode", choices=("auto", "direct", "proxy"), default="auto"
    )
    stage.add_argument("--timeout-seconds", type=int, default=60)
    stage.add_argument("--quiet", action="store_true")
    verify = commands.add_parser("verify", help="verify every staged public test input")
    verify.add_argument("--root", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "stage":
            result = stage_upstream_test_data(
                output_root=args.output_root,
                network_mode=args.network_mode,
                timeout_seconds=args.timeout_seconds,
                progress=None if args.quiet else DownloadProgress(),
            )
        else:
            result = verify_upstream_test_data(args.root)
    except (OSError, UpstreamDataError, urllib.error.URLError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
