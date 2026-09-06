"""CLI for sealing completed nf-core/rnaseq primary-processing runs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .primary import PrimaryDeliveryError, create_primary_receipt


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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rnaseq-service-primary",
        description="Validate primary outputs and write an immutable completion receipt.",
    )
    parser.add_argument("--run-plan", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--quiet", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = create_primary_receipt(
            run_plan=args.run_plan,
            output=args.output,
            progress=None if args.quiet else ArtifactProgress(),
        )
    except (PrimaryDeliveryError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
