"""CLI for producing a reviewed, non-executing primary RNA-seq run plan."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .plan import PlanError, PreflightPlanError, create_rnaseq_plan


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rnaseq-service-plan",
        description="Validate intake and write an immutable nf-core/rnaseq plan.",
    )
    parser.add_argument("--samplesheet", type=Path, required=True)
    parser.add_argument("--input-manifest", type=Path, required=True)
    parser.add_argument(
        "--genome",
        required=True,
        help="reviewed nf-core iGenomes key, for example GRCh38",
    )
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument("--contrasts", type=Path, required=True)
    parser.add_argument(
        "--workflow-lock",
        type=Path,
        default=Path("config/workflows.toml"),
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument(
        "--network-mode",
        choices=("auto", "direct", "proxy", "offline"),
        required=True,
    )
    parser.add_argument(
        "--container-engine",
        choices=("apptainer", "docker", "singularity"),
        default="apptainer",
    )
    parser.add_argument(
        "--executor",
        choices=("local", "slurm"),
        default="local",
    )
    parser.add_argument("--infrastructure-config", type=Path)
    parser.add_argument(
        "--offline-manifest",
        type=Path,
        help="sealed bundle manifest; required for offline and optional auto fallback",
    )
    parser.add_argument("--min-replicates", type=int, default=2)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.min_replicates < 1:
        raise SystemExit("--min-replicates must be at least 1")
    try:
        plan = create_rnaseq_plan(
            samplesheet=args.samplesheet,
            design=args.design,
            contrasts=args.contrasts,
            workflow_lock=args.workflow_lock,
            output=args.output,
            outdir=args.outdir,
            workdir=args.workdir,
            network_mode=args.network_mode,
            container_engine=args.container_engine,
            executor=args.executor,
            infrastructure_config=args.infrastructure_config,
            offline_manifest=args.offline_manifest,
            min_replicates=args.min_replicates,
            input_manifest=args.input_manifest,
            genome=args.genome,
        )
    except PreflightPlanError as exc:
        print(exc.report.to_json(), file=sys.stderr, end="")
        return 2
    except (PlanError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(plan, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
