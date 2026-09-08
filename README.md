# RNA-seq Service Framework

[![CI](https://github.com/Remrem08007/rnaseq-service-framework/actions/workflows/ci.yml/badge.svg)](https://github.com/Remrem08007/rnaseq-service-framework/actions/workflows/ci.yml)

A client-delivery framework for reproducible bulk RNA-seq analysis using pinned
community workflows, validated study designs, and portable HPC execution.

This repository does **not** reimplement RNA-seq alignment or differential
expression. It adds the operational layer needed to deliver those analyses
reliably: intake validation, design checks, version locking, infrastructure
profiles, offline/proxy-aware staging, provenance, failure diagnostics, and
client-facing handoff reports.

## Upstream workflows

The current lock records:

- [`nf-core/rnaseq` 3.26.0](https://nf-co.re/rnaseq/3.26.0/) for read QC,
  trimming, alignment or pseudoalignment, quantification, and expression
  matrices;
- [`nf-core/differentialabundance` 2.0.0](https://nf-co.re/differentialabundance/2.0.0/)
  for statistical contrasts, differential analysis, enrichment, and reporting.

Versions are pinned in [`config/workflows.toml`](config/workflows.toml) and are
changed only through a reviewed, retested update.

The same lock pins Nextflow 26.04.6 and nf-core/tools 4.1.0 for deterministic
staging behavior.

## M0 preflight

The first implementation validates the client intake contract before any
workflow or SLURM job can be launched:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'

rnaseq-service-preflight \
  --samplesheet examples/intake/samplesheet.csv \
  --design examples/intake/design.csv \
  --contrasts examples/intake/contrasts.csv
```

The command prints a JSON report and exits with status `2` when validation
fails. Add `--check-files` for a real engagement to require every FASTQ path to
exist. Use `--json-out preflight.json` to preserve the validation artifact.

The sample sheet follows the four-column contract documented by
`nf-core/rnaseq`: `sample`, `fastq_1`, `fastq_2`, and `strandedness`. The design
contains one row per biological sample, and contrast direction is explicit as
reference versus target.

## M1 run planning

Once real FASTQs and reviewed intake tables exist, write a non-executing,
immutable run plan:

```bash
rnaseq-service-inputs \
  --samplesheet private/intake/samplesheet.csv \
  --output private/inputs/study-001.fastq-manifest.json

rnaseq-service-plan \
  --samplesheet private/intake/samplesheet.csv \
  --input-manifest private/inputs/study-001.fastq-manifest.json \
  --genome GRCh38 \
  --design private/intake/design.csv \
  --contrasts private/intake/contrasts.csv \
  --output private/plans/study-001.json \
  --outdir results/study-001 \
  --workdir work/study-001 \
  --network-mode direct
```

The input command hashes every FASTQ with visible progress. Planning requires
that immutable manifest, validates the intake and reviewed iGenomes key, reads
the exact workflow lock, hashes every control file, and stores the canonical
Nextflow argument vector. It does not run Nextflow or contact the network.
Existing manifests and plans are never overwritten. See
[`docs/RUN_PLAN.md`](docs/RUN_PLAN.md).

## Operating modes

- **Direct:** stage pinned workflows, containers, and references through the
  available network.
- **Proxy:** inherit runtime `HTTP_PROXY`, `HTTPS_PROXY`, and `NO_PROXY`
  settings without committing or logging credentials.
- **Offline:** verify and use a fully pre-staged, checksum-bound bundle on
  restricted compute nodes.

M2 implements staging and offline verification. It downloads both workflows
and their Apptainer images with visible progress, supports standard proxy
environment variables without storing their values, seals every artifact with
SHA-256, and refuses offline planning until the bundle verifies. See
[`docs/OFFLINE_STAGING.md`](docs/OFFLINE_STAGING.md).

## M3 HPC execution

M3 validates private infrastructure settings and renders a generic
SLURM/Apptainer Nextflow configuration. It prepares an immutable resumable
controller script, reserves a receipt before submission, reports scheduler
state through `squeue`/`sacct`, and summarizes task-level Nextflow trace metrics
plus measured work/results storage with visible progress. Cluster accounts,
partitions, modules, and filesystem paths remain external. See
[`docs/HPC_EXECUTION.md`](docs/HPC_EXECUTION.md).

## M4 primary delivery

M4 verifies completed `nf-core/rnaseq` STAR/Salmon outputs against the locked
plan, successful Nextflow task states, pipeline parameters, sample-complete
count/TPM matrices, MultiQC, and software provenance. Required artifacts are
hashed with progress into an immutable `pipeline_complete_qc_pending` receipt.
Failure diagnostics classify logs without copying their content, and safe
restart preparation refuses active/successful jobs before creating a reviewed
`-resume` launcher. See [`docs/PRIMARY_DELIVERY.md`](docs/PRIMARY_DELIVERY.md).

## M5 QC acceptance

M5 reads checksum-bound parsed MultiQC evidence through a study-reviewed TOML
policy. Numeric mapping, duplication, and complexity measures can have explicit
limits and robust median/MAD outlier flags; categorical evidence supports checks
such as strandedness agreement. Missing required evidence produces a visible
stop state. A flag never removes a sample automatically: a reviewer must record
include/exclude decisions, reviewers, and reasons before the framework seals an
accepted-sample manifest and verifies every requested contrast still has enough
replicates. See [`docs/QC_ACCEPTANCE.md`](docs/QC_ACCEPTANCE.md).

## Documentation

- [`docs/SERVICE_CONTRACT.md`](docs/SERVICE_CONTRACT.md) — supported intake,
  delivery gates, responsibilities, and exclusions;
- [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) — workflow boundary, network
  modes, security, and reproducibility;
- [`docs/RUN_PLAN.md`](docs/RUN_PLAN.md) — immutable planning, hashes, command
  safety, and execution boundary;
- [`docs/OFFLINE_STAGING.md`](docs/OFFLINE_STAGING.md) — pinned downloads,
  proxy handling, progress, bundle sealing, and offline verification;
- [`docs/HPC_EXECUTION.md`](docs/HPC_EXECUTION.md) — external infrastructure
  settings, resumable SLURM launches, status, and resource observations;
- [`docs/PRIMARY_DELIVERY.md`](docs/PRIMARY_DELIVERY.md) — FASTQ provenance,
  primary completion evidence, diagnostics, and safe retries;
- [`docs/QC_ACCEPTANCE.md`](docs/QC_ACCEPTANCE.md) — configurable MultiQC
  evidence, flags, human decisions, and accepted samples;
- [`ROADMAP.md`](ROADMAP.md) — M0–M8 implementation plan.

## Data and privacy policy

Client FASTQs, alignments, count matrices, metadata, reports, work directories,
credentials, and proxy secrets are excluded from Git. Public tests use synthetic
metadata and placeholder paths only.

## Scope limitation

This framework supports research bioinformatics delivery. It does not provide
clinical diagnosis, replace experimental-design review, infer causal biology,
or automatically choose a scientifically appropriate comparison.

## License

Repository code and documentation are released under the [MIT License](LICENSE).
Upstream workflows, containers, reference resources, and client data retain
their own terms.
