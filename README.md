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

## Operating modes

- **Direct:** stage pinned workflows, containers, and references through the
  available network.
- **Proxy:** inherit runtime `HTTP_PROXY`, `HTTPS_PROXY`, and `NO_PROXY`
  settings without committing or logging credentials.
- **Offline:** verify and use a fully pre-staged, checksum-bound bundle on
  restricted compute nodes.

The staging and execution implementations arrive in later milestones; M0 locks
their security and reproducibility requirements first.

## Documentation

- [`docs/SERVICE_CONTRACT.md`](docs/SERVICE_CONTRACT.md) — supported intake,
  delivery gates, responsibilities, and exclusions;
- [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) — workflow boundary, network
  modes, security, and reproducibility;
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
