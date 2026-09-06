# RNA-seq Service Framework

A client-delivery framework for reproducible bulk RNA-seq analysis using pinned
community workflows, validated study designs, and portable HPC execution.

This repository does **not** reimplement RNA-seq alignment or differential
expression. It adds the operational layer needed to deliver those analyses
reliably: intake validation, design checks, version locking, infrastructure
profiles, offline/proxy-aware staging, provenance, failure diagnostics, and
client-facing handoff reports.

## Planned upstream workflows

- [nf-core/rnaseq](https://nf-co.re/rnaseq/) for read QC, trimming,
  alignment/pseudoalignment, quantification, and expression matrices.
- [nf-core/differentialabundance](https://nf-co.re/differentialabundance/) for
  statistical contrasts, differential analysis, enrichment, and reporting.

Exact versions are pinned in repository configuration and changed only through a
reviewed update.

## Operating modes

- direct network access;
- proxy-aware staging using environment-provided settings;
- fully pre-staged offline execution on restricted HPC compute nodes.

Proxy URLs and credentials are never committed to Git or copied into reports.

## Status

M0—service contract and preflight validation—is in progress. See
[ROADMAP.md](ROADMAP.md) once the first milestone is merged.

## Scope and limitations

This framework supports research bioinformatics delivery. It does not make
clinical diagnostic claims, replace experimental-design review, or automatically
decide which biological comparison is scientifically appropriate.
