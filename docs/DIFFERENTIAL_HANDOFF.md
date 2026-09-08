# Differential-analysis handoff

M6 converts the human-reviewed M5 sample acceptance into a pinned
`nf-core/differentialabundance` 2.0.0 DESeq2 run. It never reopens excluded
samples or silently changes a contrast.

## What enters the analysis

The planner verifies the entire chain back to the primary run, then creates
private, mode-`0600` files containing only accepted samples:

- `observations.accepted.tsv` — sample metadata in the reviewed order;
- `gene_counts.accepted.tsv` — raw Salmon gene counts;
- `gene_lengths.accepted.tsv` — matching average transcript lengths;
- `contrasts.accepted.tsv` — explicit target-versus-reference comparisons.

Gene counts are the measurements tested by DESeq2. Gene lengths let DESeq2
apply the length-aware normalization used for transcript-level quantification.
The two matrices must have identical gene IDs and row order. Values must be
finite and non-negative. Every accepted sample must occur exactly once in the
observations and both matrices.

The baseline service profile locks `differential_method=deseq2`, variance
stabilization, HTML reporting, and `functional_method=none`. Pathway enrichment
is not inferred from a contrast or silently enabled: organism, identifier,
database, and gene-set choices require a separate scientific review.

## 1. Create the differential plan

For local execution:

```bash
rnaseq-service-differential plan \
  --qc-acceptance private/qc/study-001/qc_acceptance.json \
  --handoff-dir private/differential/study-001/inputs \
  --output private/differential/study-001/run-plan.json \
  --outdir results/study-001-differential \
  --workdir /shared/project/work/study-001-differential \
  --study-name study001 \
  --network-mode direct
```

The progress bars stream the count and length matrices once. The resulting
plan records SHA-256 checksums, exact workflow/runtime revisions, accepted
samples, feature and contrast counts, reference genome, execution mode, and a
shell-safe argument vector. Planning does not run Nextflow. Existing handoff
directories and plans are never overwritten.

For SLURM, first render the infrastructure configuration as described in
[`HPC_EXECUTION.md`](HPC_EXECUTION.md), then add:

```bash
  --executor slurm \
  --infrastructure-config private/hpc/nextflow.config
```

Use `--network-mode proxy` when approved proxy variables are already present in
the environment. For a verified offline bundle, use `--network-mode offline`
and add `--offline-manifest`; proxy values and credentials are never recorded.

## 2. Launch on SLURM

The generic M3 launcher supports both the primary and differential plans:

```bash
rnaseq-service-hpc prepare \
  --run-plan private/differential/study-001/run-plan.json \
  --settings private/hpc/infrastructure.toml \
  --output private/differential/study-001/launch/controller.sbatch

rnaseq-service-hpc submit \
  --launcher private/differential/study-001/launch/controller.sbatch \
  --receipt private/differential/study-001/launch/submission.json

rnaseq-service-hpc status \
  --receipt private/differential/study-001/launch/submission.json
```

Preparation re-verifies the acceptance receipt and every derived input. Submit
reserves a receipt before calling `sbatch`, preventing an accidental duplicate
submission. Nextflow uses `-resume`; a retry must follow the reviewed restart
process rather than deleting the work directory.

## 3. Seal successful outputs

After the scheduler reports completion, run:

```bash
rnaseq-service-differential complete \
  --run-plan private/differential/study-001/run-plan.json \
  --output private/differential/study-001/completion.json
```

This command requires a successful terminal Nextflow trace and verifies:

- execution report, trace, timeline, and DAG;
- pinned pipeline parameters, accepted observations, and software versions;
- normalized-count and variance-stabilized matrices with the accepted samples;
- full, filtered, and annotated DESeq2 tables for every planned contrast;
- a volcano plot for every contrast plus dispersion, PCA, and MAD plots;
- the HTML analysis report and editable report bundle.

The workflow may remove genes under its locked minimum-abundance filter, so the
analyzed feature count can be lower than the input feature count. It may never
be higher; normalized/VST row counts must agree; each full contrast table must
have the analyzed row count; and a filtered table cannot contain more rows than
its full table.

Every required artifact is hashed with visible progress. The mode-`0600`
completion receipt is non-overwriting and ends in
`differential_complete_interpretation_pending`.

## Interpreting a contrast

For a contrast whose target is `treated` and reference is `control`, a positive
DESeq2 `log2FoldChange` means higher modeled abundance in treated samples; a
negative value means lower abundance. `padj` is the multiple-testing-adjusted
p-value used by the locked filter. Effect size, uncertainty, count support,
sample QC, design assumptions, batch variables, and biological plausibility all
remain part of human interpretation.

Pipeline completion is not proof of causality, mechanism, clinical utility, or
biological importance. The framework records no automatic biological claims.

