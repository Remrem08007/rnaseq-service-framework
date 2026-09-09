# Synthetic end-to-end validation

M7 provides a fast, deterministic check of the service framework when no real
client run exists. It can run on a laptop, login node, CI runner, or allocated
compute node because it does not invoke Nextflow, containers, STAR, Salmon, or
DESeq2.

## Run the complete contract check

From the repository virtual environment:

```bash
python -m pip install -e '.[dev]'

rnaseq-service-synthetic validate-contract \
  --output-dir results/synthetic-contract \
  --workflow-lock config/workflows.toml
```

The command has a nine-stage progress bar. A successful summary looks like:

```json
{
  "direction_validation_passed": true,
  "n_genes": 4,
  "n_samples": 6,
  "scientific_execution_performed": false,
  "status": "synthetic_contract_validation_complete"
}
```

The output directory is deliberately non-overwriting. Delete or archive an old
synthetic result yourself, or choose a new output directory for another run.

## What it tests

The validator creates three control and three treated paired-end samples and
passes them through:

1. file-aware intake preflight;
2. SHA-256 FASTQ manifest creation;
3. pinned primary run planning;
4. primary output-contract and completion-receipt validation;
5. QC assessment and explicit all-sample human-style acceptance;
6. accepted-only differential input generation;
7. differential output-contract and completion-receipt validation;
8. expected fold-change-direction checks;
9. direct, offline, SLURM, and restart/resume control checks.

The final `synthetic_contract_receipt.json` checksum-binds each intermediate
plan and receipt. It also confirms that no sample was automatically excluded.

The direct smoke requires the remote workflow name and exact pinned revision.
The offline smoke builds a deliberately non-runnable placeholder bundle, seals
all six members, verifies exact membership and checksums, plans with a local
workflow path, removes the remote `-r` lookup, and requires
`NXF_OFFLINE=true`. The restart smoke prepares a second launcher using the same
work directory and `-resume`, but never calls `sbatch`.

## Known synthetic truth

The samples have four artificial genes. Both conditions have equal
housekeeping and neutral counts. The treated condition has more
`gene_treated_marker`; the control condition has more `gene_control_marker`.
The contract check therefore requires:

- positive target-versus-reference log2 fold change for the treated marker;
- negative target-versus-reference log2 fold change for the control marker;
- near-zero fold changes for the housekeeping and neutral genes.

To create only the paired FASTQs, mini FASTA/GTF, intake tables, and truth
manifest:

```bash
rnaseq-service-synthetic create \
  --output-dir results/synthetic-study
```

The gzip streams use a fixed timestamp, paths in the truth manifest are
relative, and the same seed produces byte-identical files in different
directories.

## Validation boundary

The framework creates upstream-shaped synthetic count, MultiQC, and DESeq2
artifacts so every wrapper gate can be tested quickly. These artifacts are not
outputs of the scientific tools. Both the CLI summary and final receipt record
`scientific_execution_performed: false` and
`synthetic_upstream_artifacts: true`.

This check proves orchestration, checksum, non-overwrite, sample-boundary,
direction, offline, and resume behavior. A real pinned nf-core test-profile run
is still required to validate container execution and upstream-tool behavior on
target infrastructure. It remains the final M7 item and must not use client
data.

