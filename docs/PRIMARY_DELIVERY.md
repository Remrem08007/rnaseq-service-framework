# Primary RNA-seq delivery

M4 connects the validated intake, pinned software, staged dependencies, and HPC
launcher into a primary-processing delivery gate for the current
`nf-core/rnaseq` STAR/Salmon profile. A successful pipeline is recorded as
`pipeline_complete_qc_pending`; it is not treated as scientific QC acceptance.

## Supported reference contract

Every new CLI-generated plan must name a reviewed nf-core iGenomes key with
`--genome`. On HPC, `paths.reference_cache` in the private infrastructure
settings supplies the shared iGenomes root. The key is stored in the plan,
passed directly to nf-core, and checked against the completed pipeline's
`params.json`.

Custom FASTA/GTF and custom prebuilt-index manifests are not accepted by this
version of the service wrapper. Do not substitute custom paths behind an
iGenomes key; extending the contract requires explicit checksums and tests.

## 1. Seal the FASTQ dataset

After intake preflight and before run planning:

```bash
rnaseq-service-inputs \
  --samplesheet private/intake/samplesheet.csv \
  --output private/inputs/study-001.fastq-manifest.json
```

This command reads every FASTQ once and shows dataset-wide hashing progress. It
records each resolved file path, sample/read role, size, nanosecond modification
time, and SHA-256. The mode-`0600` manifest is never overwritten. The manifest
contains client file identities and must remain private.

Planning verifies the samplesheet checksum and current FASTQ metadata. SLURM
launcher preparation performs the metadata verification again, so a moved,
resized, or modified input is rejected before submission. Full hashes remain
in the immutable manifest without forcing another complete read on a shared
filesystem.

## 2. Create the run plan

```bash
rnaseq-service-plan \
  --samplesheet private/intake/samplesheet.csv \
  --input-manifest private/inputs/study-001.fastq-manifest.json \
  --genome GRCh38 \
  --design private/intake/design.csv \
  --contrasts private/intake/contrasts.csv \
  --output private/plans/study-001.json \
  --outdir results/study-001 \
  --workdir /shared/project/work/study-001 \
  --network-mode offline \
  --offline-manifest /shared/software/rnaseq/offline_bundle.manifest.json \
  --container-engine apptainer \
  --executor slurm \
  --infrastructure-config private/hpc/nextflow.config
```

Replace `GRCh38` with the scientifically reviewed key available in the
configured iGenomes resource. The plan binds the input manifest as a hashed
control file and adds `--genome <key>` to the canonical argv vector.

## 3. Prepare, submit, and monitor

```bash
rnaseq-service-hpc prepare \
  --run-plan private/plans/study-001.json \
  --settings private/hpc/infrastructure.toml \
  --output private/launch/study-001/attempt-1/controller.sbatch

rnaseq-service-hpc submit \
  --launcher private/launch/study-001/attempt-1/controller.sbatch \
  --receipt private/launch/study-001/attempt-1/submission.json

rnaseq-service-hpc status \
  --receipt private/launch/study-001/attempt-1/submission.json
```

The controller launches the plan's exact argv and always includes `-resume`.
Each attempt gets a distinct launcher directory and submission receipt.

## 4. Diagnose and safely prepare a retry

If the scheduler job fails or times out:

```bash
rnaseq-service-primary diagnose \
  --run-plan private/plans/study-001.json \
  --submission-receipt private/launch/study-001/attempt-1/submission.json
```

Diagnostics combine trace state with `.nextflow.log` and controller log pattern
matches. They report categories, paths, and line numbers—not raw log lines—so a
URL credential or other sensitive log text is not copied into JSON output.

After reviewing the cause and any authorized infrastructure change, prepare a
new launcher:

```bash
rnaseq-service-primary restart \
  --run-plan private/plans/study-001.json \
  --settings private/hpc/infrastructure.toml \
  --previous-receipt private/launch/study-001/attempt-1/submission.json \
  --output private/launch/study-001/attempt-2/controller.sbatch
```

Restart preparation refuses an active or successful previous job, a changed
launcher, missing work cache, changed plan control, or already-complete trace.
It creates but does not submit a fresh `-resume` launcher. Review it, then use
`rnaseq-service-hpc submit` with a new attempt-2 receipt.

## 5. Record compute observations

After a successful run, create the M3 resource summary:

```bash
rnaseq-service-hpc summarize \
  --trace results/study-001/execution/trace.tsv \
  --json-out results/study-001/execution/resource_summary.json \
  --tsv-out results/study-001/execution/resource_by_process.tsv \
  --storage work=/shared/project/work/study-001 \
  --storage results=results/study-001
```

## 6. Seal pipeline completion

```bash
rnaseq-service-primary complete \
  --run-plan private/plans/study-001.json \
  --output private/delivery/study-001.primary-completion.json
```

The command checks:

- every plan control file, including the FASTQ manifest;
- only successful terminal task states in the Nextflow trace;
- the outer Nextflow report, trace, timeline, and DAG;
- nf-core software versions, parameters, and validated samplesheet;
- the default STAR/Salmon gene-count, gene-length, and gene-TPM matrices;
- every planned biological sample in all three matrix headers;
- exactly one merged MultiQC report and exactly one non-empty parsed-data
  directory named `multiqc_report_data` (current 3.26 output) or
  `multiqc_data` (documented/legacy output);
- agreement of pipeline input, output, and genome parameters with the plan.

It then hashes each required artifact and every parsed MultiQC data file with a
visible progress bar and writes a mode-`0600`, non-overwriting schema-version-2
receipt. The receipt inventories paths, sizes, SHA-256 values, task counts,
sample/matrix checks, workflow/runtime versions, the reference key, and MultiQC
data storage.

The final state is intentionally `pipeline_complete_qc_pending`. Continue with
[`QC_ACCEPTANCE.md`](QC_ACCEPTANCE.md); pipeline success is not QC acceptance.
