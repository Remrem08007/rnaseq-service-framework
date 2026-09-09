# Immutable run plan

M1 separates planning from execution. `rnaseq-service-plan` performs the full
intake preflight, reads the reviewed workflow lock, hashes every control file,
builds a Nextflow argument vector, and writes an exclusive JSON plan. It does
not invoke Nextflow, contact a registry, submit a job, or create analysis
results.

## Create a plan

Seal the FASTQ dataset before planning:

```bash
rnaseq-service-inputs \
  --samplesheet private/intake/samplesheet.csv \
  --output private/inputs/study-001.fastq-manifest.json
```

Real FASTQ files must exist because manifest creation and planning enable
filesystem checks.

```bash
rnaseq-service-plan \
  --samplesheet private/intake/samplesheet.csv \
  --input-manifest private/inputs/study-001.fastq-manifest.json \
  --genome GRCh38 \
  --design private/intake/design.csv \
  --contrasts private/intake/contrasts.csv \
  --workflow-lock config/workflows.toml \
  --output private/plans/study-001.json \
  --outdir results/study-001 \
  --workdir work/study-001 \
  --network-mode direct \
  --container-engine apptainer \
  --executor local
```

For SLURM planning, `--executor slurm` requires an explicit
`--infrastructure-config`. M1 records and hashes that file; the generic SLURM
templates themselves arrive in M3.

## Plan contents

- exact upstream workflow name and revision;
- creation time and schema version;
- SHA-256, size, and resolved path for each control file;
- a bound FASTQ manifest and reviewed iGenomes reference key;
- complete preflight outcome;
- executor, container engine, network mode, output directory, and work directory;
- a structured `command_argv` list;
- a shell-escaped preview for human review;
- proxy-variable presence flags, never proxy values;
- explicit `planned_not_executed` status.

M4's input-manifest step reads every FASTQ once with visible aggregate progress
and records its SHA-256, size, modification time, sample, and read role.
Planning binds the manifest and checks current metadata without rereading the
entire dataset. SLURM launcher preparation repeats that quick metadata check.

## Immutability

The output is created with mode `0600` using exclusive creation. An existing
plan is never overwritten. A changed design, workflow revision, or execution
setting requires a new plan path so the earlier decision record remains intact.

## Command safety

The canonical command is stored as an argument vector. The implementation does
not concatenate user-controlled paths into a command passed to `shell=True`.
The rendered preview is informational and shell-quoted with `shlex.join`.

## Proxy and offline status

Only standard proxy environment-variable names are recorded; their values are
never serialized. M2 requires `--offline-manifest` whenever
`--network-mode offline` is selected. `--network-mode auto` optionally accepts
the same manifest and binds an offline fallback alongside the pinned online
command. The SLURM launcher prefers a verified offline bundle when supplied;
without one it resolves direct access and then an HTTPS proxy before Nextflow
starts. Planning verifies every bundled artifact;
offline execution switches from the remote workflow name to the verified local
workflow directory and records `NXF_OFFLINE=true`. See
[`OFFLINE_STAGING.md`](OFFLINE_STAGING.md) for the complete staging and
verification sequence.
