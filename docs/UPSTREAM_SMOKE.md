# Real upstream smoke validation

M7's final infrastructure check runs the official public-data test profiles for
the exact workflow revisions in [`config/workflows.toml`](../config/workflows.toml):

- `nf-core/rnaseq` 3.26.0 with profile `test`;
- `nf-core/differentialabundance` 2.0.0 with profile
  `test_rnaseq_deseq2_gsea`.

This is a real Nextflow, container, STAR/Salmon, DESeq2, and GSEA execution. It
uses upstream public test data, not client data. A pass demonstrates that the
pinned workflows and container runtime execute on the target infrastructure; it
does not establish accuracy for a future client's study or replace study-specific
QC and design review.

## 1. Install the current framework

Run on a login node after cloning or updating the repository:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[dev]'
python -m pytest
```

An existing virtual environment can be reactivated and reinstalled with the last
three commands. Hosted CI must also be green before the smoke run is interpreted.

## 2. Discover cluster settings

Alliance and NiBi are both suitable. Determine the values on the cluster where
the job will run instead of committing institution-specific settings:

```bash
module spider nextflow
module spider apptainer
sacctmgr show assoc where user="$USER" format=Account,Partition,QOS
sinfo -o '%P %a %l %c %m'
```

Set the exact values you reviewed. These examples are placeholders:

```bash
export SMOKE_ACCOUNT='replace-with-project-account'
export SMOKE_PARTITION='replace-with-partition'
export NEXTFLOW_MODULE='replace-with-nextflow-module/version'
export APPTAINER_MODULE='replace-with-apptainer-module/version'
```

The selected Nextflow module must provide the version pinned in the workflow
lock. Confirm both executables before planning:

```bash
module purge
module load "$NEXTFLOW_MODULE"
module load "$APPTAINER_MODULE"
nextflow -version
apptainer --version
```

If the cluster module is older than the lock, do not silently substitute it.
Install the exact official Nextflow release in a private executable directory,
then load only the container module in the generated launcher:

```bash
export REPO="$PWD"
mkdir -p "$HOME/.local/bin" "$HOME/.nextflow"
cd "$HOME/.local/bin"
export NXF_VER="$(python -c \
  'import sys,tomllib; print(tomllib.load(open(sys.argv[1], "rb"))["runtime"]["nextflow_version"])' \
  "$REPO/config/workflows.toml")"
curl -fsSL https://get.nextflow.io | bash
chmod 700 nextflow
export PATH="$HOME/.local/bin:$PATH"
nextflow -version
```

The launcher checks the observed Nextflow version against the lock before either
workflow starts. It also records the observed Nextflow and container-runtime
versions in `runtime/versions.tsv`; the completion receipt verifies and hashes
that file. When using a private Nextflow executable, omit the Nextflow module
from `prepare` and keep the Apptainer module.

If the site requires an outbound proxy, configure it in the submission shell.
The framework records only which proxy environment-variable names are present;
it never records their values. Automatic mode tries direct HTTPS and then the
configured proxy when no offline bundle is supplied.

For Narval or another cluster whose compute nodes cannot reach the internet,
first create and transfer the smoke-capable bundle described in
[`OFFLINE_STAGING.md`](OFFLINE_STAGING.md). The bundle must include the optional
`upstream-test-data` component; a software-only bundle is deliberately rejected
for an offline smoke run.

## 3. Create a unique non-executing plan

The run root must not already exist. This prevents stale results from being
mistaken for evidence from the planned run.

```bash
export SMOKE_ID="$(date -u +%Y%m%dT%H%M%SZ)"
export SMOKE_PRIVATE="$PWD/private/upstream-smoke/$SMOKE_ID"
export SMOKE_RUN="$PWD/results/upstream-smoke/$SMOKE_ID"

mkdir -p "$SMOKE_PRIVATE"
chmod 700 "$SMOKE_PRIVATE"

rnaseq-service-upstream-smoke plan \
  --workflow-lock config/workflows.toml \
  --output "$SMOKE_PRIVATE/run-plan.json" \
  --run-root "$SMOKE_RUN" \
  --network-mode auto \
  --offline-manifest /shared/path/offline-bundle/offline_bundle.manifest.json \
  --container-engine apptainer
```

Planning does not run Nextflow, download data, or submit a job. It checksum-binds
the workflow lock and records shell-safe argument arrays for both stages. When a
verified bundle is supplied, `auto` immediately selects it: the commands use the
local workflows and local public inputs, omit remote `-r` lookups, set
`NXF_OFFLINE=true`, and do not call `curl`. Bundle integrity is checked again
when the launcher is prepared.

If no bundle is available on a connected cluster, omit `--offline-manifest`.
The generated launcher probes direct GitHub/Quay access, then a configured proxy,
and exits before Nextflow if neither works. Explicit `direct`, `proxy`, and
`offline` modes remain available; `offline` requires the manifest.

## 4. Prepare and review the SLURM launcher

```bash
rnaseq-service-upstream-smoke prepare \
  --run-plan "$SMOKE_PRIVATE/run-plan.json" \
  --output "$SMOKE_PRIVATE/launch/smoke.sbatch" \
  --account "$SMOKE_ACCOUNT" \
  --partition "$SMOKE_PARTITION" \
  --time 08:00:00 \
  --cpus 8 \
  --memory-gb 32 \
  --module "$NEXTFLOW_MODULE" \
  --module "$APPTAINER_MODULE"

sed -n '1,240p' "$SMOKE_PRIVATE/launch/smoke.sbatch"
bash -n "$SMOKE_PRIVATE/launch/smoke.sbatch"
```

Preparation creates a mode-`0700` script but does not call `sbatch`. The two
stages run sequentially and both use `-resume`. A one-minute heartbeat identifies
the active stage and elapsed time. The selected route is recorded in
`runtime/network_selection.tsv` without recording proxy values.

## 5. Submit exactly once

```bash
rnaseq-service-hpc submit \
  --launcher "$SMOKE_PRIVATE/launch/smoke.sbatch" \
  --receipt "$SMOKE_PRIVATE/submission.json"
```

The submission command reserves a mode-`0600` receipt before calling `sbatch`,
then records the launcher checksum and scheduler job ID. It refuses to overwrite
the receipt, preventing an unnoticed duplicate submission.

## 6. Monitor scheduler and heartbeat progress

```bash
rnaseq-service-hpc status \
  --receipt "$SMOKE_PRIVATE/submission.json"

export SMOKE_JOB_ID="$(python -c \
  'import json,os; print(json.load(open(os.environ["SMOKE_PRIVATE"] + "/submission.json"))["job_id"])')"

squeue -j "$SMOKE_JOB_ID" -o '%.18i %.10T %.10M %.10l %.6D %R'
tail -f "$SMOKE_PRIVATE/launch/logs/smoke-$SMOKE_JOB_ID.out" \
        "$SMOKE_PRIVATE/launch/logs/smoke-$SMOKE_JOB_ID.err"
```

The controller log ends with `[upstream-smoke] COMPLETE` only after both
Nextflow commands return success. Scheduler success alone is insufficient: the
receipt gate below also validates the traces and scientific output contract.

## 7. Validate and seal the evidence

After `rnaseq-service-hpc status` reports a terminal successful job:

```bash
rnaseq-service-upstream-smoke complete \
  --run-plan "$SMOKE_PRIVATE/run-plan.json" \
  --output "$SMOKE_PRIVATE/completion.json"
```

The command shows artifact-by-artifact hashing progress. It refuses completion
unless both traces contain only successful terminal tasks, both `params.json`
files point to the planned output directories, and all of the following exist:

- Nextflow trace, report, timeline, DAG, parameters, and software versions for
  both stages;
- RNA-seq MultiQC, merged Salmon gene counts, and gene TPMs;
- differential-analysis report, normalized and VST matrices, DESeq2 full and
  filtered tables, volcano plot, and GSEA report.

The non-overwriting mode-`0600` receipt stores SHA-256 hashes of this evidence
and explicitly records that real scientific execution occurred on public
upstream test data with no client data. A missing or failed artifact leaves M7
incomplete; do not edit the receipt or copy old outputs into the run root.

## Upstream evidence sources

The required outputs follow the official pinned test configuration and snapshot
contracts:

- [`nf-core/rnaseq` 3.26.0 test configuration](https://github.com/nf-core/rnaseq/blob/3.26.0/conf/test.config)
- [`nf-core/rnaseq` 3.26.0 default test snapshot](https://github.com/nf-core/rnaseq/blob/3.26.0/tests/default.nf.test.snap)
- [`nf-core/differentialabundance` 2.0.0 RNA-seq DESeq2 snapshot](https://github.com/nf-core/differentialabundance/blob/2.0.0/tests/test_rnaseq_deseq2.nf.test.snap)
