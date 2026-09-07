# SLURM and Apptainer execution

M3 turns a reviewed M1 run plan into a resumable SLURM controller job while
keeping cluster-specific accounts, partitions, modules, and paths outside the
public repository.

## Execution model

The controller itself is one small `sbatch` job. Inside that allocation,
Nextflow runs with `process.executor = 'slurm'` and submits each pipeline task
as its own scheduler job. The generated configuration also enables Apptainer,
sets nf-core task ceilings, and limits the number of simultaneously queued
tasks.

This separation matters:

- launcher resources cover the Nextflow controller only;
- task resources continue to come from the nf-core process labels, bounded by
  the configured maxima;
- the partition is rendered through Nextflow's `queue` directive;
- the project account is rendered as a SLURM `clusterOptions` value;
- `-resume` is mandatory, so a new controller job can reuse successful task
  outputs in the same work directory.

## 1. Create private infrastructure settings

Copy the generic example into an ignored private directory:

```bash
mkdir -p private/hpc
cp examples/hpc/infrastructure.example.toml private/hpc/infrastructure.toml
chmod 600 private/hpc/infrastructure.toml
```

Edit every placeholder. The settings define:

- SLURM account, partition, and controller job name;
- controller CPUs, memory, and wall time;
- maximum CPUs, memory, and wall time for individual nf-core tasks;
- Nextflow queue size and per-CPU memory behavior;
- absolute shared work, Apptainer cache, and optional iGenomes reference paths;
- module names required on the target cluster.

The validator accepts only constrained scheduler/module tokens and absolute
paths. It never guesses an institution, account, partition, or filesystem.

## 2. Render the Nextflow configuration

```bash
rnaseq-service-hpc render-config \
  --settings private/hpc/infrastructure.toml \
  --output private/hpc/nextflow.config
```

The output is created with mode `0600` and is never overwritten. Review it
before planning. To change infrastructure settings, create a new settings and
output path so the earlier execution contract remains attributable.

## 3. Create the SLURM run plan

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

For a direct or proxy-aware run, omit `--offline-manifest` and select the
corresponding network mode. The work directory must be on storage visible to
the controller and every compute node.

## 4. Prepare the immutable controller script

```bash
rnaseq-service-hpc prepare \
  --run-plan private/plans/study-001.json \
  --settings private/hpc/infrastructure.toml \
  --output private/launch/study-001/controller.sbatch
```

Preparation re-hashes every plan control file, confirms the Nextflow config is
the exact rendering of the supplied settings, requires a SLURM plan containing
`-resume`, checks that the work directory is under `work_root`, and enforces
`NXF_OFFLINE=true` for offline plans. The script and log directory are created
before submission; the script is mode `0700`.

Proxy values are never written into the launcher. `#SBATCH --export=ALL`
allows a deliberately configured submission environment to be inherited. For
offline operation, only the non-secret bundle paths and offline flag are
explicitly exported.

## 5. Submit exactly once

```bash
rnaseq-service-hpc submit \
  --launcher private/launch/study-001/controller.sbatch \
  --receipt private/launch/study-001/submission.json
```

The command reserves the mode-`0600` receipt before calling
`sbatch --parsable`. It records the launcher SHA-256 and scheduler job ID. An
existing receipt is refused, including a receipt from a failed submission, so a
second submission cannot happen silently.

## 6. Check status

```bash
rnaseq-service-hpc status \
  --receipt private/launch/study-001/submission.json
```

Active jobs are reported from `squeue` with state, elapsed time, limit, and
node/reason. Once absent from the queue, the command uses `sacct` to report the
terminal state, exit code, elapsed time, controller CPU time, and maximum RSS.
This controller accounting does not replace task-level Nextflow metrics.

## 7. Summarize task resources

The M1 plan enables a Nextflow trace at
`results/study-001/execution/trace.tsv`. After completion:

```bash
rnaseq-service-hpc summarize \
  --trace results/study-001/execution/trace.tsv \
  --json-out results/study-001/execution/resource_summary.json \
  --tsv-out results/study-001/execution/resource_by_process.tsv \
  --storage work=/shared/project/work/study-001 \
  --storage results=results/study-001
```

The command shows row progress and groups tasks by process. It reports task,
completed, cached, and failure counts; maximum observed resident memory; total
real and scheduler duration; and mean observed CPU percentage. The JSON also
records logical and filesystem-allocated bytes, file/directory counts, and
ignored symlink counts for each requested storage root. Trace rows have a
percentage bar; storage walks show a live file counter. Outputs are mode `0600`
and never overwritten.

Nextflow trace metrics are estimates. Very short tasks and task images missing
standard process-inspection tools can have incomplete measurements. Missing
metric counts are therefore explicit instead of being silently converted to
zero.

Storage measurements do not follow symlinks and can count hard-linked content
more than once. Allocated bytes come from filesystem block accounting and may
differ from storage-system quotas or snapshots.

## Reference reuse

When `paths.reference_cache` is present, the renderer sets nf-core's
`params.igenomes_base` to that shared absolute path. M4 requires an explicit,
reviewed `--genome` key in every new CLI-generated plan and checks the completed
pipeline used the same key. Custom FASTA/GTF and custom index manifests are not
yet supported by the wrapper.
