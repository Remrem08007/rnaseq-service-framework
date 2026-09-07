# Direct, proxy, and offline staging

M2 separates network access from client-data processing. Staging downloads the
two pinned nf-core workflows and their Apptainer images on an approved
networked system. Sealing then hashes the resulting files. Restricted compute
nodes only receive and verify that sealed bundle.

## Locked tools and workflows

`config/workflows.toml` pins four independent versions:

- `nf-core/rnaseq` 3.26.0;
- `nf-core/differentialabundance` 2.0.0;
- Nextflow 26.04.6;
- nf-core/tools 4.1.0.

The staging runner checks both installed tool versions before downloading. A
workflow or tool update requires a reviewed lock change and new bundle.

## 1. Install the staging tools

Use a networked login or transfer node where institutional policy permits
downloads. Java 17 or newer, Apptainer (or Singularity), Nextflow 26.04.6, and
nf-core/tools 4.1.0 must be available.

```bash
python -m venv .staging-venv
source .staging-venv/bin/activate
python -m pip install 'nf-core==4.1.0'

nextflow -version
nf-core --version
apptainer --version
```

Install the official `nextflow-26.04.6-dist` standalone executable when the
target environment cannot bootstrap Nextflow from the internet.

## 2. Create a non-executing staging plan

Direct network:

```bash
rnaseq-service-stage plan \
  --workflow-lock config/workflows.toml \
  --bundle-dir "$PWD/private/offline-bundle" \
  --output "$PWD/private/staging-plan.json" \
  --network-mode direct \
  --parallel-downloads 4
```

Proxy-aware network:

```bash
export HTTPS_PROXY='http://USER:PASSWORD@proxy.example:8080'
export HTTP_PROXY="$HTTPS_PROXY"
export NO_PROXY='localhost,127.0.0.1'

rnaseq-service-stage plan \
  --workflow-lock config/workflows.toml \
  --bundle-dir "$PWD/private/offline-bundle" \
  --output "$PWD/private/staging-plan.json" \
  --network-mode proxy \
  --parallel-downloads 4
```

The plan records only which standard proxy variable names are present. It
never records their values, host names, usernames, or passwords. Keep proxy
configuration outside the repository and unset it after staging.

## 3. Execute staging once

Review `private/staging-plan.json`, then perform the read-only tool and network
check:

```bash
rnaseq-service-stage check \
  --plan "$PWD/private/staging-plan.json" \
  --timeout-seconds 10
```

This verifies the exact nf-core/tools and Nextflow versions, finds Apptainer or
Singularity, and checks reachability of the GitHub API, Quay, and Docker
registry. HTTP authentication responses count as reachable; connection and
proxy failures do not. The JSON report contains no proxy values.

Once `"ready": true`, run:

```bash
rnaseq-service-stage run \
  --plan "$PWD/private/staging-plan.json" \
  --heartbeat-seconds 30
```

The command uses argv vectors rather than a shell, displays the nf-core download
progress plus an elapsed-time heartbeat, redacts exact proxy values from
forwarded output, and creates `private/staging-plan.receipt.json`. A receipt is
reserved before downloads begin; an existing target or receipt is refused.

Both workflows share
`private/offline-bundle/containers` through `NXF_SINGULARITY_CACHEDIR`, which
avoids downloading duplicate images. Their pinned Nextflow plugins are staged
under `private/offline-bundle/plugins` through `NXF_PLUGINS_DIR`. The completed
receipt records the exact downloaded workflow directories discovered from
`main.nf` and `nextflow.config`.

## 4. Seal the completed bundle

Use the two `workflow_dir` values from the receipt:

```bash
rnaseq-service-bundle seal \
  --bundle-dir "$PWD/private/offline-bundle" \
  --workflow-lock config/workflows.toml \
  --rnaseq-workflow /absolute/path/from/rnaseq/receipt \
  --differential-workflow /absolute/path/from/differential/receipt \
  --container-root "$PWD/private/offline-bundle/containers" \
  --plugin-root "$PWD/private/offline-bundle/plugins"
```

Sealing shows byte-level progress while hashing large images and writes
`offline_bundle.manifest.json` with mode `0600`. It refuses symbolic links,
missing workflow entry points, an empty image cache, or an existing manifest.
The manifest declares that staging artifacts contain neither client data nor
secrets.

## 5. Transfer and verify offline

Transfer the entire directory without adding or removing files. On the target
cluster:

```bash
rnaseq-service-bundle verify \
  --manifest /path/to/offline-bundle/offline_bundle.manifest.json \
  --workflow-lock config/workflows.toml
```

Verification checks exact file membership, every byte count and SHA-256 digest,
the component layout, and all workflow/runtime pins. Missing, unexpected, or
modified files fail verification.

## 6. Create an offline analysis plan

```bash
rnaseq-service-plan \
  --samplesheet private/intake/samplesheet.csv \
  --input-manifest private/inputs/study-001.fastq-manifest.json \
  --genome GRCh38 \
  --design private/intake/design.csv \
  --contrasts private/intake/contrasts.csv \
  --output private/plans/study-001.json \
  --outdir results/study-001 \
  --workdir work/study-001 \
  --network-mode offline \
  --offline-manifest /path/to/offline-bundle/offline_bundle.manifest.json
```

Offline planning refuses to proceed without successful bundle verification. It
uses the local RNA-seq workflow directory, omits the remote `-r` lookup, and
records the required `NXF_OFFLINE=true`, `NXF_SINGULARITY_CACHEDIR`, and
`NXF_PLUGINS_DIR` settings.
Nextflow documents that `NXF_OFFLINE=true` blocks remote project updates and
plugin downloads; therefore plugin versions in staged workflows must also be
explicit.

## Security boundary

- Never commit a staging plan, receipt, bundle, proxy configuration, or client
  path.
- Do not place FASTQs, sample metadata, reference genomes, or results inside a
  software bundle.
- A successful checksum verification proves bundle integrity, not scientific
  validity or QC acceptance.
- If staging fails after the receipt is reserved, inspect the cause and create
  a fresh plan and bundle path. Do not erase the failure boundary silently.
