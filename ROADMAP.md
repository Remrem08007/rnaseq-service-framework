# RNA-seq Service Framework Roadmap

Status: **M0–M5 complete; M6 planned**

## Goal

Deliver a generic, client-ready operational layer around established nf-core
RNA-seq workflows. Success means a study can move from validated intake to a
reproducible handoff on real HPC infrastructure without pretending that this
repository invented the underlying RNA-seq methods.

## Locked principles

- reuse and pin maintained community workflows;
- validate the scientific and file contracts before spending compute;
- separate pipeline completion from QC acceptance;
- preserve explicit contrast direction and sample-exclusion history;
- support local, SLURM, proxy-aware, and offline execution;
- keep credentials, institutional accounts, client paths, and client data out
  of the public repository;
- use synthetic fixtures for public CI;
- report limitations and avoid clinical claims.

## Milestones

### M0 — Contract and intake preflight

- service scope and non-goals;
- architecture and network/security policies;
- sample-sheet, design, and contrast schemas;
- dependency-light preflight CLI;
- unit tests and Python 3.11/3.12 CI.

### M1 — Reproducible launch plan

- parse the pinned workflow lock;
- generate inspectable nf-core commands without shell interpolation hazards;
- record parameters, versions, input hashes, and execution mode;
- support dry-run planning before submission.

### M2 — Direct, proxy, and offline staging

- preflight network and registry access;
- inherit but never persist proxy configuration;
- download pinned pipelines and Apptainer images;
- build a checksum-bound offline bundle;
- refuse missing artifacts and silent network fallback.

### M3 — HPC execution profiles

- generic SLURM/Apptainer configuration;
- external account, partition, storage, and resource settings;
- reference-index reuse;
- resumable driver and status commands;
- measured CPU, memory, time, and storage summaries.

### M4 — Primary RNA-seq delivery

- validated `nf-core/rnaseq` launch;
- completion and artifact checks;
- MultiQC collection;
- machine-readable provenance and run receipt;
- safe restart and failure diagnostics.

### M5 — QC acceptance gate

- configured mapping, strandedness, complexity, duplication, and outlier checks;
- explicit accept/review/stop decisions;
- immutable record of excluded samples and reasons;
- no automatic exclusion based only on an arbitrary threshold.

### M6 — Differential analysis handoff

- lock the accepted sample set and design;
- convert validated contrasts to pinned `nf-core/differentialabundance` inputs;
- run appropriate supported analysis profiles;
- verify expected tables, plots, enrichment, and report artifacts.

### M7 — Synthetic end-to-end validation

- deterministic paired-end fixture with known expression differences;
- direct test profile;
- pre-staged offline smoke test;
- restart/resume test;
- expected-direction and output-contract assertions.

### M8 — Client handoff and portfolio release

- client intake checklist;
- methods, QC, limitations, and delivery templates;
- example report using only public synthetic data;
- HPC runbook and troubleshooting guide;
- release tag and profile update.

## Definition of success

1. Invalid sample or contrast metadata fails before workflow launch.
2. Every delivered run names exact workflow and software revisions.
3. A restricted-network cluster can run from a verified offline bundle.
4. Pipeline success cannot be confused with QC acceptance.
5. Differential analysis uses a reviewed, locked design and sample set.
6. Public CI exercises synthetic data only.
7. A client receives results, provenance, exclusions, limitations, and rerun
   instructions as one coherent handoff.
