# Service contract

## Purpose

The framework prepares and supervises reproducible research-use bulk RNA-seq
delivery. It uses pinned community workflows for primary processing and
differential analysis, and supplies the validation, infrastructure, provenance,
quality-control, and handoff layer required for a client engagement.

It does not claim that a successful workflow run proves that the study design is
scientifically valid. Automated checks support, but do not replace, analyst and
investigator review.

## Initial supported analysis

- bulk RNA-seq from an organism with a reference genome and annotation;
- single-end or paired-end FASTQ input, including multiple sequencing lanes;
- one row per lane in the upstream-compatible sample sheet;
- categorical fixed-effect comparisons with at least two biological replicates
  in each requested level;
- pinned `nf-core/rnaseq` primary processing;
- pinned `nf-core/differentialabundance` downstream analysis;
- local or SLURM execution with Apptainer;
- direct, proxy-aware staging, and fully pre-staged offline operation.

Complex blocking, paired, longitudinal, interaction, continuous-covariate, and
mixed-effects designs require explicit statistical review. They are not silently
reduced to a simple two-group comparison.

## Required intake artifacts

### RNA-seq sample sheet

CSV with the four upstream-required columns:

| Column | Meaning |
|---|---|
| `sample` | Stable biological sample identifier; repeated rows represent lanes |
| `fastq_1` | Read 1 or single-end FASTQ |
| `fastq_2` | Read 2 FASTQ; blank for single-end data |
| `strandedness` | `auto`, `forward`, `reverse`, or `unstranded` |

The framework rejects reused FASTQs, duplicate rows, invalid identifiers, mixed
single/paired layout within one sample, and missing files when filesystem checks
are enabled.

### Design table

CSV with exactly one row per biological sample. The first column is `sample` and
remaining columns are analysis variables such as `condition`, `batch`, or
`sex`. All RNA-seq samples must appear exactly once, and unknown samples are
rejected.

### Contrast table

CSV with:

| Column | Meaning |
|---|---|
| `contrast_id` | Stable output-safe comparison identifier |
| `variable` | Design-table column being compared |
| `reference` | Denominator/reference level |
| `target` | Numerator/target level |

Reference and target direction is explicit. The framework never infers the
scientific comparison from alphabetical factor ordering.

## Delivery gates

1. **Intake gate:** identifiers, files, read layout, metadata, and contrasts pass
   preflight.
2. **Design gate:** an analyst confirms experimental unit, replication,
   covariates, confounding, and contrast direction.
3. **Primary-processing gate:** workflow completion and MultiQC outputs exist.
4. **QC gate:** mapping, library complexity, strandedness, sample similarity,
   and outlier findings are reviewed before differential analysis.
5. **Differential-analysis gate:** the locked design and accepted sample set are
   recorded before modeling.
6. **Handoff gate:** versions, parameters, checksums, reports, exclusions,
   limitations, and reproducibility instructions are delivered together.

## Client inputs and analyst responsibilities

The client supplies data-access authorization, sample metadata, experimental
context, intended comparisons, and any required data-retention constraints. The
analyst validates the computational inputs, reviews the design, records all
exclusions and transformations, executes the pinned workflows, interprets QC,
and communicates limitations.

## Exclusions

The initial service does not provide clinical diagnosis, clinical-grade variant
interpretation, automatic causal claims, unsupervised removal of inconvenient
samples, or guaranteed statistical power. Human-subject and controlled-access
data remain subject to their governing agreements.
