# QC acceptance gate

M5 converts parsed MultiQC evidence into an explicit, reviewable sample decision
without treating a threshold as scientific truth. The output of M4 is
`pipeline_complete_qc_pending`; only a completed M5 review is
`qc_accepted_for_differential_analysis`.

## Safety properties

- The M4 schema-version-2 receipt stores size and SHA-256 for every parsed
  MultiQC file. M5 refuses changed evidence and does not consume a file added
  after completion.
- The QC policy is study-specific and checksum-bound into the assessment.
- Missing required files, columns, or sample values become
  `stop_missing_evidence`, never zero or pass.
- Threshold and outlier results are flags only. `automatic_exclusions` is
  always zero.
- Every sample needs an explicit human decision and reviewer. Every exclusion
  and every flagged inclusion needs a reason.
- Finalization rechecks the checksum-bound design and contrasts. It refuses an
  accepted set that falls below the run plan's minimum replicate count.
- Assessment and acceptance directories are private, non-overwriting outputs.

## 1. Create the M4 receipt

After a successful primary run:

```bash
rnaseq-service-primary complete \
  --run-plan private/plans/study-001.json \
  --output private/delivery/study-001.primary-completion-v2.json
```

An older schema-version-1 receipt is not sufficient because it did not bind
every parsed MultiQC file. Keep the old receipt as history and create a new,
distinct v2 filename; completion refuses overwrites.

## 2. Inspect the available evidence

The exact columns depend on enabled nf-core modules and MultiQC. Inspect the
HTML report scientifically, then inspect exported headers before choosing the
policy:

```bash
MQC=results/study-001/rnaseq/multiqc/star_salmon/multiqc_report_data

head -1 "$MQC/multiqc_general_stats.txt" | \
  tr '\t' '\n' | nl -ba

head -5 "$MQC/multiqc_strand_check_summary_table.txt" | \
  column -s $'\t' -t

head -5 "$MQC/multiqc_dupradar.txt" | \
  column -s $'\t' -t
```

If the run uses the legacy directory name, set `MQC` to `multiqc_data` instead.

## 3. Review a study policy

Start from the example but keep the engagement-specific copy private:

```bash
mkdir -p private/qc
cp examples/qc/policy.example.toml private/qc/study-001.policy.toml
```

Each `[[metrics]]` rule names an exact parsed TSV, sample column, and value
column. Numeric rules may set `lower`, `upper`, and `outlier = true`.
Categorical rules name `allowed_values` and use `outlier = false`.

The example numbers are operational examples, not universal biological
cutoffs. Choose mapping, duplication, complexity, and strandedness rules after
reviewing library type, organism, assay, sequencing depth, experimental design,
and the complete MultiQC report.

For numeric metric-profile outliers, M5 computes
`0.67448975 * (value - median) / MAD` across samples when enough values exist.
This robust univariate signal is useful for review but is not an expression-PCA
outlier test and does not prove a sample is bad.

## 4. Assess without excluding

```bash
rnaseq-service-qc assess \
  --completion-receipt private/delivery/study-001.primary-completion-v2.json \
  --policy private/qc/study-001.policy.toml \
  --outdir private/qc/study-001-assessment
```

The progress bar advances by configured metric. Outputs are:

- `qc_assessment.json`: evidence hashes, policy, flags, missing values, and
  study state;
- `qc_sample_metrics.tsv`: compact sample-by-metric review table;
- `qc_decisions.tsv`: template requiring human decisions.

Inspect the state and flags:

```bash
python - <<'PY'
import json
from pathlib import Path

p = Path("private/qc/study-001-assessment/qc_assessment.json")
q = json.loads(p.read_text())
print("status:", q["status"])
print("samples:", q["sample_count"])
print("flagged:", q["flagged_sample_count"])
print("missing required:", q["missing_required_sample_count"])
print("automatic exclusions:", q["automatic_exclusions"])
PY

column -s $'\t' -t \
  private/qc/study-001-assessment/qc_sample_metrics.tsv | less -S
```

Do not finalize `stop_missing_evidence`. Correct the pipeline/module or policy
mismatch, preserve the stopped assessment, then create a new assessment
directory.

## 5. Record human decisions

Edit `qc_decisions.tsv` so every row has:

| Column | Required value |
|---|---|
| `sample` | unchanged assessed sample ID |
| `decision` | exactly `include` or `exclude` |
| `reason` | required for exclusions and flagged inclusions |
| `reviewer` | non-empty reviewer identity for every sample |

The tool never infers a decision from the flags. Review MultiQC, sample history,
experimental design, and any client-provided facts before deciding.

## 6. Seal the accepted sample set

```bash
rnaseq-service-qc finalize \
  --assessment private/qc/study-001-assessment/qc_assessment.json \
  --decisions private/qc/study-001-assessment/qc_decisions.tsv \
  --outdir private/qc/study-001-accepted
```

The command creates:

- `qc_acceptance.json`: immutable provenance, all decisions and reasons,
  evidence/design/policy hashes, and per-contrast replicate counts;
- `accepted_samples.tsv`: the only sample set M6 may hand to differential
  analysis.

If an exclusion leaves too few reference or target replicates for any planned
contrast, finalization stops. That is a scientific-design problem requiring an
explicit revised scope or new data—not a reason to silently weaken the check.
