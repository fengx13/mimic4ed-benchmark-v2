# MIMIC-IV-ED 5-Outcome Benchmark

Reproduces the **MultiAutoScore** master ED dataset — the same cohort, the
same 269 predictor columns (demographics, triage and serial vital signs,
five-year comorbidity history, laboratory panel, severity scores) — trimmed
to five clinical outcomes instead of the full project's forty-nine:

| Outcome | Column | Definition |
|---|---|---|
| Acute coronary syndrome | `outcome_cp_acs` | ICD-9/10 prefix match (STEMI + NSTEMI + unstable angina + other acute ischemic heart disease) |
| Sepsis | `outcome_sepsis` | ICD-9/10 exact match |
| Acute heart failure (AHF) | `outcome_hf_acute` | ICD-9/10 exact match (acute / acute-on-chronic only) |
| Stroke | `outcome_stroke` | ICD-9/10 exact match (non-traumatic cerebrovascular events) |
| Pulmonary embolism | `outcome_cp_pe` | ICD-9/10 prefix match |

Every outcome is positive if the code appears in the ED visit's own diagnosis
list **or** the linked hospital admission's discharge diagnosis list — see
[Outcome definitions](#outcome-definitions) for the exact code lists and the
clinical reasoning behind exact-vs-prefix matching.

**This repository's output was verified column-by-column against the full
49-outcome master dataset: all 269 shared predictor columns and all 5
outcome columns are value-identical.** The pipeline here is not an
approximation — it is the same code, with the outcome-generation step cut
down to five conditions.

## What you get

One row per MIMIC-IV-ED stay, **274 columns**:
demographics and timing, core ED outcomes (hospitalization, ICU transfer,
mortality, 3-day revisit), prior utilization counts, triage vitals (cleaned,
median-imputed, with missingness flags) plus first/last serial vitals during
the ED stay, ten chief-complaint flags, 48 Charlson/Elixhauser comorbidity
flags from five years of prior admissions, six bedside severity scores
(CART, REMS, NEWS, NEWS2, MEWS, Charlson index), a ~90-assay laboratory
panel, five-year diagnosis history, and the five outcomes above.

## Data access

You need **credentialed PhysioNet access** to both source datasets — this
repository contains no patient data, only code:

1. **MIMIC-IV v3.1** — <https://physionet.org/content/mimiciv/3.1/>
   (`hosp/` and `icu/` modules)
2. **MIMIC-IV-ED v2.2** — <https://physionet.org/content/mimic-iv-ed/2.2/>
   (`ed/` module)

Both require a signed PhysioNet Data Use Agreement and a completed human
subjects research training course. Row-level output from this pipeline is
equally restricted — do not share it outside your own DUA, and do not send
it to any third-party service.

## Setup

```
mimic-iv-ed-benchmark/
├── data/
│   ├── ed/      <- MIMIC-IV-ED v2.2 files go here (edstays.csv.gz, triage.csv.gz, ...)
│   ├── hosp/    <- MIMIC-IV v3.1 hosp/ files go here
│   └── icu/     <- MIMIC-IV v3.1 icu/ files go here
├── derived/     <- pipeline writes its outputs here
└── scripts/
```

Download the two datasets and either copy or symlink their `ed/`, `hosp/`
and `icu/` directories into `data/`:

```bash
ln -s /path/to/your/mimic-iv-3.1/ed   data/ed
ln -s /path/to/your/mimic-iv-3.1/hosp data/hosp
ln -s /path/to/your/mimic-iv-3.1/icu  data/icu
```

Install dependencies:

```bash
pip install -r requirements.txt
```

## Running the pipeline

Stages run in order; each reads the previous stage's checkpoint from `derived/`.

```bash
python3 scripts/01_cohort.py       # cohort spine, core outcomes, utilization
python3 scripts/02_vitals.py       # triage + serial vitals, chief complaints
python3 scripts/03_comorbidity.py  # Charlson/Elixhauser, 5-year dx history
python3 scripts/04_labs_extract.py # scans labevents.csv.gz (~2.4 GB, ~2 min)
python3 scripts/05_outcomes.py     # the 5 benchmark outcomes + severity scores
python3 scripts/06_labs_impute.py  # lab panel pivot, vital imputation, rescoring
python3 scripts/07_assemble.py     # final column order, write outputs
```

`04_labs_extract.py` and `05_outcomes.py` are independent of each other and
can run in either order or in parallel; `06_labs_impute.py` needs both to
have finished. Total runtime is roughly 5–7 minutes on a laptop, dominated by
the lab extraction scan.

**Output**: `derived/master_benchmark5.csv.gz` (full), `master_benchmark5_slim.csv.gz`
(without the four `dx_*_5y` history strings — loads faster), `master_benchmark5.pkl`
(preserves dtypes), and `column_dictionary.csv` (dtype/missingness/range per column).

## Outcome definitions

### ACS and PE — prefix match (`scripts/cp_codes.py`)

A diagnosis matches if its dot-less ICD code **starts with** any listed
prefix. Every prefix was individually expanded against MIMIC's
`hosp/d_icd_diagnoses.csv.gz` and read for clinical fit.

| | ICD-9-CM prefixes | ICD-10-CM prefixes |
|---|---|---|
| ACS | `410`, `4111`, `4118` | `I21`, `I22`, `I200`, `I248`, `I249` |
| PE | `4151`, `6732` | `I26`, `O882` |

### Sepsis, AHF, stroke — exact match (`scripts/outcome_codes_exact.py`)

Only the literal codes listed count — no prefix expansion. This mirrors how
these three were originally defined: `icd_codes.csv` in the source project
enumerates every code it intends for these three phenotypes, and exact
matching reproduces that file's own recorded event counts.

- **Sepsis**: 21 ICD-9 codes / 29 ICD-10 codes. Septicaemia, sepsis, severe
  sepsis, septic shock, and obstetric/abortion-related sepsis.
- **AHF**: 6 ICD-9 codes / 8 ICD-10 codes. Only codes explicitly labelled
  *acute* or *acute-on-chronic* heart failure — chronic/unspecified heart
  failure is deliberately excluded.
- **Stroke**: 15 ICD-9 codes / 161 ICD-10 codes. Non-traumatic cerebrovascular
  events (subarachnoid hemorrhage, intracerebral hemorrhage, cerebral
  infarction). All 161 ICD-10 codes were checked by automated keyword match
  against their official titles with zero misses.

### A coding error this benchmark does NOT inherit

An earlier version of the chest-pain prefix set included ICD-9 codes `8602`
and `8603` ("traumatic hemothorax") under a *pneumothorax* outcome — not one
of this benchmark's five, but discovered during an audit of the same
`cp_codes.py` file this benchmark imports. ICD-9-CM groups hemothorax and
pneumothorax into the same numeric category (860) for historical reasons,
but hemothorax is blood in the pleural space with no air involved — a
different condition. The fix (removing those two codes from the
pneumothorax definition) is already applied in `scripts/cp_codes.py` as
shipped here; it does not affect ACS or PE, which were separately audited
and had no such issue.

## Known limitations

- **GCS does not exist anywhere in MIMIC-IV-ED.** Verified across every ED
  table. `triage_gcs_first`/`triage_gcs_last` are 100% NaN; REMS and MEWS
  assume an alert patient; `score_NEWS2` is numerically identical to
  `score_NEWS` for the same reason (no supplemental-oxygen data either).
- **Labs are not triage-time information.** Median first lab result lands
  about an hour after ED arrival. Every `lab_*` column is post-triage; a
  model trained on them is not a triage-moment model.
- **Troponin and D-dimer are sparse by clinical design, not by accident.**
  These are ordered only when a clinician already suspects ACS or PE — in
  this cohort, `lab_troponin_t_blood` is non-null for about 5.6% of all
  stays. Treat "was this lab drawn" as an informative feature in its own
  right, not a value to impute.
- **ACS and PE use the broader cp_codes.py composite**, not the narrower
  exact-match single-code definitions used for sepsis/AHF/stroke — this is
  deliberate (see above), not an inconsistency to fix.
- **Comorbidity flags are inpatient-admission-only.** A zero Charlson/
  Elixhauser flag means "no prior inpatient record in this database," not
  "no comorbidity." About 52% of stays have no prior 5-year inpatient
  history at all.

## Reproducibility

This repository's output was verified **column-by-column against the full
49-outcome master dataset**: all 269 shared predictor columns and all 5
outcome columns are value-identical (0 mismatches). Given the same MIMIC-IV
v3.1 + MIMIC-IV-ED v2.2 files, anyone running this pipeline should get the
same `master_benchmark5.csv.gz` you would get from the full build, restricted
to these 5 outcomes.

Two things that *could* affect reproducibility if not matched exactly:

- **Dataset version.** This was tested against MIMIC-IV v3.1 and
  MIMIC-IV-ED v2.2 specifically. An earlier (or later) release can have
  different row counts, column names, or even data-correctness fixes (for
  example, MIMIC-IV v2.0 fixed an ED-length-of-stay bug present in v1.0) —
  using a different version will not silently fail, but it will not produce
  identical numbers either.
- **Package versions.** `requirements.txt` pins the exact pandas/numpy
  versions this was built and verified against. Several lab/vital
  extraction steps sort rows with duplicate keys (1.7% of lab rows share an
  exact `stay_id`/`itemid`/`timestamp`) before taking the first value; this
  was confirmed safe because pandas dispatches multi-column `sort_values`
  to `numpy.lexsort`, which is stable regardless of the `kind` argument or
  pandas version — ties resolve to "preserve the source CSV's row order,"
  which is identical for any official download of the same file. Using an
  unpinned, much older or newer pandas is still expected to work, but
  hasn't been separately tested.

## Citation

If you use this benchmark, cite the underlying datasets:

> Johnson, A., Bulgarelli, L., Pollard, T., Horng, S., Celi, L. A., & Mark,
> R. (2024). MIMIC-IV (version 3.1). PhysioNet.
>
> Johnson, A., Bulgarelli, L., Pollard, T., Horng, S., Celi, L. A., & Mark,
> R. (2023). MIMIC-IV-ED (version 2.2). PhysioNet.
