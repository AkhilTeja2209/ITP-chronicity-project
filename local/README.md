# Local analyses (run on your own computer, never in CI)

Two analyses need data that must not leave your institution. They run locally with Python and
write their outputs to `local/output/`. Both `local/data/` and `local/output/` are git-ignored:
**never commit patient data**, and use only de-identified records covered by your ethics approval.

Setup (once):

```bash
python -m pip install --user -r local/requirements.txt
```

## 1. Time-to-remission cohort (`cohort_analysis.py`)

Put a de-identified spreadsheet (`.csv` or `.xlsx`) at `local/data/cohort.csv`, one row per child,
with the columns in `cohort_template.csv`:

| Column | Meaning | Required |
|---|---|---|
| `id` | study number (not a hospital number) | yes |
| `sex` | `M` / `F` | yes |
| `age_at_dx_years` | age at diagnosis, years (decimals fine) | yes |
| `dx_date` | date of diagnosis, `YYYY-MM-DD` | yes |
| `symptom_days_before_dx` | days of bleeding symptoms before diagnosis | yes |
| `plt_dx` | platelet count at diagnosis, ×10⁹/L | yes |
| `plt_1mo` | platelet count about one month after diagnosis, ×10⁹/L | if known |
| `infection_4wk` | viral illness in the 4 weeks before onset, `1`/`0` | if known |
| `vaccine_4wk` | vaccination in the 4 weeks before onset, `1`/`0` | if known |
| `first_treatment` | `none`, `steroid`, `ivig`, `steroid+ivig`, `other` | yes |
| `remission_date` | first date of remission (platelets ≥ 100 ×10⁹/L without treatment, IWG 2009); blank if never | yes |
| `last_followup_date` | last visit date | yes |

```bash
python local/cohort_analysis.py local/data/cohort.csv
```

It reports: cumulative remission at 3, 6 and 12 months (Kaplan–Meier, censored at last follow-up);
a Cox model of time to remission on age < 10, abrupt onset (≤ 14 days), platelet count at
diagnosis, sex, preceding infection and first treatment; a one-month landmark analysis of
platelets > 100 at one month (children already in remission or without one month of follow-up
are excluded, to avoid immortal-time bias); and a comparison with Cheng et al. 2023
(73% remission by 3 months, 21% chronic). Output: `local/output/cohort_report.html` plus tables.

To test the code before you have data: `python local/cohort_analysis.py --demo` runs it on
**simulated** children and labels every output as simulated.

## 2. Validate candidate genes (`validate_candidates.py`)

Checks whether the T-cell signature and the robust newly measurable genes replicate in an
independent expression dataset (RNA-seq or arrays, ideally sorted T cells).

* `local/data/expression.csv`: genes in rows (first column gene symbol or Ensembl ID), samples in
  columns, log-scale values (log2 CPM / log2 intensity).
* `local/data/samples.csv`: columns `sample`, `group` (`ND` and `chronic`, or `ND` and
  `remission`), optional `pair` for paired designs (same child at two time points).

```bash
python local/validate_candidates.py local/data/expression.csv local/data/samples.csv
```

It compares each candidate's direction with the discovery result, runs a sign test, scores the
discovery signature in the new samples (AUC with bootstrap CI) and writes
`local/output/validation_report.html`.

To test it without data: `--demo` plants 40% of the discovery effect in simulated samples (it should
replicate), and `--demo --demo-effect 0` plants none (it should not). Both outputs are labelled simulated.
