#!/usr/bin/env python3
"""Time-to-remission analysis of a local pediatric ITP cohort.

Usage:
    python local/cohort_analysis.py local/data/cohort.csv     # your de-identified data
    python local/cohort_analysis.py --demo                    # simulated children, to test the code

Outputs go to local/output/ (git-ignored). See local/README.md for the input columns.
Runs entirely on this computer; nothing is uploaded.
"""
from __future__ import annotations

import argparse
import base64
import html
import io
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import statsmodels.api as sm  # noqa: E402
from lifelines import CoxPHFitter, KaplanMeierFitter  # noqa: E402
from lifelines.statistics import logrank_test, proportional_hazard_test  # noqa: E402
from scipy import stats  # noqa: E402

OUT = Path(__file__).resolve().parent / "output"
REQUIRED = ["id", "sex", "age_at_dx_years", "dx_date", "symptom_days_before_dx", "plt_dx",
            "first_treatment", "remission_date", "last_followup_date"]
OPTIONAL = ["plt_1mo", "infection_4wk", "vaccine_4wk"]
# Cheng et al. 2023 (Diagnostics 13:341), 497 children with severe ITP
REFERENCE = {"remitted_by_3mo": 0.732, "chronic_12mo": 0.209}


# ---------------------------------------------------------------- input
def load(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path) if path.suffix.lower() in (".xlsx", ".xls") else pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        sys.exit(f"Missing required columns: {', '.join(missing)}. See local/cohort_template.csv.")
    for c in OPTIONAL:
        if c not in df.columns:
            df[c] = np.nan
    return df


def prepare(df: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    problems: list[str] = []
    d = df.copy()
    for c in ["dx_date", "remission_date", "last_followup_date"]:
        d[c] = pd.to_datetime(d[c], errors="coerce")
    d["event"] = d["remission_date"].notna().astype(int)
    end = d["remission_date"].where(d["event"] == 1, d["last_followup_date"])
    d["days"] = (end - d["dx_date"]).dt.days
    bad = d["days"].isna() | (d["days"] < 0)
    if bad.any():
        problems.append(f"{int(bad.sum())} rows dropped: missing dates, or remission/follow-up before diagnosis "
                        f"(ids: {', '.join(map(str, d.loc[bad, 'id'].head(10)))})")
    d = d[~bad].copy()
    d["days"] = d["days"].clip(lower=0.5)  # remission on the day of diagnosis
    d["age_lt10"] = (d["age_at_dx_years"] < 10).astype(int)
    d["abrupt_onset"] = (d["symptom_days_before_dx"] <= 14).astype(int)
    d["male"] = (d["sex"].astype(str).str.upper().str[0] == "M").astype(int)
    d["plt_dx_per10"] = d["plt_dx"] / 10
    d["plt_1mo_gt100"] = np.where(d["plt_1mo"].isna(), np.nan, (d["plt_1mo"] > 100).astype(float))
    d["first_treatment"] = d["first_treatment"].fillna("unknown").astype(str).str.lower().str.strip()
    # chronic at 12 months: no remission within 365 days, among children followed >= 365 days or remitted earlier
    d["chronic_12mo"] = np.select([(d["event"] == 1) & (d["days"] <= 365), d["days"] > 365], [0, 1], default=np.nan)
    known = d["chronic_12mo"].notna()
    if (~known).any():
        problems.append(f"{int((~known).sum())} children have < 12 months of follow-up without remission: "
                        "excluded from the 12-month chronicity analysis, kept (censored) in the time-to-event analysis")
    return d, problems


# ---------------------------------------------------------------- demo data
def simulate(n: int = 300, seed: int = 2026) -> pd.DataFrame:
    """Simulated children with hazards loosely based on published predictors. NOT real data."""
    rng = np.random.default_rng(seed)
    age = np.clip(rng.gamma(2.2, 2.6, n), 0.3, 17.9)
    symptom = rng.choice([1, 3, 7, 10, 21, 35, 60], n, p=[.25, .25, .2, .1, .1, .06, .04])
    plt_dx = np.clip(rng.lognormal(2.4, .7, n), 1, 29).round()
    infection = rng.binomial(1, .45, n)
    treat = rng.choice(["none", "ivig", "steroid", "steroid+ivig"], n, p=[.3, .35, .25, .1])
    # built-in effects (so the demo shows whether the models recover them): younger age and abrupt
    # onset raise the chance of remitting and shorten the time to it
    lin = 1.2 * (age < 10) + 1.0 * (symptom <= 14) + 0.3 * infection
    cure = rng.random(n) < 1 / (1 + np.exp(-(-0.4 + lin)))          # who ever remits within ~2 years
    t_rem = np.where(cure, rng.weibull(0.9, n) * 90 / np.exp(0.5 * lin), np.inf)
    followup = rng.uniform(200, 1500, n)
    dx = pd.Timestamp("2018-01-01") + pd.to_timedelta(rng.integers(0, 1500, n), unit="D")
    remitted = t_rem <= followup
    plt1 = np.where(t_rem < 30, rng.normal(180, 40, n), np.where(t_rem < 90, rng.normal(90, 40, n), rng.normal(45, 25, n)))
    return pd.DataFrame({
        "id": [f"SIM-{i:04d}" for i in range(n)], "sex": rng.choice(["M", "F"], n),
        "age_at_dx_years": age.round(1), "dx_date": dx.strftime("%Y-%m-%d"), "symptom_days_before_dx": symptom,
        "plt_dx": plt_dx, "plt_1mo": np.clip(plt1, 2, 400).round(), "infection_4wk": infection,
        "vaccine_4wk": rng.binomial(1, .05, n), "first_treatment": treat,
        "remission_date": [(d + pd.Timedelta(days=float(t))).strftime("%Y-%m-%d") if r else "" for d, t, r in zip(dx, t_rem, remitted)],
        "last_followup_date": [(d + pd.Timedelta(days=float(f))).strftime("%Y-%m-%d") for d, f in zip(dx, followup)],
    })


# ---------------------------------------------------------------- analyses
def km_table(d: pd.DataFrame) -> tuple[pd.DataFrame, KaplanMeierFitter]:
    km = KaplanMeierFitter().fit(d["days"], d["event"], label="all")
    rows = []
    for t, lab in [(90, "3 months"), (180, "6 months"), (365, "12 months")]:
        s = float(km.survival_function_at_times(t).iloc[0])
        ci = km.confidence_interval_survival_function_
        lo_s, hi_s = (float(np.interp(t, ci.index, ci.iloc[:, j])) for j in (0, 1))
        rows.append({"time": lab, "cumulative_remission": 1 - s, "ci_low": 1 - hi_s, "ci_high": 1 - lo_s,
                     "at_risk": int((d["days"] >= t).sum())})
    return pd.DataFrame(rows), km


def cox(d: pd.DataFrame) -> tuple[pd.DataFrame, str, int]:
    cols = ["days", "event", "age_lt10", "abrupt_onset", "plt_dx_per10", "male"]
    x = d[cols + ["infection_4wk", "first_treatment"]].copy()
    if x["infection_4wk"].notna().mean() >= 0.8:
        cols.append("infection_4wk")
    tdum = pd.get_dummies(x["first_treatment"], prefix="tx", drop_first=False).astype(float)
    if "tx_none" in tdum:
        tdum = tdum.drop(columns="tx_none")  # reference: no treatment
    else:
        tdum = tdum.iloc[:, 1:]
    x = pd.concat([x[cols], tdum], axis=1).dropna()
    x = x.loc[:, (x.std() > 0) | x.columns.isin(["days", "event"])]
    cph = CoxPHFitter(penalizer=0.01).fit(x, "days", "event")
    s = cph.summary[["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]].reset_index()
    s.columns = ["covariate", "hazard_ratio_remission", "ci_low", "ci_high", "p"]
    try:
        ph = proportional_hazard_test(cph, x, time_transform="rank").summary
        ph_note = "Proportional-hazards check (Schoenfeld residuals): " + ", ".join(
            f"{i} p={r['p']:.3f}" for i, r in ph.iterrows()) + ". p < 0.05 means that covariate's effect changes over time."
    except Exception as e:  # pragma: no cover
        ph_note = f"Proportional-hazards check could not run: {e}"
    return s, ph_note, len(x)


def landmark(d: pd.DataFrame) -> tuple[pd.DataFrame | None, str]:
    lm = d[(d["days"] > 30) & d["plt_1mo_gt100"].notna()].copy()
    if lm["event"].sum() < 10 or lm["plt_1mo_gt100"].nunique() < 2:
        return None, "Not enough children with a one-month platelet count still unremitted at one month."
    lm["days_from_1mo"] = lm["days"] - 30
    cph = CoxPHFitter().fit(lm[["days_from_1mo", "event", "plt_1mo_gt100"]], "days_from_1mo", "event")
    s = cph.summary[["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]].reset_index()
    s.columns = ["covariate", "hazard_ratio_remission", "ci_low", "ci_high", "p"]
    return s, f"{len(lm)} children still unremitted at one month with a platelet count recorded then."


def auc(y: np.ndarray, score: np.ndarray) -> float:
    pos, neg = score[y == 1], score[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    return float(stats.mannwhitneyu(pos, neg, alternative="two-sided").statistic / (len(pos) * len(neg)))


def chronic_logistic(d: pd.DataFrame) -> tuple[pd.DataFrame | None, str]:
    c = d[d["chronic_12mo"].notna()].dropna(subset=["age_lt10", "abrupt_onset", "plt_dx_per10"])
    if c["chronic_12mo"].sum() < 5:
        return None, "Fewer than 5 chronic children with 12 months of follow-up: logistic model not fitted."
    X = sm.add_constant(c[["age_lt10", "abrupt_onset", "plt_dx_per10"]].astype(float))
    fit = sm.Logit(c["chronic_12mo"].astype(float), X).fit(disp=0)
    ci = np.exp(fit.conf_int())
    s = pd.DataFrame({"covariate": X.columns, "odds_ratio_chronic": np.exp(fit.params).values,
                      "ci_low": ci[0].values, "ci_high": ci[1].values, "p": fit.pvalues.values})
    # leave-one-out AUC so the figure is not optimistic
    pred = np.empty(len(c))
    for i in range(len(c)):
        m = np.ones(len(c), bool); m[i] = False
        try:
            fi = sm.Logit(c["chronic_12mo"].astype(float)[m], X[m]).fit(disp=0)
            pred[i] = float(fi.predict(X.iloc[[i]]).iloc[0])
        except Exception:
            pred[i] = np.nan
    ok = ~np.isnan(pred)
    a = auc(c["chronic_12mo"].values[ok].astype(int), pred[ok])
    return s[s.covariate != "const"], f"{len(c)} children, {int(c['chronic_12mo'].sum())} chronic. Leave-one-out AUC = {a:.2f} (0.5 = chance)."


def compare_reference(d: pd.DataFrame, kmt: pd.DataFrame) -> pd.DataFrame:
    rows = []
    c = d[d["chronic_12mo"].notna()]
    k, n = int(c["chronic_12mo"].sum()), len(c)
    if n:
        rows.append({"measure": "chronic at 12 months", "this_cohort": k / n, "n": n,
                     "cheng_2023": REFERENCE["chronic_12mo"], "binomial_p": stats.binomtest(k, n, REFERENCE["chronic_12mo"]).pvalue})
    r3 = kmt.loc[kmt.time == "3 months", "cumulative_remission"].iloc[0]
    rows.append({"measure": "remitted by 3 months (Kaplan-Meier)", "this_cohort": r3, "n": len(d),
                 "cheng_2023": REFERENCE["remitted_by_3mo"], "binomial_p": np.nan})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- report
def fig_to_b64(fig) -> str:
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=150, bbox_inches="tight"); plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def km_figure(d: pd.DataFrame, simulated: bool) -> tuple[str, str]:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    notes = []
    for ax, (col, labels, title) in zip(axes, [("age_lt10", ("age ≥ 10", "age < 10"), "By age at diagnosis"),
                                                ("abrupt_onset", ("onset > 14 days", "abrupt onset ≤ 14 days"), "By onset")]):
        for v, lab in zip((0, 1), labels):
            g = d[d[col] == v]
            if len(g) < 5:
                continue
            km = KaplanMeierFitter().fit(g["days"], g["event"], label=f"{lab} (n={len(g)})")
            ax.step(km.survival_function_.index, 1 - km.survival_function_.iloc[:, 0], where="post", label=km._label)
        lr = logrank_test(d.loc[d[col] == 0, "days"], d.loc[d[col] == 1, "days"],
                          d.loc[d[col] == 0, "event"], d.loc[d[col] == 1, "event"])
        notes.append(f"{title}: log-rank p = {lr.p_value:.3g}")
        ax.set_title(title); ax.set_xlim(0, 730); ax.set_xlabel("days from diagnosis"); ax.axvline(365, ls=":", c="grey")
        ax.legend(frameon=False, fontsize=8)
    axes[0].set_ylabel("cumulative remission")
    if simulated:
        fig.suptitle("SIMULATED DATA (code test only)", color="red")
    return fig_to_b64(fig), "; ".join(notes)


def table_html(df: pd.DataFrame | None) -> str:
    if df is None or df.empty:
        return "<p><em>Not available.</em></p>"
    f = df.copy()
    for c in f.columns:
        if f[c].dtype.kind == "f":
            f[c] = f[c].map(lambda v: "" if pd.isna(v) else (f"{v:.3g}" if abs(v) < 0.001 or c == "p" or c.endswith("_p") else f"{v:.2f}"))
    return f.to_html(index=False, escape=True, border=0, classes="t")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data", nargs="?", type=Path, help="cohort .csv or .xlsx")
    ap.add_argument("--demo", action="store_true", help="run on simulated children to test the code")
    a = ap.parse_args()
    if not a.demo and not a.data:
        ap.error("give a data file, or --demo")
    simulated = a.demo
    raw = simulate() if simulated else load(a.data)
    d, problems = prepare(raw)
    prefix = "demo_" if simulated else ""
    OUT.mkdir(exist_ok=True)

    kmt, _ = km_table(d)
    cox_t, ph_note, n_cox = cox(d)
    lm_t, lm_note = landmark(d)
    lg_t, lg_note = chronic_logistic(d)
    ref_t = compare_reference(d, kmt)
    img, lr_note = km_figure(d, simulated)
    events = int(d["event"].sum())
    epv = events / max(1, len(cox_t))
    if epv < 10:
        problems.append(f"Only {epv:.1f} remission events per Cox covariate (10 or more is advisable): treat hazard ratios as unstable.")

    for name, t in [("km", kmt), ("cox", cox_t), ("landmark", lm_t), ("logistic_chronic", lg_t), ("reference", ref_t)]:
        if t is not None:
            t.to_csv(OUT / f"{prefix}cohort_{name}.csv", index=False)

    banner = ('<p class="sim">SIMULATED DATA. These numbers test the code and describe no real children. '
              'Do not cite or present them.</p>') if simulated else ""
    page = f"""<!doctype html><meta charset="utf-8"><title>{'SIMULATED ' if simulated else ''}ITP cohort report</title>
<style>body{{font:15px/1.5 Calibri,Arial,sans-serif;max-width:60rem;margin:2rem auto;padding:0 1rem;color:#1f1a1c}}
h1{{font-family:Cambria,serif}} .sim{{background:#fde8e8;color:#9e1b32;padding:.6rem 1rem;font-weight:bold}}
table.t{{border-collapse:collapse;margin:.5rem 0 1rem}} table.t td,table.t th{{border-bottom:1px solid #ddd;padding:.25rem .6rem;text-align:left}}
.note{{color:#555}}</style>
<h1>{'SIMULATED: ' if simulated else ''}Time to remission in pediatric ITP</h1>{banner}
<p>{len(raw)} children in the file, {len(d)} analysable, {events} remissions, median follow-up {d['days'].median():.0f} days.</p>
{''.join(f'<p class="note">⚠ {html.escape(p)}</p>' for p in problems)}
<h2>Cumulative remission</h2>{table_html(kmt)}
<img src="data:image/png;base64,{img}" style="max-width:100%"><p class="note">{html.escape(lr_note)}</p>
<h2>Comparison with Cheng et al. 2023</h2>{table_html(ref_t)}
<h2>Cox model of time to remission</h2><p class="note">Hazard ratio &gt; 1 means faster remission. {n_cox} children with complete covariates. Reference treatment: none.</p>
{table_html(cox_t)}<p class="note">{html.escape(ph_note)}</p>
<h2>One-month landmark: platelets &gt; 100 ×10⁹/L at one month</h2><p class="note">{html.escape(lm_note)}</p>{table_html(lm_t)}
<h2>Chronic at 12 months (logistic)</h2><p class="note">{html.escape(lg_note)}</p>{table_html(lg_t)}
"""
    (OUT / f"{prefix}cohort_report.html").write_text(page, encoding="utf-8")
    print(f"Wrote {OUT / (prefix + 'cohort_report.html')}")
    print(kmt.to_string(index=False))


if __name__ == "__main__":
    main()
