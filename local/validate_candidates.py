#!/usr/bin/env python3
"""Check the discovery T-cell signature and the robust newly measurable genes in an independent dataset.

Usage:
    python local/validate_candidates.py local/data/expression.csv local/data/samples.csv
    python local/validate_candidates.py --demo          # simulated samples, to test the code

expression.csv : first column gene (symbol or Ensembl ID), one column per sample, log-scale values
samples.csv    : columns sample, group (ND + chronic, or ND + remission), optional pair
Discovery results are read from results/tables/ in this repository. Outputs go to local/output/.
"""
from __future__ import annotations

import argparse
import html
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
TABLES = ROOT / "results" / "tables"
OUT = Path(__file__).resolve().parent / "output"
CONTRAST = {frozenset({"ND", "chronic"}): ("GSE46922", "chronic_vs_ND", "chronic"),
            frozenset({"ND", "remission"}): ("GSE56232", "remission_vs_ND", "remission")}


def discovery(contrast: str, cohort: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    de = pd.read_csv(TABLES / "de_adjusted_all.tsv", sep="\t")
    de = de[(de.cohort == cohort) & (de.contrast == contrast) & (de.model == "adjusted")]
    nov = pd.read_csv(TABLES / "novel_gene_hits.tsv", sep="\t")
    nov = nov[(nov.contrast == contrast) & (nov.robust == True)]  # noqa: E712
    return de, nov


def harmonise(expr: pd.DataFrame, de: pd.DataFrame) -> pd.DataFrame:
    """Index the new dataset by discovery gene_id, whether it uses symbols or Ensembl IDs."""
    ids = expr.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    if (ids.str.startswith("ENSG")).mean() > 0.5:
        expr.index = ids
    else:
        sym2id = de.dropna(subset=["symbol"]).drop_duplicates("symbol").set_index("symbol")["gene_id"]
        expr.index = ids.map(sym2id)
    expr = expr[expr.index.notna()]
    return expr.groupby(level=0).mean()


def per_gene(expr: pd.DataFrame, s: pd.DataFrame, case: str) -> pd.DataFrame:
    a = expr[s.loc[s.group == case, "sample"]]
    b = expr[s.loc[s.group == "ND", "sample"]]
    if "pair" in s and s["pair"].notna().all():
        pairs = s.pivot(index="pair", columns="group", values="sample").dropna()
        a, b = expr[pairs[case]], expr[pairs["ND"]]
        diff = a.values - b.values
        t, p = stats.ttest_1samp(diff, 0, axis=1, nan_policy="omit")
        lfc = np.nanmean(diff, axis=1)
    else:
        t, p = stats.ttest_ind(a, b, axis=1, equal_var=False, nan_policy="omit")
        lfc = a.mean(axis=1).values - b.mean(axis=1).values
    return pd.DataFrame({"logFC_new": lfc, "t_new": np.asarray(t), "p_new": np.asarray(p)}, index=expr.index)


def score_auc(expr: pd.DataFrame, s: pd.DataFrame, case: str, up: list[str], down: list[str], rng) -> dict:
    z = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1).replace(0, np.nan), axis=0)
    up, down = [g for g in up if g in z.index], [g for g in down if g in z.index]
    sc = (z.loc[up].mean() if up else 0) - (z.loc[down].mean() if down else 0)
    keep = s[s.group.isin([case, "ND"])]
    y = (keep.group == case).astype(int).values
    x = sc[keep["sample"]].values

    def auc(yy, xx):
        pos, neg = xx[yy == 1], xx[yy == 0]
        return stats.mannwhitneyu(pos, neg).statistic / (len(pos) * len(neg))

    a = auc(y, x)
    boot = []
    for _ in range(2000):
        i = rng.integers(0, len(y), len(y))
        if 0 < y[i].sum() < len(y):
            boot.append(auc(y[i], x[i]))
    perm = np.mean([auc(rng.permutation(y), x) >= a for _ in range(2000)])
    return {"genes_up": len(up), "genes_down": len(down), "AUC": a, "ci_low": np.percentile(boot, 2.5),
            "ci_high": np.percentile(boot, 97.5), "permutation_p": perm, "n_case": int(y.sum()), "n_ND": int((1 - y).sum())}


def simulate(de: pd.DataFrame, case: str, rng) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Simulated dataset: 8 vs 8 samples, the discovery effects present at 40% of their size. NOT real data."""
    genes = de.gene_id.values
    n = 8
    base = rng.normal(6, 1.5, len(genes))[:, None] + rng.normal(0, 0.6, (len(genes), 2 * n))
    eff = np.where(de["adj.P.Val"].values < 0.1, 0.4 * de.logFC.values, 0)
    base[:, :n] += eff[:, None]
    cols = [f"SIM_{case}_{i}" for i in range(n)] + [f"SIM_ND_{i}" for i in range(n)]
    expr = pd.DataFrame(base, index=genes, columns=cols)
    samples = pd.DataFrame({"sample": cols, "group": [case] * n + ["ND"] * n})
    return expr, samples


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("expression", nargs="?", type=Path)
    ap.add_argument("samples", nargs="?", type=Path)
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    rng = np.random.default_rng(2026)
    simulated = a.demo
    if simulated:
        cohort, contrast, case = CONTRAST[frozenset({"ND", "chronic"})]
        de, nov = discovery(contrast, cohort)
        expr, s = simulate(de, case, rng)
    else:
        if not (a.expression and a.samples):
            ap.error("give expression.csv and samples.csv, or --demo")
        s = pd.read_csv(a.samples)
        if not {"sample", "group"} <= set(s.columns):
            sys.exit("samples.csv needs columns sample and group")
        key = frozenset(s.group.unique())
        if key not in CONTRAST:
            sys.exit(f"groups must be ND + chronic or ND + remission, got {sorted(key)}")
        cohort, contrast, case = CONTRAST[key]
        de, nov = discovery(contrast, cohort)
        expr = pd.read_csv(a.expression, index_col=0)
        missing = set(s["sample"]) - set(expr.columns)
        if missing:
            sys.exit(f"samples not in expression.csv: {sorted(missing)[:10]}")
        expr = harmonise(expr[s["sample"]], de)

    pg = per_gene(expr, s, case)
    m = de.set_index("gene_id")[["symbol", "logFC", "t", "adj.P.Val"]].join(pg, how="inner")
    m.index.name = "gene_id"
    sig = m[m["adj.P.Val"] < 0.1]
    rows = []
    for name, sub in [(f"discovery T-cell signature (FDR < 0.1, {cohort})", sig),
                      ("robust newly measurable genes", m.loc[m.index.intersection(nov.gene_id)])]:
        sub = sub.dropna(subset=["logFC_new"])
        k = int((np.sign(sub.logFC) == np.sign(sub.logFC_new)).sum())
        rows.append({"gene_set": name, "genes_in_discovery": len(de[de["adj.P.Val"] < 0.1]) if "signature" in name else len(nov),
                     "genes_measured_here": len(sub), "same_direction": k,
                     "fraction_same_direction": k / len(sub) if len(sub) else np.nan,
                     "sign_test_p": stats.binomtest(k, len(sub), 0.5, alternative="greater").pvalue if len(sub) else np.nan,
                     "replicated_p05_same_direction": int(((sub.p_new < 0.05) & (np.sign(sub.logFC) == np.sign(sub.logFC_new))).sum())})
    conc = pd.DataFrame(rows)
    rho, rho_p = stats.spearmanr(m["t"], m["t_new"], nan_policy="omit")
    sc = score_auc(expr, s, case, list(sig[sig.logFC > 0].index), list(sig[sig.logFC < 0].index), rng)

    OUT.mkdir(exist_ok=True)
    prefix = "demo_" if simulated else ""
    conc.to_csv(OUT / f"{prefix}validation_concordance.csv", index=False)
    m.reset_index().to_csv(OUT / f"{prefix}validation_per_gene.csv", index=False)
    pd.DataFrame([sc]).to_csv(OUT / f"{prefix}validation_signature_auc.csv", index=False)
    novel = m.loc[m.index.intersection(nov.gene_id)].reset_index()

    banner = '<p class="sim">SIMULATED DATA. These numbers test the code only. Do not cite them.</p>' if simulated else ""
    f = lambda df: df.to_html(index=False, border=0, classes="t", float_format=lambda v: f"{v:.3g}")  # noqa: E731
    page = f"""<!doctype html><meta charset="utf-8"><title>{'SIMULATED ' if simulated else ''}Candidate validation</title>
<style>body{{font:15px/1.5 Calibri,Arial,sans-serif;max-width:60rem;margin:2rem auto;padding:0 1rem}} .sim{{background:#fde8e8;color:#9e1b32;padding:.6rem 1rem;font-weight:bold}}
table.t{{border-collapse:collapse}} table.t td,table.t th{{border-bottom:1px solid #ddd;padding:.25rem .6rem}}</style>
<h1>{'SIMULATED: ' if simulated else ''}Validation of the {contrast.replace('_', ' ')} signature</h1>{banner}
<p>{sc['n_case']} {case} vs {sc['n_ND']} newly diagnosed samples; {len(m):,} genes matched to the discovery analysis.
Genome-wide Spearman correlation of t statistics with discovery: rho = {rho:.3f} (p = {rho_p:.2g}).</p>
<h2>Direction concordance</h2>{f(conc)}
<h2>Signature score</h2><p>Discovery FDR genes scored as mean z(up) − mean z(down); AUC with bootstrap 95% CI.</p>{f(pd.DataFrame([sc]))}
<h2>Robust newly measurable genes</h2>{f(novel[['gene_id', 'symbol', 'logFC', 'logFC_new', 'p_new']]) if len(novel) else '<p>None measured in this dataset.</p>'}
"""
    (OUT / f"{prefix}validation_report.html").write_text(page, encoding="utf-8")
    print(conc.to_string(index=False)); print(pd.DataFrame([sc]).to_string(index=False))
    print(f"Wrote {OUT / (prefix + 'validation_report.html')}")


if __name__ == "__main__":
    main()
