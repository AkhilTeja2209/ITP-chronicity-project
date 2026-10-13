# Pediatric ITP chronicity: integrated variant × expression analysis

[![pipeline](../../actions/workflows/pipeline.yml/badge.svg)](../../actions/workflows/pipeline.yml)

Most children with immune thrombocytopenia (ITP) recover within a year, but 20–30% develop chronic
disease, and nothing measured at diagnosis reliably tells them apart. This repository is a fully
reproducible pipeline that asks two questions:

1. **Expression layer:** which genes and co-expression modules in children's blood differ between
   self-limiting and chronic ITP, *after* removing artefacts (platelet-RNA loss, steroid exposure,
   cell composition, sex, batch)?
2. **Variant layer:** which of those genes also carry inherited evidence of involvement in ITP,
   through TWAS and colocalisation of ITP GWAS variants with blood and T-cell eQTLs?

It then tests whether the two layers converge more than chance, and benchmarks the resulting
signature against a baseline on a held-out cohort from a different platform.

**Latest results:** [`results/REPORT.md`](results/REPORT.md) (tables and figures alongside), rebuilt and committed by GitHub Actions on every run.

## Why this design (short version)

* The public pediatric ITP datasets are **microarrays**. They measure mRNA abundance and never read DNA,
  so GATK-style variant calling is impossible. Germline sequence of ITP patients exists only under
  controlled access (dbGaP/EGA). The variant layer therefore uses published **GWAS summary statistics**.
* Adults have the genetic data but not the clinical question (adult ITP is chronic by default).
  Children have the question but not the genetic data. Integrating the two layers is the workable middle.

Full rationale, methods and anticipated questions: [`docs/PROJECT_PLAN.md`](docs/PROJECT_PLAN.md).

## Data (all public, downloaded automatically)

| Layer | Source | What it contributes |
|---|---|---|
| Expression | [GSE23754](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23754): whole blood, Stanford cDNA array (GPL10845), 8 acute / 4 progression-to-chronic / 14 chronic | the only samples drawn **before** the outcome was known; used as the held-out evaluation cohort |
| Expression | [GSE46922](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46922): sorted T cells, GPL570, 7 newly diagnosed / 6 chronic | T-cell-specific chronic vs newly diagnosed contrast |
| Expression | [GSE56232](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56232): sorted T cells, GPL570, 6 children at diagnosis and at remission 12 months later | does the signal reverse on recovery (paired)? |
| Variant | [GWAS Catalog GCST90018865](https://www.ebi.ac.uk/gwas/studies/GCST90018865): ITP GWAS (Sakaue et al. 2021), 843 cases / 668k participants | inherited susceptibility |
| Translator | [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/): GTEx v10 whole blood, BLUEPRINT CD4⁺ T cells | links non-coding risk variants to the genes they regulate (colocalisation) |
| Translator | [PredictDB](https://predictdb.org/) GTEx v8 MASHR whole-blood models | genetically predicted expression (TWAS) |
| Pathways | MSigDB Hallmark + Reactome | enrichment |

### Data issues found by the audit step (not reported in earlier analyses)

* **GSE23754 is not Affymetrix.** It is a two-colour cDNA array (patient vs Universal Reference RNA). It
  cannot be RMA-normalised or ComBat-merged with GPL570, so it is analysed separately and used for
  cross-platform validation.
* **GSE46922 and GSE56232 share raw arrays.** Several "newly diagnosed" CEL files in GSE56232 are
  byte-identical to files in GSE46922, deposited under different patient numbers. The pipeline
  checksums every array and counts each physical array once (`results/tables/audit_duplicate_arrays.tsv`).
* GCST90018865 is a European + East Asian meta-analysis (843 cases), not the 675-case European-only figure
  that is often quoted.

## Pipeline

```
00 audit ── 01 expression ──┬─ 02 differential expression ── 03 modules + enrichment ─┐
            (RMA, QC,       │                                                           ├─ 20 integration ── 21 benchmark ── 30 report
             re-annotation, └─ 10 GWAS loci ──┬─ 11 TWAS ──────────────────────────────┤
             confounders)                     └─ 12 colocalisation ─────────────────────┘
```

| Step | Script | Key methods |
|---|---|---|
| 00 | `R/00_audit.R` | GEOquery; metadata-completeness table; platform technology; MD5 duplicate detection |
| 01 | `R/01_expression.R` | RMA (affy); RLE QC; GPL570 probes re-mapped with current Ensembl BioMart vs legacy `hgu133plus2.db`; cDNA clones re-mapped via Entrez/HGNC; platelet, glucocorticoid, interferon, cell-type and sex scores; ComBat only when batch is not confounded with group |
| 02 | `R/02_differential_expression.R` | limma, unadjusted / confounder-adjusted / steroid-sensitivity models; paired design for GSE56232; platelet-artefact check |
| 03 | `R/03_modules_enrichment.R` | WGCNA signed network; module scores projected into every cohort; fgsea / ORA |
| 04 | `R/04_novel_genes.R` | genes measurable only after re-annotation: hits, lncRNAs, robustness (probe support, expression, steroid model, reversal) |
| 10 | `R/10_gwas.R` | streamed QC, λ<sub>GC</sub>, distance-clumped loci, Manhattan/QQ |
| 11 | `R/11_twas.R` | S-PrediXcan formula on GTEx v8 MASHR whole-blood weights |
| 12 | `R/12_coloc.R` | `coloc.abf` per locus × gene × eQTL dataset, regions fetched with tabix |
| 20 | `R/20_integration.R` | four evidence filters per gene; Fisher overlap, cameraPR, Spearman and a-priori tests of convergence |
| 21 | `R/21_benchmark.R` | leave-one-cohort-out AUC on GSE23754 vs baseline; DeLong CIs, label-permutation and random-gene-set nulls; nested-LOOCV elastic net |
| 30 | `R/30_report.R` | writes `results/REPORT.md` |

All parameters are in [`config/config.yaml`](config/config.yaml); marker gene sets in
[`config/gene_sets.yaml`](config/gene_sets.yaml).

## Running it

**On GitHub (default).** Every push to `main` that touches code runs `.github/workflows/pipeline.yml`
(about 30 minutes on the first run, mostly package installation and downloads; both are cached afterwards).
To re-run by hand, use *Actions → pipeline → Run workflow*, optionally listing step ids such as `20 21`.
Outputs are uploaded as an artifact and committed to `results/` on `main`.

**Locally** (Linux/macOS; on Windows use WSL). Needs R ≥ 4.4 and Python 3 with `pysam` (`pip install pysam`), used for region queries of remote eQTL files:

```bash
Rscript -e 'install.packages("pak"); pak::local_install_deps()'
```

```bash
Rscript run_pipeline.R
```

About 3 GB of disk is needed for downloads. Arms are independent: if the variant arm fails (e.g. a
download server is down), the expression arm still completes and the report shows which steps ran.

## Interpreting the output honestly

Sample sizes are small (4 progressors in the only predictive contrast). The report generates
**[auto]** sentences from the numbers. When a test does not beat chance, the report says so. The
deliverable is a documented, confounder-aware, convergent-evidence candidate list plus an honest
benchmark, **not** a clinical predictor.

## Citations / acknowledgements

Data: GEO series GSE23754 (Stanford; PMID 21325602), GSE46922 and GSE56232 (Gothenburg; cite the
publications listed on each GEO record); Sakaue S. et al., *Nat Genet* 53:1415 (2021); Kerimov N. et al.,
*Nat Genet* 53:1290 (2021) (eQTL Catalogue); GTEx Consortium, *Science* 369:1318 (2020); Barbeira A.N.
et al., *Genome Biol* 22:49 (2021) (MASHR models) and *Nat Commun* 9:1825 (2018) (S-PrediXcan);
Giambartolomei C. et al., *PLoS Genet* 10:e1004383 (2014) (coloc); Liberzon A. et al. (MSigDB).

MIT licence (code only; each dataset keeps its own terms).
