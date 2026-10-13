# Project plan and rationale

## 1. Question

At the time a child is diagnosed with ITP, is there a molecular signal in blood that distinguishes
those who will recover (≈70–80%) from those who will become chronic? And do the genes carrying that
signal also carry *inherited* evidence of involvement in ITP?

**Deliverables**
1. A documented audit of whether public data can answer the question, including what cannot be done.
2. A confounder-aware set of chronicity-associated genes and modules.
3. A test of whether inherited ITP susceptibility converges on those genes, which helps separate cause from consequence.
4. An honest benchmark against a baseline on a held-out cohort.

## 2. Why the original "variant analysis on GSE46922" plan changed

| Fact | Consequence |
|---|---|
| GSE46922 / GSE56232 are Affymetrix GPL570 expression arrays; GSE23754 is a two-colour cDNA array | None of them contains DNA sequence, so variants cannot be called (no GATK, no Funcotator) |
| Germline sequence of patients is identifiable and is held under controlled access (dbGaP/EGA), not in GEO | No public pediatric ITP dataset with chronicity labels has callable reads (audit result: **0 datasets**) |
| ITP GWAS summary statistics are open | The variant layer starts from already-called variants; tools are TWAS + colocalisation, not read processing |

## 3. The two layers

**Expression layer: what the immune cells are doing now.**
Normalise each platform separately → re-annotate probes to the current genome → estimate confounders
from expression → limma per cohort → WGCNA modules → project modules across cohorts.

**Variant layer: what the child was born with.**
GWAS loci → TWAS (genetically predicted blood expression ↔ ITP risk) → colocalisation (is it the
*same* causal variant driving both the GWAS and the eQTL signal, or two neighbours in LD?).

**Intersection.** Four filters per gene:

| Filter | Dataset | Question |
|---|---|---|
| F1 | GSE23754 | Does it differ *before* the outcome is known (progression vs acute, both <6 months)? |
| F2 | GSE46922 | Does it differ inside T cells (chronic vs newly diagnosed)? |
| F3 | GSE56232 | Does it reverse when the same children remit? |
| F4 | GWAS + eQTL | Is there inherited evidence (coloc PP.H4 > 0.5 or TWAS FDR < 0.1)? |

DNA is fixed at conception, so it cannot have been changed by the illness. A chronicity gene that
also passes F4 is more plausibly upstream of the disease course than a reaction to it.

## 4. Confounders, and why they are modelled

| Confounder | Mechanism | Handling |
|---|---|---|
| Platelet RNA | Whole blood carries platelet mRNA; thrombocytopenia lowers it, so PF4/PPBP look "down" by circular reasoning | platelet transcript score is a covariate; artefact check on top genes |
| Steroids | Glucocorticoid receptor directly induces FKBP5, TSC22D3, DDIT4… | steroid score; sensitivity model |
| Infection | Interferon signature; a preceding infection predicts recovery | interferon score reported |
| Cell composition | Shifts in cell proportions look like regulation | marker-based cell scores (contamination in sorted sets) |
| Sex | Possible sex-specific chronicity markers | inferred from XIST / chrY genes |
| Batch | Lab / run differences | ComBat only if not confounded with group; otherwise cross-cohort validation instead of pooling |
| Duplicate arrays | Same physical arrays deposited in two series | MD5 checksums; each array counted once |

## 5. Benchmark

Signatures are derived only from the GPL570 T-cell cohorts and evaluated on GSE23754 (a different
tissue and platform): a true leave-one-cohort-out test. The baseline uses an expression-derived platelet score
and inferred sex, because GEO provides no per-sample age, sex or platelet count. Each AUC is reported
with a DeLong CI, a label-permutation p-value and a random-gene-set p-value. With 4 progressors, a CI
spanning 0.5 is the expected outcome. Reporting it is part of the result, not a failure.

## 6. Anticipated questions

**"Are you separating children's variants from adults', or chronicity variants from non-chronicity variants?"**
Neither, directly. No public dataset has chronicity-labelled genotypes. The design asks whether
*susceptibility* genes (from adult-dominated GWAS) overlap *chronicity-associated expression* in children.
That is a proxy, and the report says so.

**"Why not do it in adults, where variant data exists?"**
Adults have the data but not the question: adult ITP is chronic by default, with roughly 9% spontaneous remission.
Children have the question but not the data. The only adult longitudinal transcriptome (GSE112278)
measures eltrombopag response, which has already been analysed with DEG/WGCNA/hub-gene methods.

**"Is the prediction-model idea new?"**
Clinical predictors are established: abrupt onset, age < 10, platelets > 100×10⁹/L at one month. So are flow-cytometry
models (108 children). VNN1 and a CNR2 variant have been proposed and never validated. Genome-wide transcriptomic prediction
with an integrated genetic layer has not been published.

**"A parent had childhood ITP that resolved; could the child inherit it?"**
That is a familial question. It belongs to the inherited thrombocytopenia field (ThromboGenomics / NGS gene panels,
WES in families), because up to 40% of "chronic ITP" with a family history turns out to be an inherited
thrombocytopenia. It needs family DNA sequencing and is out of scope. Individual risk is a question for a clinician.

**"Is n enough?"** Not for a clinical model. It is enough for hypothesis-generating signature discovery with
cross-cohort validation and an explicit limitation.

**"What if nothing converges?"** The audit, the duplicate-array finding, the confounder analysis and a
well-powered negative benchmark are still a reportable methods-and-limitations result.

## 7. Future work (needs new data)

* Time-to-remission (Cox) modelling: needs transcriptomes linked to follow-up dates. A retrospective
  chart review plus prospective sampling at a tertiary haematology centre could build this cohort.
* Sustained remission off TPO-RA in adults: explicitly unsolved; flow cytometry at discontinuation.
* FinnGen R12 ITP endpoint as GWAS replication once access is granted (placeholder `gwas.replication_url` in config; not yet implemented).
