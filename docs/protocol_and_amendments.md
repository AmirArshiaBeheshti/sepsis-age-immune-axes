# Analysis protocol and amendments (verbatim texts from the analysis code)

The authoritative registered files are the `.md` files and their `.json` records (UTC time stamp and SHA-256 hash) that the scripts write to `SEPSIS_immuneage/preregistration/`. Copy them into this folder (see `protocol/README.md`) and check them with `tools/verify_protocol_hashes.py`. The texts below are reproduced from the code for reading convenience.

Notes: protocol v1 kept its original title ("Pre-registration v1"); it was not deposited in a public registry before the analyses. Hypotheses H1 and H3 of v1 (bulk immune-age clocks) were withdrawn by Amendment 3 before any outcome analysis, and v1 H2 became H1 of the manuscript. No Amendment 2 governs any reported analysis. Amendment 5 (raw single-cell counts) was superseded by Amendment 6 before any drug result was computed.

## Protocol v1

```text
# Pre-registration v1 - Sepsis immune-age decomposition study

## Rationale and prior art (checked before registration)
TraMA (Klopack et al., Aging 2025) applied a bulk transcriptomic mortality-age clock to one sepsis cohort
(GSE185263): sepsis patients appeared ~10 years older, and the association with mortality weakened to P = 0.09
after adjustment for cell type. Epigenetic PhenoAge acceleration was associated with sepsis non-survival (Crit Care
Explor 2024). Single-cell immune ageing clocks (Nature Aging 2025; IDEAL-Age, Genome Biol 2026) were not applied to
sepsis or to mortality. No study has separated cell-intrinsic immune ageing from leukocyte composition and
inflammation in multiple independent sepsis cohorts.

## Cohorts
  - GSE65682: role=discovery; array; whole blood; outcome=28-day mortality; MARS; ICU admission; includes non-infectious ICU
  - E-MTAB-4421: role=validation; array; leukocytes; outcome=28-day mortality; GAinS CAP discovery
  - E-MTAB-4451: role=validation; array; leukocytes; outcome=28-day mortality; GAinS CAP validation
  - E-MTAB-5273: role=validation; array; leukocytes; outcome=28-day mortality; GAinS faecal peritonitis
  - E-MTAB-5274: role=validation; array; leukocytes; outcome=28-day mortality; GAinS faecal peritonitis
  - GSE63042: role=validation; rnaseq; whole blood; outcome=28-day mortality; CAPSOD; ED presentation
  - GSE95233: role=validation; array; whole blood; outcome=hospital survival; septic shock; ICU day 1
  - GSE54514: role=validation; array; whole blood; outcome=28-day mortality; daily samples; first sample used
  - GSE185263: role=validation; rnaseq; whole blood; outcome=in-hospital mortality; already analysed by TraMA (Aging 2025)
  - GSE33828: role=healthy_reference; array; whole blood; outcome=none; Rotterdam Study healthy reference (age in supplementary file)

## Inclusion rules (applied identically to every cohort)
1. Adult patients (age >= 18) with sepsis / septic shock / sepsis due to CAP or faecal peritonitis as defined
   by the source study; healthy and non-infectious ICU samples are kept as separate groups.
2. One sample per patient: the earliest available time point (ties: first sample in the series).
3. Known age and known death label (0/1). Outcome definitions differ (28-day / hospital / in-hospital) and are kept
   per cohort; cohort-specific estimates are pooled by random-effects meta-analysis.
4. A validation cohort is analysed if it has >= 30 eligible patients and >= 10 deaths.
5. Expression: probes -> HGNC approved symbols (probes mapping to >1 gene removed; highest-mean probe per gene);
   Ensembl ids -> HGNC; arrays on log2 scale; RNA-seq counts -> log2(CPM + 1). All scores use genes standardised
   within each cohort among the eligible sepsis patients.

## Pre-specified scores
Gene sets (frozen; a score = mean within-cohort z of available genes; computed only if >= 60% of the
genes are measured):
  - immature_neutrophil: CEACAM8, MMP8, LCN2, LTF, CAMP, OLFM4, DEFA4, BPI, CRISP3, CD24
  - mature_neutrophil: FCGR3B, CXCR2, CSF3R, FPR1, MME, CMTM2, PROK2, FFAR2
  - ms1_monocyte: RETN, ALOX5AP, IL1R2, CD163, SLC39A8
  - classical_monocyte: CD14, LYZ, FCN1, VCAN, CSF1R, CD68
  - t_cell: CD3E, CD3D, CD3G, CD247, CD2, LCK, CD28, IL7R, CCR7, TRAT1
  - b_cell: MS4A1, CD79A, CD79B, CD19, CD22, BANK1, FCRLA, TCL1A
  - nk_cell: GNLY, NKG7, KLRD1, PRF1, KLRF1, GZMB, SH2D1B
  - mhc_class_ii: HLA-DRA, HLA-DRB1, HLA-DPA1, HLA-DPB1, HLA-DMA, HLA-DMB, CD74
  - inflammation: MSigDB Hallmark 'Inflammatory Response' (Enrichr library MSigDB_Hallmark_2020, frozen copy saved).
Composite composition index (H2): myelopoiesis-lymphopenia index = z(immature_neutrophil) - z(t_cell).

## Immune-age clocks (Stage 2, trained only on healthy donors)
- Single-cell clocks: OneK1K (training) and AIDA (external test); donor-level pseudobulk per cell type
  (CD4 T, CD8 T, NK, B, CD14 monocyte, CD16 monocyte); elastic net (reference) and an attention-based deep
  multiple-instance model over cells (comparison). Prior gene set = genes age-associated within a cell type with
  the same sign in OneK1K and AIDA (FDR < 0.05 in both).
- Bulk clocks trained in the healthy whole-blood reference (GSE33828), 10-fold cross-validation:
  'naive clock' = elastic net on all genes; 'intrinsic clock' = elastic net restricted to the single-cell prior
  genes after regressing every gene on the composition scores (cell-intrinsic ageing).
- Fidelity rule: a bulk clock is used only if its cross-validated Pearson r with age >= 0.50 in GSE33828;
  otherwise every hypothesis that uses it is reported as 'not testable'.
- Age acceleration in a sepsis cohort = residual of the clock prediction regressed on chronological age within the
  eligible patients of that cohort, standardised (per SD).

## Primary hypotheses (two-sided alpha 0.05, Holm correction over H1-H3, untestable hypotheses count as P = 1)
H1  Intrinsic immune-age acceleration is associated with higher death risk, adjusted for chronological age and sex.
    Model: Cox (time-to-event where available) or logistic regression (binary label). Discovery GSE65682 sepsis
    patients; validation = random-effects meta-analysis (DerSimonian-Laird) of the validation cohorts.
    Supported if the validation meta-estimate is > 1 with P < 0.05 (Holm) and the discovery estimate has the same sign.
H2  The myelopoiesis-lymphopenia index is associated with higher death risk, adjusted for chronological age and sex
    (same models, same supported rule).
H3  Composition explains the naive transcriptomic-age signal: the log effect of naive-clock age acceleration on death
    is attenuated by >= 50% after adding the composition scores (immature neutrophil, mature neutrophil, MS1,
    classical monocyte, T, B, NK, MHC-II). Estimate in the pooled validation data (cohort-stratified);
    95% CI by 2000 patient-level bootstrap. Supported if the attenuation estimate >= 50% and its 95% CI excludes 0.

## Secondary / descriptive analyses (no multiplicity claim)
S1  Decomposition: joint model with intrinsic age acceleration, composition scores, MHC-II and inflammation;
    likelihood-ratio contribution of each block (discovery and pooled validation).
S2  Specificity: the H1 model in GSE65682 non-infectious ICU patients (if death labels exist).
S3  Comparison with source-study endotypes (MARS1-4, SRS1/2) where deposited.
S4  Cell-type attribution in sepsis single-cell data (GSE167363 survivors vs non-survivors; COMBAT; Reyes 2020):
    single-cell clocks per cell type.
S5  Drug repurposing (exploratory). Signature = genes associated with death in GSE65682 adjusted for age, sex and
    the composition scores (FDR < 0.05), replicated in direction in the pooled validation cohorts (P < 0.05).
    Fidelity gates (decided now, after the TGCT experience): >= 10 genes per direction; >= 50% of genes per
    direction measured directly in L1000 (landmark genes); signature-level specificity: the number of reversing
    signatures must exceed the 95th percentile of 100 size-matched random signatures - otherwise the analysis is
    reported as 'not interpretable'. Queries: SigCom LINCS l1000_cp (all cell lines and an immune/haematopoietic
    subset: THP1, HL60, U937, JURKAT, NALM6), L1000CDS2, Enrichr drug libraries. Candidate = >= 2 reversing
    signatures in >= 2 cell lines and empirical BH-FDR < 0.10. Drugs tested in sepsis trials (hydrocortisone,
    dexamethasone, simvastatin, ascorbic acid, thiamine) are reported descriptively, not as positive controls.

## Rules for deviations
Corrections of metadata parsing are allowed only before any outcome analysis (Stage 3), as a dated amendment file.
All negative and untestable results are reported.
```

## Amendment 1

```text
# Amendment 1 to Pre-registration v1 - cohort-specific metadata rules

Written after the metadata field audit of Stage 1 and BEFORE any outcome analysis (only field names, value counts and
numbers of deaths per cohort had been seen). Reason: the deposited fields differ between cohorts; the generic
rules of v1 misclassified groups. Hypotheses, models, gene sets and gates of v1 are unchanged.

1. GSE65682 (MARS): sepsis patients = samples with a MARS endotype (endotype_class Mars1-4), i.e. the sepsis cohort
   of Scicluna et al. (Lancet Respir Med 2017), for which 28-day mortality and time-to-event are deposited.
   Healthy = icu_acquired_infection 'healthy'. Non-infectious ICU = abdominal_sepsis_and_controls 'ctrl_GI' or
   pneumonia diagnoses 'no-cap' (no deposited outcome -> secondary analysis S2 not testable).
2. GAinS: E-MTAB-5273/5274 contain community-acquired-pneumonia patients that are also in E-MTAB-4421/4451 (same
   patient identifiers). To keep validation cohorts independent, only faecal-peritonitis patients are used from
   E-MTAB-5273/5274; any faecal-peritonitis identifier also present in E-MTAB-4421/4451 is removed. Patient = sample
   identifier without the '.B.<day>' suffix; time point = the day suffix (first sample used).
   E-MTAB-4451 sampled survivors and non-survivors in equal numbers (outcome-dependent sampling): logistic models
   only (odds ratios remain valid), no absolute risks.
3. GSE95233: sepsis = samples with survival 'Survivor' / 'Non Survivor'; healthy = no survival label; patient = the
   identifier in the sample title (e.g. 'SS_500'); time point = 'time point' (D01 < D02 < D03); sex 'H' (French
   'homme') = male.
4. GSE54514: patient = group_id; time point = day in group_day; outcome = disease status (sepsis survivor /
   nonsurvivor).
5. GSE33828 (healthy reference): expression from the deposited raw matrix (RS3-ExpressionMatrix-881-raw-data),
   log2 and quantile normalisation; age / sex from the deposited RS3-881-age-gender file.
```

## Amendment 3

```text
# Amendment 3 to Pre-registration v1 - revised primary hypotheses (written before ANY outcome analysis)

Status when written: Stage-1 data preparation finished (cohort sizes, numbers of deaths, age distributions and
field audits seen). No association between any expression score and death, and no clock, had been computed.
Reason: the healthy whole-blood reference GSE33828 has a narrow age distribution (IQR 54-63 years), so the
pre-registered bulk clocks (v1 H1, H3) were judged infeasible; they are withdrawn before any training (Amendment 2,
if present, is superseded). A literature review (Nov 2026) found no study decomposing the age-mortality gradient in
sepsis through the blood transcriptomic immune state (closest: plasma-protein mediation in COVID-19, Michels et al.,
Eur Respir J 2023).

Cohorts, inclusion rules, gene sets, scores and the myelopoiesis-lymphopenia (ML) index are unchanged (v1, Amendment 1).
Scores: within-cohort z of the gene-set score among eligible patients; ML index = z(immature_neutrophil) - z(t_cell),
re-standardised within cohort. Age per 10 years.

PRIMARY (Holm over 2 hypotheses; v1 H1 and H3 withdrawn untested)
H1 (= v1 H2, unchanged): the ML index is associated with higher death risk, adjusted for age and sex.
   Discovery GSE65682: Cox (28-day time-to-event); validation cohorts: logistic regression; DerSimonian-Laird
   random-effects meta-analysis of validation log-ORs. Supported if the validation estimate > 1 with Holm P < 0.05
   and the discovery estimate has the same sign.
H2: part of the effect of older age on death is mediated by the blood immune state.
   Mediators (joint): ML index, MHC-II score, inflammation (Hallmark) score. Regression-based g-computation
   (VanderWeele): linear mediator models M ~ age + sex; logistic outcome model death ~ age + M + age x M + sex
   (age x M terms omitted in cohorts with < 40 deaths to keep >= 5 events per parameter); bootstrap fits that do
   not converge are discarded and counted;
   counterfactual mediators = observed value + age coefficient x 10 years; natural direct / indirect effects and
   total effect on the risk-difference scale for 10 years of age; 1000 patient-level bootstraps per cohort.
   E-MTAB-4451 (outcome-dependent sampling) is excluded from risk-scale mediation. Discovery = GSE65682; validation
   = random-effects meta-analysis of the natural indirect effect (NIE). Supported if the validation NIE > 0 with
   Holm P < 0.05 and the discovery NIE has the same sign. Proportion mediated = pooled NIE / pooled total effect
   (95% CI by Monte-Carlo from the pooled estimates). The decomposition is statistical (single time point);
   reverse causation cannot be excluded and is stated as a limitation.

SECONDARY / SENSITIVITY (descriptive)
- each immune score separately as H1-type predictor and as single mediator;
- age x ML-index interaction on death; age -> score associations (meta-analysed);
- leave-one-cohort-out; leave-E-MTAB-4451-out; platform-stratified; Cox HR in GSE65682 vs logistic OR;
- severity-adjusted models where severity is deposited (GSE185263 SOFA, GSE54514 APACHE II), noting that severity
  may lie on the causal path;
- comparison with deposited endotypes (MARS1-4, SRS1/2).

LATER STAGES (pre-specified now, run after this stage)
- Healthy single-cell immune ageing (OneK1K, AIDA; donor x cell-type pseudobulk, OLS age models): concordance of
  the healthy-ageing direction of the mediator genes with their sepsis age and mortality directions.
- Ligand-receptor analysis (LIANA+, donor-level, consensus ranking) on public sepsis single-cell data (Kwok et al.
  Nat Immunol 2023 whole blood incl. neutrophils = primary; Reyes et al. Nat Med 2020; GSE167363 = corroboration).
  Pre-specified axes: HLA class II-CD4, MIF-CD74, CD80/CD86-CD28/CTLA4, CD274-PDCD1, IL10-IL10RA/B,
  S100A8/S100A9-TLR4, ARG1 (neutrophil suppression). Pair scores transferred to the bulk cohorts and tested against
  age and death (same models and meta-analysis as H1).
- Pathway enrichment: per-cohort gene-level age-adjusted death statistics -> random-effects meta-analysis ->
  pre-ranked GSEA (Hallmark, Reactome); leading-edge genes of pathways concordant for death and age form the
  reversal signature.
- Drug repurposing: reversal signature vs cell-type-specific PBMC drug-perturbation signatures (OP3 / NeurIPS 2023,
  GSE279945) and DGIdb / Open Targets for the ligand-receptor targets; positive controls: glucocorticoids
  (dexamethasone / hydrocortisone) for the inflammation axis and IFN-gamma / GM-CSF (sargramostim) for the MHC-II
  axis; LINCS L1000 only as a secondary, caveated check.
```

## Amendment 4

```text
# Amendment 4 - specification of analysis stage 2 (pathways, ligand-receptor axes in bulk)

Written before stage 2 was run. Already seen at this time: the primary results of analysis stage 1 (H1 supported,
H2 not supported). Nothing below was computed before this amendment. The plan follows Amendment 3; this file
fixes the remaining operational details.

Gene-level statistics: in every analysed cohort, genes are standardised within the cohort; death ~ gene + age +
sex (logistic, all cohorts incl. GSE65682 with the binary 28-day label) and gene ~ age + sex (OLS). Genes measured
in all analysed cohorts are meta-analysed (DerSimonian-Laird); the meta z (estimate / SE) ranks the genes.

Pathways: pre-ranked GSEA (Subramanian 2005, weight p = 1) on the death ranking and, separately, the age ranking;
MSigDB Hallmark 2020 and Reactome (Enrichr library; frozen copies saved); sets with 15-500 measured genes;
10,000 gene-set permutations; NES = ES / mean of same-sign null ES; nominal P from the permutation null; FDR =
the NES-based GSEA q-value (Subramanian et al. 2005) within each collection and ranking.

Reversal signature for drug repurposing (as in Amendment 3): leading-edge genes of pathways with death FDR < 0.05
whose age NES has the same sign with P < 0.05 ('death-and-age concordant'), restricted to genes with meta
P(death) < 0.05; UP = higher in non-survivors, DOWN = lower. If either direction has < 20 genes, the
pre-specified signature is reported as not usable and an EXPLORATORY 'death-only' signature (same rule without
the age condition) is reported, clearly labelled as such.

Ligand-receptor axes (pre-specified in Amendment 3): pair score = mean of the within-cohort z of the ligand side
and of the receptor side (multi-gene side = mean z); death ~ pair + age + sex (logistic) and pair ~ age + sex (OLS)
in every cohort, random-effects meta-analysis, BH-FDR over the axes. Secondary: death ~ pair + ML index + age + sex
(is the axis informative beyond the composition index?). Single-cell inference of these axes follows in stage 3.
```

## Amendment 5

```text
# Amendment 5 - specification of analysis stage 3 (drug repurposing); written before any perturbation data
# were downloaded or opened

Signature: S20_reversal_signature_prespecified.tsv (death-and-age concordant leading-edge genes, Amendment 4);
UP = higher in non-survivors and with age, DOWN = lower.

Known limitation stated in advance: part of the UP side reflects emergency erythropoiesis / granulopoiesis
(erythroid and neutrophil genes). Human PBMC cultures contain neither cell type, so these genes cannot be tested
for reversal; only programmes expressed by PBMC cell types are testable. This is handled by the expression gate below.

Data: OP3 / NeurIPS 2023 single-cell perturbation atlas (raw counts; 144 compounds; human PBMCs, 3 donors, 24 h;
DMSO control). Pseudobulk = sum of raw counts per plate x well x cell type (>= 10 cells). Normalisation: TMM;
log-CPM. Model per cell type (as in the OP3 reference pipeline): ~ 0 + compound + plate; contrasts compound - DMSO;
empirical-Bayes variance moderation with a mean-variance trend (limma-trend); moderated t clipped at +-10.
A compound is tested in a cell type if >= 2 pseudobulk samples exist; cytotoxic flag: mean cells of that type per
well < 50% of the DMSO wells on the same plates.

Expression gate (per cell type): a gene is expressed if CPM >= 1 in >= 70% of the DMSO pseudobulks and detected in
>= 5% of the DMSO cells. The reversal analysis is run in a cell type only if >= 20 expressed signature genes remain
in each direction; coverage >= 50% per direction = 'full', otherwise 'partial' (reported). Analysis gene universe =
expressed genes of that cell type.

Primary readout 1 - signature reversal: R = mean t(expressed UP genes) - mean t(expressed DOWN genes) (R < 0 =
reversal). Null: 10,000 random gene sets of the same sizes from the expressed non-signature genes, matched on mean
expression (10 bins); z = (R - null mean) / null SD; two-sided P from the normal approximation (empirical P also
reported); BH-FDR within cell type. Specificity: R among the lowest 10% of all
tested compounds in that cell type (CMap-tau-like).
Primary readout 2 - MHC class II restoration (the strongest protective, age-independent axis of stage 2; HLA class II
- CD4 OR 0.61): log2 fold-change of the summed MHC class II transcripts (frozen pre-registration gene set) per cell
type, same model, t-test, BH-FDR within cell type; a change counts if FDR < 0.05 and |log2FC| >= 0.5.
Secondary readout: S100A8/S100A9 (alarmin) decrease (same rule).
Primary cell type: myeloid cells (carrier of the death-associated myeloid programme); other cell types secondary.

Candidates (not cytotoxic): Tier 1 = reversal (FDR < 0.05, specific) AND MHC-II increase (FDR < 0.05) in myeloid
cells; Tier 2 = only one of the two in myeloid cells, or reversal in >= 2 cell types. Caution list = compounds that
mimic the signature (R > 0, FDR < 0.05, top 10%) or lower myeloid MHC-II (FDR < 0.05).

Positive controls (validity checks, not candidate selection):
(a) technical: glucocorticoids present in OP3 induce the canonical glucocorticoid-receptor targets (TSC22D3, FKBP5,
    ZBTB16, PER1): >= 2 with t > 0 and gene-level FDR < 0.05 in >= 1 cell type;
(b) biological: glucocorticoids lower myeloid MHC class II (known suppression of monocytic HLA-DR): log2FC < 0, P < 0.05;
(c) IFN-gamma and GM-CSF (biologics, not in OP3) are expected to reverse the signature in the Enrichr 'Ligand
    Perturbations from GEO' libraries (adjusted P < 0.05);
(d) target layer: DGIdb lists sargramostim / molgramostim or interferon gamma for the GM-CSF or IFN-gamma receptors.
Target layer: DGIdb interactions for the death-associated ligand-receptor genes, annotated with the desired direction.
LINCS L1000 is not used: < 50% of the signature genes can be landmark genes (lesson from a previous project).
```

## Amendment 6

```text
# Amendment 6 - change of the stage-3 data source; supersedes Amendment 5
# written before any drug-signature result (reversal, MHC class II, candidates) was computed

Amendment 5 specified the OP3 raw single-cell counts (a file of several GB) with our own pseudobulk model. Because
of the download size, the official OP3 differential-expression files (same experiment; limma pseudobulk models of
the OP3 authors) are used instead. No reversal, MHC class II or candidate result had been computed under Amendment 5;
the script checked this before writing this amendment (audit in AMENDMENT_6.json). The specification below
replaces Amendment 5 in full.

## Specification of analysis stage 3 (drug repurposing)

Signature: S20_reversal_signature_prespecified.tsv (death-and-age concordant leading-edge genes, Amendment 4);
UP = higher in non-survivors and with age, DOWN = lower.

Data: official differential-expression results of the OP3 / NeurIPS 2023 PBMC perturbation atlas (Szalata et al.,
NeurIPS 2024; 144 LINCS compounds, human PBMCs from 3 donors, 24 h; limma on pseudobulks of each compound x cell
type vs DMSO, plate-adjusted; train + test files combined = all compound x cell-type profiles). Statistic per gene =
sign(logFC) x -log10(P), clipped at +-4 (the benchmark's own clipping).

Known limitation stated in advance: part of the UP side reflects emergency erythropoiesis / granulopoiesis.
PBMC cultures contain neither erythroid cells nor neutrophils, and OP3 keeps only genes expressed in all PBMC cell
types (edgeR filterByExpr), so these genes are absent. Expression gate: the reversal analysis is run only if
>= 20 signature genes per direction are in the OP3 gene universe; coverage >= 50% per direction = 'full', otherwise
'partial' (reported with the lineage of the missing genes). Cell counts are not part of these files, so a
cytotoxicity filter is not possible (stated as a limitation).

Primary readout 1 - signature reversal per compound x cell type: R = mean statistic(UP genes) - mean statistic(DOWN
genes) (R < 0 = reversal). Null: 10,000 random gene sets of the same sizes drawn from the non-signature OP3 genes
(the same draws for every profile); z = (R - null mean) / null SD; two-sided P from the normal approximation
(empirical P also reported); BH-FDR within cell type. Specificity: R among the lowest 10% of all compounds in that
cell type (CMap-tau-like).
Primary readout 2 - MHC class II restoration (strongest protective, age-independent axis of stage 2; HLA class II
- CD4 OR 0.61): mean statistic of the frozen MHC class II gene set vs random gene sets of the same size (z, normal
P, BH-FDR within cell type); counts as a change if FDR < 0.05 and >= 4 of the MHC class II genes individually change
in the same direction (P < 0.05). Secondary readout: S100A8/S100A9 (alarmin) decrease (FDR < 0.05).
Primary cell type: myeloid cells; other cell types secondary.

Candidates: Tier 1 = reversal (FDR < 0.05, specific) AND MHC class II increase in myeloid cells; Tier 2 = only one of
the two in myeloid cells, or reversal in >= 2 cell types. Caution list = compounds that mimic the signature
(R > 0, FDR < 0.05, top 10%) or lower myeloid MHC class II.

Positive controls (validity checks, not candidate selection):
(a) technical: glucocorticoids present in OP3 induce the canonical glucocorticoid-receptor targets (TSC22D3, FKBP5,
    ZBTB16, PER1): >= 2 with logFC > 0 and gene-level FDR < 0.05 in >= 1 cell type;
(b) biological: glucocorticoids lower myeloid MHC class II (known suppression of monocytic HLA-DR): z < 0, P < 0.05;
(c) IFN-gamma and GM-CSF (biologics, not in OP3) reverse the signature in the Enrichr 'Ligand Perturbations from
    GEO' libraries (adjusted P < 0.05);
(d) target layer: DGIdb lists sargramostim / molgramostim or interferon gamma for the GM-CSF or IFN-gamma receptors.
Target layer: DGIdb interactions for the death-associated ligand-receptor genes, annotated with the desired direction.
LINCS L1000 is not used: < 50% of the signature genes can be landmark genes (lesson from a previous project).
```

## Amendment 7

```text
# Amendment 7 - specification of the closing analyses (stage 05); written before any of them was run

All analyses use the eligible patients and the cohort scores of Amendments 1, 3 and 4 (discovery GSE65682 and the
validation cohorts). Logistic models (death ~ predictor + age per 10 y + sex) in every cohort, DerSimonian-Laird
random-effects pooling, as before.

(A) Severity. In every cohort with a severity score (SOFA / APACHE / SAPS) for >= 70% of the patients: odds ratio of
the ML index, of MHC class II and of the glucocorticoid-response score, adjusted for age, sex and the severity score
(z-scored), next to the unadjusted estimate in the same complete-case sample; pooled over these cohorts.

(B) Incremental discrimination (external): logistic models fitted in the discovery cohort - base (age + sex) and
base + ML index - applied unchanged to every validation cohort; AUC of each and the paired difference with the DeLong
standard error; differences pooled. Within each severity cohort: AUC of (age + sex + severity) vs (+ ML index),
optimism-corrected by 500 bootstrap resamples.

(C) Established signature. The Sweeney mortality score (DEFA4, CD163, RGS1, PER1, HIF1A, SELENOP, IFTAP, CIT up;
LY86, TST, KCNJ2 down; z-score difference; computed if >= 5 UP and >= 2 DOWN genes are measured): association with
death, correlation with the ML index, the ML index adjusted for the Sweeney score, and external incremental AUC of
the ML index over (age + sex + Sweeney score). The SRS classifier is not used (its model is not reproducible from
public information without the original reference data).

(D) In-vivo glucocorticoid-response (GC) score - the link between the drug stage and the patients. Genes: from the
OP3 myeloid profiles of ALL glucocorticoids in OP3 (stage 3), genes with statistic >= +3 (UP) or <= -3 (DOWN) in every
glucocorticoid; restricted to genes measured in all analysed cohorts; genes of the death-age signature (S20) and of
the frozen composition / MHC class II gene sets excluded (no overlap with what is tested); top 50 per direction by the
mean statistic; >= 10 per direction required. Score = mean z(UP) - mean z(DOWN), standardised within cohort.
Primary tests (all cohorts pooled; the gene set was derived outside the patient data, so every cohort is an
independent test; Holm over the two): (D1) death ~ GC score + age + sex; (D2) GC score ~ age + sex.
Secondary: death ~ GC score adjusted for the ML index and MHC class II; Spearman correlations of the GC score with
MHC class II and with the ML index (Fisher-z pooled). Interpretation: in vivo the score can also reflect endogenous
cortisol and leukocyte composition; this is stated, and the adjusted model addresses composition.

Not done (reasons stated in the paper): single-cell validation (the only public PBMC data set with survivors and
non-survivors, GSE167363, has 2 non-survivors); healthy-ageing anchor (needs multi-GB single-cell downloads).
```

## Amendment 8

```text
# Amendment 8 - external single-cell validation (stage 06); written before the single-cell data were downloaded

Data: GSE167363 (Qiu et al., J Leukoc Biol 2021): 10x PBMC scRNA-seq of 5 patients with Gram-negative sepsis
(3 survivors, 2 non-survivors) at sepsis recognition (T0) and 6 h later (T6), and 2 non-septic donors; processed
count matrices as deposited. Patient-level inference only; cells are never used as independent replicates for a
group comparison. With 2 non-survivors and 2 controls, outcome and case-control contrasts are descriptive.

Processing (pre-specified): cells with >= 200 genes and <= 20% mitochondrial reads; genes in >= 3 cells; library-
size normalisation (10,000) and log1p; 2,000 highly variable genes (per-sample batches); 30 PCs; Harmony integration
by sample if harmonypy is installed, otherwise none (reported); kNN graph, UMAP, Leiden clusters; clusters labelled
by the highest mean z-scored marker score among 12 frozen marker sets (listed in the code, class A5).
Cell scores (scanpy score_genes): glucocorticoid-response (GC) score = score(S27 UP) - score(S27 DOWN), the gene
set frozen in Amendment 7 / S27; MHC class II score = frozen pre-registration set; death-age signature score =
score(S20 UP) - score(S20 DOWN).

Primary test (sc1) - coupling in patients' monocytes: within the CD14 monocytes of each septic patient at T0
(>= 100 cells), partial Spearman correlation of the GC score with the MHC class II score, adjusted for log total
counts; Fisher-z random-effects pooling over patients. Supported if the pooled rho < 0 with P < 0.05 AND >= 4 of 5
patients have rho < 0. Sensitivity: T6 samples; healthy donors reported alongside.
Descriptive (sc2) - localisation: donor-level mean GC score per cell type in septic patients (T0) minus healthy
donors; the bulk/OP3 data predict the largest increase in monocytes.
Descriptive (sc3): death-age signature score per cell type (the erythroid / heme part is expected in erythroid
cells); donor-level monocyte GC and MHC class II scores by outcome group (no test).
```
