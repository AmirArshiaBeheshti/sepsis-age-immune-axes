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
