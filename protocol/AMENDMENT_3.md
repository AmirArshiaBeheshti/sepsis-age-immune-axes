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
