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
