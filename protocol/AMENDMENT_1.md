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
