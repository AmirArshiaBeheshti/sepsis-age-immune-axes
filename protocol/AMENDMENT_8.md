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
