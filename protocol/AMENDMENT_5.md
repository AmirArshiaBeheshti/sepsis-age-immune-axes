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
