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
