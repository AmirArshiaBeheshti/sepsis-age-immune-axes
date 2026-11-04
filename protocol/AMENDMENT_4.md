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
