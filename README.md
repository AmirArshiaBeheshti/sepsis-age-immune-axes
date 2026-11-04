# Age and blood immune dysregulation are separable mortality risk axes in sepsis: a multi-cohort transcriptomic study

Analysis code, protocol and supplementary results for the article of the same title (*Heliyon*, [year; DOI to be added]).

## Summary

We analyzed seven public blood transcriptome cohorts of adults with sepsis (1375 patients, 315 deaths). The analysis followed a time-stamped, hash-locked protocol that was not deposited in a public registry. All amendments are provided here.

The two primary hypotheses were:

- **H1:** a myelopoiesis–lymphopenia (ML) index (immature-neutrophil minus T-cell gene-set score) is associated with short-term mortality.
- **H2:** blood immune scores mediate the effect of age on death.

H1 was supported (validation pooled OR 1.41 per SD, 95% CI 1.18–1.68). H2 was not supported (5% mediated, −11% to 31%).

Secondary analyses covered:

- pathways and ligand–receptor axes;
- severity adjustment and comparison with an established 11-gene mortality score;
- drug-perturbation connectivity in human PBMCs (Open Problems / NeurIPS 2023 data set);
- single-cell validation (GSE167363);
- additional post hoc sensitivity analyses.

All data are public. This repository does not redistribute raw data; the scripts download them from the original repositories.

## Repository structure

```
code/                         analysis scripts 01-08 (run in this order)
tools/finalize_repository.py      one cell: copies and verifies the protocol, records versions, writes the Methods 2.14 text
tools/verify_protocol_hashes.py   checks the registered protocol files against their SHA-256 records
tools/record_environment.py       writes the software and hardware versions used for a run
protocol/                     registered protocol and amendments (.md) with their .json records and REGISTRATION_NOTE.md
logs/                         run logs of every stage (they document the registration times)
docs/protocol_and_amendments.md   the same texts, reproduced from the code for reading
results/supplementary/        Supplementary_Data_S0-S45.xlsx (all result tables of the article)
results/figures/              Figures S1-S2 (Figures 1-7 are written by the scripts)
requirements.txt, environment.yml, CITATION.cff, LICENSE
```

## Data sources

| Data set | Role | Source |
| --- | --- | --- |
| GSE65682 (MARS) | Discovery cohort | GEO |
| E-MTAB-4421, E-MTAB-4451, E-MTAB-5273 (GAinS) | Validation cohorts | ArrayExpress / BioStudies |
| GSE95233, GSE54514, GSE185263 | Validation cohorts | GEO |
| E-MTAB-5274, GSE63042 | Screened, not analyzed (event threshold; no count table) | ArrayExpress, GEO |
| GSE33828 (Rotterdam Study) | Healthy reference for the withdrawn clock hypotheses; not used in reported analyses | GEO |
| GSE167363 | Single-cell PBMC validation | GEO |
| Open Problems single-cell perturbation data set (NeurIPS 2023), differential-expression files `2023-09-12_de_by_cell_type_train.h5ad` and `..._test.h5ad` | Drug perturbation | openproblems-bio public S3 bucket (URLs in the header of `code/04_drug_repurposing.py`) |
| MSigDB Hallmark 2020 and Reactome libraries (Enrichr) | Gene sets | Enrichr |
| HGNC complete set | Gene symbols | HGNC |
| DGIdb 5.0 (GraphQL API) | Drug–gene interactions | dgidb.org |

All public data were retrieved on 2 November 2026 (run logs; `environment_versions.txt`).

## Requirements

Python 3.10 or later. Install the packages with:

```
pip install -r requirements.txt
```

or create the conda environment with `conda env create -f environment.yml`.

- Scripts 01–03, 05 and 08 need only the scientific Python stack.
- Script 04 additionally needs `anndata` and `h5py`.
- Script 06 needs `scanpy`, `umap-learn`, `python-igraph`, `leidenalg` and `harmonypy`.

The published results were produced with Python 3.12.13 on Linux, with NumPy 2.0.2, pandas 2.3.3, SciPy 1.16.3, statsmodels 0.14.6, scikit-learn 1.9.0, matplotlib 3.10.9, Scanpy 1.12.4, anndata 0.12.19, harmonypy 2.0.0, umap-learn 0.5.12, igraph 1.0.0, leidenalg 0.12.0 and requests 2.34.2. The full list is in `environment_versions.txt`. You can record your own environment with `python tools/record_environment.py`.

## How to run

Each script is self-contained. Run it either by pasting the whole file into one Jupyter cell, or from the repository root with `python code/<script>.py`.

All output goes to `./SEPSIS_immuneage`, or to the folder named in the environment variable `SEPSIS_PROJECT_DIR`. Run the scripts in order; each one checks that the previous stages finished and that their quality gates passed.

| Script | Content | Main outputs | Approximate run time |
| --- | --- | --- | --- |
| `01_data_preparation.py` | Writes protocol v1, downloads cohorts, metadata audit, harmonization, gene mapping | S0 (stage 1), S1–S4, Fig. 1 | Depends on download speed |
| `02_analysis_H1_H2.py` | Amendment 3; immune scores; H1 meta-analysis; H2 mediation (1000 bootstraps per cohort) | S10–S17, Figs. 2–3 | 10–30 min |
| `03_analysis_pathways_LR.py` | Amendment 4; gene-level meta-analysis; pre-ranked GSEA; death–age signature; ligand–receptor axes | S18–S21b, Fig. 4 | 5–15 min |
| `04_drug_repurposing.py` | Amendments 5/6; OP3 reversal and MHC class II scoring; positive controls; DGIdb | S22–S26, Fig. 5 | Minutes |
| `05_closing_analyses.py` | Amendment 7; severity adjustment; external AUC; 11-gene score; glucocorticoid-response score | S27–S32, Fig. 6 | Minutes |
| `06_single_cell_validation.py` | Amendment 8; GSE167363 QC, Harmony, UMAP, Leiden; monocyte coupling test | S33–S36, Fig. 7 | 5–15 min |
| `07_figure_polish.py` | Redraws panels of Figs. 4, 5 and 7 from saved tables only (no result changes) | Figs. 4, 5, 7 | < 1 min |
| `08_revision_analyses.py` | Additional post hoc sensitivity analyses | S37–S45 | About 2 min |

**Manual steps.**

- **Script 04:** if the automatic download of the two OP3 files fails, download them from the URLs in the script header and place them in `SEPSIS_immuneage/data/raw/op3/`.
- **Script 06:** the GSE167363 download (321 MB) resumes if interrupted. Run the cell again.
- **Script 08:** reads only saved tables. If `SEPSIS_immuneage/data/raw/GSE54514/GSE54514_series_matrix.txt.gz` is present (it is downloaded by script 01), it also writes S45, the comparison with the measured neutrophil proportion.

## Protocol and amendments

On first run, each script writes its protocol text to `SEPSIS_immuneage/preregistration/` as a `.md` file. Each text is written before the analysis it governs, with a `.json` record holding the UTC time stamp and SHA-256 hash. If a registered text is later changed, the scripts stop.

To publish the records, run `tools/finalize_repository.py` once, in the folder that contains both `SEPSIS_immuneage` and this repository. It copies the files into `protocol/`, verifies every hash, records the software versions and writes the text for Methods 2.14. To verify the records again later:

```
python tools/verify_protocol_hashes.py protocol
```

Table 1 of the article lists which analyses are confirmatory (H1, H2), secondary, exploratory or post hoc. In brief:

- Amendment 3 withdrew the original clock-based hypotheses before any outcome analysis.
- Amendment 4 and the later amendments were written after the primary results were known.
- Script 08 contains post hoc analyses.

## Reproducibility notes

- The random seed is 0 throughout.
- Script 08 first checks that it reproduces the per-cohort H1 estimates exactly.
- Point estimates are deterministic.
- The number of bootstrap resamples that converge in the smallest cohorts can differ by a few between `statsmodels` versions, for example 872 versus 873 of 1000 in GSE54514. Confidence limits can change in the last digit as a result.
- Every script writes its quality gates (S0 tables) and a log to `SEPSIS_immuneage/logs/`.

## Supplementary tables

`results/supplementary/Supplementary_Data_S0-S45.xlsx` contains every result table of the article, one sheet per table. The first sheet is an index. Tables S5–S9 were reserved for the withdrawn clock analyses and do not exist.

## Citation

Please cite the article (see `CITATION.cff`) and this repository (https://github.com/AmirArshiaBeheshti/sepsis-age-immune-axes).

## License

The code is released under the MIT License (see `LICENSE`). The public data sets remain under the terms of their original repositories.

## Use of generative AI

Generative AI (Claude, Anthropic) assisted with writing analysis code. All code was reviewed and run by the authors, who take full responsibility for it.
