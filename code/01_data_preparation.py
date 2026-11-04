# =====================================================================================
# SEPSIS IMMUNE-AGE PROJECT - STAGE 1: pre-registration, data acquisition, metadata audit and harmonisation
#
# Study: "Decomposing the 28-day mortality signal of the sepsis blood transcriptome into cell-intrinsic immune
#         ageing, leukocyte composition and inflammation - single-cell-informed clocks, multi-cohort validation
#         and drug repurposing"
#   (working title when this stage was written; published title: "Age and blood immune dysregulation are
#   separable mortality risk axes in sepsis: a multi-cohort transcriptomic study")
#
# HOW TO RUN: paste the whole file into ONE Jupyter cell and run it. Everything is written to ./SEPSIS_immuneage
#   (or to the folder in the environment variable SEPSIS_PROJECT_DIR). Downloads are cached; re-running only
#   repeats what is missing. Packages: numpy, pandas, matplotlib, requests (no extra installs).
#
# Steps of Stage 1 (no outcome association is computed in Stage 1):
#   0  pre-registration: hypotheses, models, gates and gene sets are written ONCE with a UTC time stamp and a
#      SHA-256 hash, BEFORE any outcome data are parsed. A later change of the text stops the run (amendments are
#      written as separate files).
#   1  downloads: GEO series matrices / supplementary count tables / platform tables, ArrayExpress (BioStudies)
#      SDRF + processed data, HGNC gene table, MSigDB Hallmark (via Enrichr)
#   2  metadata audit (every field of every cohort, S1) and harmonised clinical table (S2/S3) with the
#      pre-specified parsing and inclusion rules
#   3  expression: probe/Ensembl -> HGNC symbol, log2 scale, one sample per patient; gene-set coverage (S4)
#   4  quality gates (S0) and Fig 1 (cohort overview)
#
# OUTPUTS TO CHECK: S0_stage1_quality_gates.tsv, S1_metadata_field_audit.tsv, S2_cohort_summary.tsv and the log.
#   The field audit is checked before Stage 2; parsing corrections are allowed only as a written amendment and
#   only before any outcome analysis (Stage 3).
# =====================================================================================
from __future__ import annotations

import datetime
import gzip
import hashlib
import io
import json
import logging
import os
import re
import sys
import time
import zipfile
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests


# =====================================================================================
# CONFIGURATION
# =====================================================================================
class C:
    PROJECT_DIR = Path(os.environ.get("SEPSIS_PROJECT_DIR", Path.cwd() / "SEPSIS_immuneage")).resolve()
    RAW = PROJECT_DIR / "data" / "raw"
    PROC = PROJECT_DIR / "data" / "processed"
    TAB = PROJECT_DIR / "results" / "tables"
    FIG = PROJECT_DIR / "results" / "figures"
    PREREG = PROJECT_DIR / "preregistration"
    LOGS = PROJECT_DIR / "logs"
    SEED = 0

    GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo"
    GEO_ACC = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
    BIOSTUDIES_INFO = "https://www.ebi.ac.uk/biostudies/api/v1/studies/{acc}/info"
    HGNC_URLS = ["https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
                 "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt"]
    ENRICHR_LIB = ("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020")

    # role: discovery | validation | healthy_reference
    # outcome: definition of the death label in the source publication (harmonised as 'death' 0/1)
    COHORTS = {
        "GSE65682": dict(source="geo", role="discovery", data="array", tissue="whole blood",
                         outcome="28-day mortality", note="MARS; ICU admission; includes non-infectious ICU"),
        "E-MTAB-4421": dict(source="ae", all_sepsis=True, role="validation", data="array", tissue="leukocytes",
                            outcome="28-day mortality", array_gpl="GPL10558", note="GAinS CAP discovery"),
        "E-MTAB-4451": dict(source="ae", all_sepsis=True, role="validation", data="array", tissue="leukocytes",
                            outcome="28-day mortality", array_gpl="GPL10558", note="GAinS CAP validation"),
        "E-MTAB-5273": dict(source="ae", all_sepsis=True, role="validation", data="array", tissue="leukocytes",
                            outcome="28-day mortality", array_gpl="GPL10558", note="GAinS faecal peritonitis"),
        "E-MTAB-5274": dict(source="ae", all_sepsis=True, role="validation", data="array", tissue="leukocytes",
                            outcome="28-day mortality", array_gpl="GPL10558", note="GAinS faecal peritonitis"),
        "GSE63042": dict(source="geo", role="validation", data="rnaseq", tissue="whole blood",
                         outcome="28-day mortality", note="CAPSOD; ED presentation"),
        "GSE95233": dict(source="geo", role="validation", data="array", tissue="whole blood",
                         outcome="hospital survival", note="septic shock; ICU day 1"),
        "GSE54514": dict(source="geo", role="validation", data="array", tissue="whole blood",
                         outcome="28-day mortality", note="daily samples; first sample used"),
        "GSE185263": dict(source="geo", role="validation", data="rnaseq", tissue="whole blood",
                          outcome="in-hospital mortality", note="already analysed by TraMA (Aging 2025)"),
        "GSE33828": dict(source="geo", role="healthy_reference", data="array", tissue="whole blood",
                         outcome="none", note="Rotterdam Study healthy reference (age in supplementary file)"),
    }
    ILLUMINA_GPLS = ("GPL10558", "GPL6947", "GPL6883", "GPL6884")   # ILMN_ probe ids are shared

    # inclusion / gates (pre-specified)
    MIN_AGE = 18
    DISCOVERY_MIN_N, DISCOVERY_MIN_EVENTS = 300, 60
    VALIDATION_MIN_N, VALIDATION_MIN_EVENTS = 30, 10
    VALIDATION_MIN_TOTAL_EVENTS, VALIDATION_MIN_COHORTS = 100, 2
    HEALTHY_MIN_N = 500
    MIN_GENES = 10000
    MIN_SET_COVERAGE = 0.60

    # gene sets (frozen in the pre-registration). Composition surrogates are marker sets from published blood
    # single-cell atlases; MS1 = sepsis-expanded CD14+ monocyte state (Reyes et al., Nat Med 2020).
    GENE_SETS = {
        "immature_neutrophil": ["CEACAM8", "MMP8", "LCN2", "LTF", "CAMP", "OLFM4", "DEFA4", "BPI", "CRISP3", "CD24"],
        "mature_neutrophil": ["FCGR3B", "CXCR2", "CSF3R", "FPR1", "MME", "CMTM2", "PROK2", "FFAR2"],
        "ms1_monocyte": ["RETN", "ALOX5AP", "IL1R2", "CD163", "SLC39A8"],
        "classical_monocyte": ["CD14", "LYZ", "FCN1", "VCAN", "CSF1R", "CD68"],
        "t_cell": ["CD3E", "CD3D", "CD3G", "CD247", "CD2", "LCK", "CD28", "IL7R", "CCR7", "TRAT1"],
        "b_cell": ["MS4A1", "CD79A", "CD79B", "CD19", "CD22", "BANK1", "FCRLA", "TCL1A"],
        "nk_cell": ["GNLY", "NKG7", "KLRD1", "PRF1", "KLRF1", "GZMB", "SH2D1B"],
        "mhc_class_ii": ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DMA", "HLA-DMB", "CD74"],
    }
    HALLMARK_TERM = "Inflammatory Response"


# =====================================================================================
# INFRASTRUCTURE: logging, gates, style
# =====================================================================================
for _d in (C.RAW, C.PROC, C.TAB, C.FIG, C.PREREG, C.LOGS):
    _d.mkdir(parents=True, exist_ok=True)


def get_logger(name: str) -> logging.Logger:
    lg = logging.getLogger(name)
    lg.setLevel(logging.INFO)
    lg.propagate = False
    if not lg.handlers:
        fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s", "%Y-%m-%d %H:%M:%S")
        sh = logging.StreamHandler(sys.stdout)
        sh.setFormatter(fmt)
        fh = logging.FileHandler(C.LOGS / f"{name}.log")
        fh.setFormatter(fmt)
        lg.addHandler(sh)
        lg.addHandler(fh)
    return lg


LOG = get_logger("sepsis_stage1")
GATES_FILE = C.TAB / "S0_stage1_quality_gates.tsv"
_GATES: list = []


def record_gates(step: str, rows: list) -> None:
    for r in rows:
        r = {"step": step, **r}
        r["status"] = "PASS" if r["passed"] else ("FAIL" if r["severity"] == "critical" else "WARN")
        _GATES.append(r)
        LOG.info(f"QUALITY GATE [{r['status']}] {r['gate']}: {r['value']} (criterion: {r['criterion']})")
    pd.DataFrame(_GATES)[["step", "gate", "value", "criterion", "passed", "severity", "status"]].to_csv(
        GATES_FILE, sep="\t", index=False)


def gate_summary() -> str:
    df = pd.DataFrame(_GATES)
    lines = ["", "=" * 90, "STAGE 1 QUALITY GATES", "=" * 90]
    lines += [f"[{r.status}] {r.step:<14} {r.gate}: {r.value}  (criterion: {r.criterion})" for r in df.itertuples()]
    lines.append("-" * 90)
    fails = df[df["status"] == "FAIL"]
    lines.append("ALL CRITICAL GATES PASSED - Stage 2 can start after the field audit has been checked."
                 if fails.empty else f"{len(fails)} CRITICAL GATE(S) FAILED - see above; do not start Stage 2.")
    lines.append("=" * 90)
    return "\n".join(lines)


MM = 1 / 25.4
FIG_W = 183 * MM
OKABE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]


def set_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 6, "axes.titlesize": 6.5, "axes.labelsize": 6, "xtick.labelsize": 5.5, "ytick.labelsize": 5.5,
        "legend.fontsize": 5.5, "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 2, "ytick.major.size": 2, "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "savefig.dpi": 600, "figure.dpi": 150})


def panel_label(ax, letter: str, x: float = -0.12, y: float = 1.04) -> None:
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="left")


def save_figure(fig, stem: str) -> None:
    for ext in ("pdf", "png", "tiff"):
        kw = {"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tiff" else {}
        fig.savefig(C.FIG / f"{stem}.{ext}", bbox_inches="tight", dpi=600, **kw)
    plt.close(fig)
    LOG.info(f"Saved figure {C.FIG / stem}.pdf/.png/.tiff")


def utc_now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


# =====================================================================================
# STEP 0 - PRE-REGISTRATION
# =====================================================================================
def preregistration_text() -> str:
    gs = "\n".join(f"  - {k}: {', '.join(v)}" for k, v in C.GENE_SETS.items())
    coh = "\n".join(f"  - {k}: role={v['role']}; {v['data']}; {v['tissue']}; outcome={v['outcome']}; {v['note']}"
                    for k, v in C.COHORTS.items())
    return f"""# Pre-registration v1 - Sepsis immune-age decomposition study

## Rationale and prior art (checked before registration)
TraMA (Klopack et al., Aging 2025) applied a bulk transcriptomic mortality-age clock to one sepsis cohort
(GSE185263): sepsis patients appeared ~10 years older, and the association with mortality weakened to P = 0.09
after adjustment for cell type. Epigenetic PhenoAge acceleration was associated with sepsis non-survival (Crit Care
Explor 2024). Single-cell immune ageing clocks (Nature Aging 2025; IDEAL-Age, Genome Biol 2026) were not applied to
sepsis or to mortality. No study has separated cell-intrinsic immune ageing from leukocyte composition and
inflammation in multiple independent sepsis cohorts.

## Cohorts
{coh}

## Inclusion rules (applied identically to every cohort)
1. Adult patients (age >= {C.MIN_AGE}) with sepsis / septic shock / sepsis due to CAP or faecal peritonitis as defined
   by the source study; healthy and non-infectious ICU samples are kept as separate groups.
2. One sample per patient: the earliest available time point (ties: first sample in the series).
3. Known age and known death label (0/1). Outcome definitions differ (28-day / hospital / in-hospital) and are kept
   per cohort; cohort-specific estimates are pooled by random-effects meta-analysis.
4. A validation cohort is analysed if it has >= {C.VALIDATION_MIN_N} eligible patients and >= {C.VALIDATION_MIN_EVENTS} deaths.
5. Expression: probes -> HGNC approved symbols (probes mapping to >1 gene removed; highest-mean probe per gene);
   Ensembl ids -> HGNC; arrays on log2 scale; RNA-seq counts -> log2(CPM + 1). All scores use genes standardised
   within each cohort among the eligible sepsis patients.

## Pre-specified scores
Gene sets (frozen; a score = mean within-cohort z of available genes; computed only if >= {int(C.MIN_SET_COVERAGE * 100)}% of the
genes are measured):
{gs}
  - inflammation: MSigDB Hallmark '{C.HALLMARK_TERM}' (Enrichr library MSigDB_Hallmark_2020, frozen copy saved).
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
"""


def write_preregistration() -> None:
    p = C.PREREG / "PREREG_v1.md"
    meta = C.PREREG / "PREREG_v1.json"
    text = preregistration_text()
    h = hashlib.sha256(text.encode()).hexdigest()
    if p.exists():
        old = json.loads(meta.read_text())
        if hashlib.sha256(p.read_text().encode()).hexdigest() != old["sha256"]:
            raise SystemExit("PREREG_v1.md was edited after registration - restore it or write an amendment file.")
        if h != old["sha256"]:
            raise SystemExit("The pre-registration text in the code differs from the registered version "
                             f"({old['registered_utc']}). Keep the code text unchanged or write an amendment.")
        LOG.info(f"Pre-registration present: registered {old['registered_utc']}, SHA-256 {old['sha256'][:16]}... (unchanged)")
        return
    p.write_text(text)
    meta.write_text(json.dumps({"registered_utc": utc_now(), "sha256": h, "file": p.name}, indent=1))
    LOG.info(f"Pre-registration written {p} (SHA-256 {h[:16]}...). Outcome data have not been parsed yet.")


AMENDMENT_1 = """# Amendment 1 to Pre-registration v1 - cohort-specific metadata rules

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
"""


def write_amendment() -> None:
    p, meta = C.PREREG / "AMENDMENT_1.md", C.PREREG / "AMENDMENT_1.json"
    h = hashlib.sha256(AMENDMENT_1.encode()).hexdigest()
    if p.exists():
        old = json.loads(meta.read_text())
        if hashlib.sha256(p.read_text().encode()).hexdigest() != old["sha256"] or h != old["sha256"]:
            raise SystemExit("AMENDMENT_1.md differs from its registered version - do not edit it; add a new amendment.")
        LOG.info(f"Amendment 1 present (registered {old['registered_utc']}, unchanged)")
        return
    p.write_text(AMENDMENT_1)
    meta.write_text(json.dumps({"registered_utc": utc_now(), "sha256": h, "file": p.name}, indent=1))
    LOG.info(f"Amendment 1 written {p} (SHA-256 {h[:16]}...) - before any outcome analysis")


# =====================================================================================
# STEP 1 - DOWNLOADS
# =====================================================================================
def file_is_valid(p: Path) -> bool:
    if not p.exists() or p.stat().st_size == 0:
        return False
    try:
        if p.suffix == ".gz":
            with gzip.open(p, "rb") as fh:
                while fh.read(1 << 20):
                    pass
        elif p.suffix == ".zip":
            with zipfile.ZipFile(p) as z:
                return z.testzip() is None
        return True
    except (OSError, EOFError, zipfile.BadZipFile):
        return False


HTTP = requests.Session()
HTTP.headers.update({"User-Agent": "Mozilla/5.0 (academic research script; sepsis immune-age study)"})
_LAST_NCBI = [0.0]


def http_get(url: str, stream: bool = False, timeout: int = 300, tries: int = 5):
    """GET with a polite pause for NCBI and back-off on 403/429/5xx (NCBI answers 403 when it rate-limits)."""
    last = None
    for k in range(tries):
        if "ncbi.nlm.nih.gov" in url:
            wait = 1.0 - (time.time() - _LAST_NCBI[0])
            if wait > 0:
                time.sleep(wait)
            _LAST_NCBI[0] = time.time()
        try:
            r = HTTP.get(url, stream=stream, timeout=timeout)
            if r.status_code == 404:
                r.raise_for_status()
            if r.status_code in (403, 429) or r.status_code >= 500:
                last = requests.HTTPError(f"{r.status_code} for {url}")
                r.close()
                time.sleep(15 * (k + 1))
                continue
            r.raise_for_status()
            return r
        except requests.HTTPError:
            raise
        except requests.RequestException as exc:
            last = exc
            time.sleep(10 * (k + 1))
    raise requests.HTTPError(f"failed after {tries} attempts: {last}")


def fetch(url: str, dest: Path, retries: int = 3) -> Path:
    if file_is_valid(dest):
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    last = None
    for k in range(retries):
        try:
            with http_get(url, stream=True) as r:
                with open(tmp, "wb") as fh:
                    for chunk in r.iter_content(1 << 20):
                        fh.write(chunk)
            tmp.replace(dest)
            if not file_is_valid(dest):
                raise OSError("downloaded file is corrupt")
            LOG.info(f"Downloaded {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)")
            return dest
        except (requests.RequestException, OSError) as exc:
            last = exc
            if tmp.exists():
                tmp.unlink()
            if "404" in str(exc):          # the file does not exist: no point in retrying
                break
            time.sleep(5 * (k + 1))
    raise RuntimeError(f"download failed: {url} ({last})")


def list_dir(url: str) -> list:
    """File names in an HTTP directory index (NCBI / EBI FTP over https)."""
    r = http_get(url, timeout=120, tries=3)
    names = re.findall(r'href="([^"?/][^"]*)"', r.text)
    return sorted({n.rstrip("/") for n in names if not n.startswith(("http", "#"))})


def geo_page_files(acc: str) -> dict:
    """Fallback when the FTP index is refused: file names linked from the GEO accession page."""
    from urllib.parse import unquote
    html = http_get(f"{C.GEO_ACC}?acc={acc}", timeout=120).text
    names = {unquote(x) for x in re.findall(r"file=([^\"&]+)", html)}
    names |= {unquote(x) for x in re.findall(r"/suppl/([^\"/]+)\"", html)}
    return {"suppl": sorted(n for n in names if n.startswith(acc) and "series_matrix" not in n)}


def geo_stub(acc: str) -> str:
    num = re.sub(r"\D", "", acc)
    prefix = re.sub(r"\d", "", acc)
    return f"{prefix}{num[:-3]}nnn" if len(num) > 3 else f"{prefix}nnn"


def download_geo(acc: str) -> dict:
    base = f"{C.GEO_FTP}/series/{geo_stub(acc)}/{acc}"
    out = C.RAW / acc
    files = {"matrix": [], "suppl": []}
    try:                                      # the usual single-platform file name: no directory listing needed
        files["matrix"].append(fetch(f"{base}/matrix/{acc}_series_matrix.txt.gz", out / f"{acc}_series_matrix.txt.gz"))
    except RuntimeError:
        try:
            for name in [n for n in list_dir(f"{base}/matrix/") if n.endswith("_series_matrix.txt.gz")]:
                files["matrix"].append(fetch(f"{base}/matrix/{name}", out / name))
        except (requests.RequestException, RuntimeError) as exc:     # supplementary files are still fetched
            LOG.warning(f"{acc}: series matrix not available ({exc})")
    try:
        sup = list_dir(f"{base}/suppl/")
    except requests.RequestException:
        try:
            sup = geo_page_files(acc)["suppl"]
        except requests.RequestException:
            sup = []
    sup = [n for n in sup if not n.endswith(("_RAW.tar", ".tar", ".CEL.gz", ".idat.gz"))]
    for name in sup:
        if re.search(r"count|expr|fpkm|tpm|cpm|norm|matrix|age|gender|sex|clinical|pheno|sample|meta", name, re.I):
            files["suppl"].append(fetch(f"{base}/suppl/{name}", out / "suppl" / name))
    return files


def download_gpl(gpl: str) -> Path:
    dest = C.RAW / "platforms" / f"{gpl}.txt"
    if file_is_valid(dest):
        return dest
    try:                                  # complete annotation table (GEO 'full table' view)
        return fetch(f"{C.GEO_ACC}?acc={gpl}&targ=self&form=text&view=data", dest)
    except RuntimeError:
        annot = C.RAW / "platforms" / f"{gpl}.annot.gz"
        return fetch(f"{C.GEO_FTP}/platforms/{geo_stub(gpl)}/{gpl}/annot/{gpl}.annot.gz", annot)


def download_ae(acc: str) -> dict:
    out = C.RAW / acc
    j = http_get(C.BIOSTUDIES_INFO.format(acc=acc), timeout=120).json()
    link = j.get("httpLink") or j.get("ftpLink", "").replace("ftp://", "https://")
    if not link:
        raise RuntimeError(f"no file link for {acc} in BioStudies info: {list(j)}")
    names = list_dir(link.rstrip("/") + "/Files/")
    files = {"sdrf": [], "processed": []}
    for n in names:
        if n.endswith("sdrf.txt"):
            files["sdrf"].append(fetch(f"{link}/Files/{n}", out / n))
        elif re.search(r"processed|normali[sz]ed|matrix|expression", n, re.I) and not n.endswith((".idf.txt",)):
            files["processed"].append(fetch(f"{link}/Files/{n}", out / n))
    return files


def download_hgnc() -> Path:
    dest = C.RAW / "reference" / "hgnc_complete_set.txt"
    for u in C.HGNC_URLS:
        try:
            return fetch(u, dest)
        except RuntimeError as exc:
            LOG.warning(f"HGNC source failed ({exc}); trying the next one")
    raise SystemExit("HGNC gene table could not be downloaded.")


def download_hallmark() -> Path:
    dest = C.RAW / "reference" / "MSigDB_Hallmark_2020.txt"
    return fetch(C.ENRICHR_LIB, dest)


def manual_instructions(acc: str, cfg: dict) -> str:
    if cfg["source"] == "geo":
        return (f"download {acc}_series_matrix.txt.gz (and, for RNA-seq, the count table from 'Supplementary file') "
                f"from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc} into {C.RAW / acc}")
    return (f"download the SDRF and the processed data files from https://www.ebi.ac.uk/biostudies/arrayexpress/"
            f"studies/{acc} into {C.RAW / acc}")


def step1_main() -> dict:
    status = {}
    download_hgnc()
    try:
        download_hallmark()
    except RuntimeError as exc:
        LOG.warning(f"MSigDB Hallmark (Enrichr) not downloaded ({exc}); inflammation score will be unavailable")
    for acc, cfg in C.COHORTS.items():
        try:
            status[acc] = download_geo(acc) if cfg["source"] == "geo" else download_ae(acc)
            LOG.info(f"{acc}: files ready")
        except Exception as exc:          # one unreachable cohort must not stop the others
            local = sorted((C.RAW / acc).glob("*")) if (C.RAW / acc).exists() else []
            if local:
                LOG.warning(f"{acc}: remote listing failed ({exc}); using {len(local)} local file(s)")
                status[acc] = {"local": local}
            else:
                LOG.warning(f"{acc}: NOT downloaded ({exc}). Manual: {manual_instructions(acc, cfg)}")
                status[acc] = None
    return status


# =====================================================================================
# PARSERS
# =====================================================================================
def read_series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (sample metadata with one column per field, expression table probes x GSM; may be empty)."""
    meta, begin = {}, None
    with gzip.open(path, "rt", errors="replace") as fh:
        for i, line in enumerate(fh):
            if line.startswith("!series_matrix_table_begin"):
                begin = i
                break
            if not line.startswith("!Sample_"):
                continue
            key, *vals = line.rstrip("\n").split("\t")
            vals = [v.strip().strip('"') for v in vals]
            key = key[len("!Sample_"):]
            if key.startswith("characteristics"):
                for j, v in enumerate(vals):
                    if ":" in v:
                        k, val = v.split(":", 1)
                        meta.setdefault(f"ch:{k.strip().lower()}", [""] * len(vals))[j] = val.strip()
                    elif v:
                        meta.setdefault("ch:unlabelled", [""] * len(vals))[j] = v
            else:
                if key in meta:                                   # repeated field (e.g. several descriptions)
                    meta[key] = [f"{a} | {b}" if b else a for a, b in zip(meta[key], vals)]
                else:
                    meta[key] = vals
    md = pd.DataFrame(meta)
    md.index = md["geo_accession"]
    expr = pd.DataFrame()
    if begin is not None:
        expr = pd.read_csv(path, sep="\t", skiprows=begin + 1, index_col=0, low_memory=False)
        expr = expr[~expr.index.astype(str).str.startswith("!")]
        expr = expr.apply(pd.to_numeric, errors="coerce").dropna(how="all")
        expr.index = expr.index.astype(str)
    return md, expr


def read_platform(path: Path, id_col: str | None = None) -> pd.Series:
    """Probe id (or another id column, e.g. Illumina Array_Address_Id) -> gene symbol (multi-gene probes removed)."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", errors="replace") as fh:
        lines = fh.read().splitlines()
    start = next((i for i, l in enumerate(lines) if l.startswith("!platform_table_begin")), None)
    body = lines[start + 1:] if start is not None else [l for l in lines if not l.startswith(("#", "!", "^"))]
    body = [l for l in body if not l.startswith("!platform_table_end")]
    tab = pd.read_csv(io.StringIO("\n".join(body)), sep="\t", dtype=str, low_memory=False)
    col = next((c for c in ["Gene Symbol", "Gene symbol", "Symbol", "ILMN_Gene", "GENE_SYMBOL", "gene_symbol"]
                if c in tab.columns), None)
    if col is None:
        raise RuntimeError(f"no gene symbol column in {path.name}: {list(tab.columns)[:15]}")
    if id_col is not None and id_col not in tab.columns:
        raise RuntimeError(f"{id_col} not in {path.name}")
    s = tab.set_index(id_col or tab.columns[0])[col].dropna().astype(str).str.strip()
    s = s[(s != "") & ~s.str.contains(r"///|//|,|;")]
    s.index = s.index.astype(str)
    return s


class HGNC:
    def __init__(self, path: Path):
        t = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
        t = t[t["status"].fillna("Approved").str.contains("Approved")]
        self.approved = set(t["symbol"])
        self.ens = {e: s for e, s in zip(t["ensembl_gene_id"], t["symbol"]) if isinstance(e, str)}
        alias = {}
        for col in ("prev_symbol", "alias_symbol"):
            for s, a in zip(t["symbol"], t[col]):
                if isinstance(a, str):
                    for x in a.split("|"):
                        alias.setdefault(x.strip(), set()).add(s)
        self.alias = {k: next(iter(v)) for k, v in alias.items() if len(v) == 1 and k not in self.approved}

    def map(self, ids) -> pd.Series:
        out = []
        for g in map(str, ids):
            g0 = g.split(".")[0] if g.startswith("ENSG") else g.strip()
            out.append(self.ens.get(g0) if g0.startswith("ENSG") else
                       (g0 if g0 in self.approved else self.alias.get(g0)))
        return pd.Series(out, index=list(ids), dtype=object)


_ILMN_MAPS: dict = {}


def illumina_maps(hg: "HGNC") -> dict:
    """Illumina HumanHT-12 (GPL10558, GPL6947) probe id and array-address id -> HGNC symbol."""
    if not _ILMN_MAPS:
        for name, idc in (("ILMN probe id", None), ("Illumina array address id", "Array_Address_Id")):
            parts = []
            for g in C.ILLUMINA_GPLS[:2]:
                try:
                    parts.append(read_platform(download_gpl(g), idc))
                except RuntimeError:
                    pass
            if parts:
                sm = pd.concat(parts)
                sm = sm[~sm.index.duplicated()]
                sm.index = sm.index.str.replace(r"\.0$", "", regex=True)
                _ILMN_MAPS[name] = pd.Series(hg.map(sm.values).values, index=sm.index).dropna()
    return _ILMN_MAPS


def map_features(ids, hg: "HGNC", extra: dict | None = None) -> tuple[pd.Series, str]:
    """Feature ids -> HGNC symbols; tries Illumina probe / address ids, Ensembl / symbols and annotation columns,
    and keeps the mapping that covers most features."""
    ids = pd.Index([str(i).replace(".0", "") if re.fullmatch(r"\d+\.0", str(i)) else str(i) for i in ids])
    cands = {"HGNC symbol / Ensembl": hg.map(ids)}
    if ids.str.startswith("ILMN_").mean() > 0.3 or ids.str.fullmatch(r"\d+").mean() > 0.3:
        for name, m in illumina_maps(hg).items():
            cands[name] = pd.Series(m.reindex(ids).values, index=ids)
    for name, col in (extra or {}).items():
        cands[f"annotation column '{name}'"] = pd.Series(hg.map(col.astype(str)).values, index=ids)
    best = max(cands, key=lambda k: cands[k].notna().sum())
    return cands[best].dropna(), f"{best} ({cands[best].notna().sum()}/{len(ids)} features)"


def quantile_normalise(x: pd.DataFrame) -> pd.DataFrame:
    ranks = x.rank(axis=0, method="first").astype(int) - 1
    ref = np.sort(x.to_numpy(), axis=0).mean(axis=1)
    return pd.DataFrame(ref[ranks.to_numpy()], index=x.index, columns=x.columns)


def collapse(expr: pd.DataFrame, symbols: pd.Series) -> pd.DataFrame:
    """Rows -> HGNC symbols; several rows per gene: the row with the highest mean expression is kept."""
    sym = symbols.reindex(expr.index)
    e = expr[sym.notna()].copy()
    e["__g"] = sym[sym.notna()].values
    e["__m"] = e.drop(columns="__g").mean(axis=1)
    e = e.sort_values("__m", ascending=False).drop_duplicates("__g")
    return e.set_index("__g").drop(columns="__m").rename_axis("gene")


def to_log2(expr: pd.DataFrame, kind: str) -> tuple[pd.DataFrame, str]:
    v = expr.to_numpy(dtype=float)
    finite = v[np.isfinite(v)]
    if kind == "rnaseq":
        integer_like = np.mean(np.isclose(finite, np.round(finite))) > 0.99
        if integer_like:
            lib = expr.sum(axis=0)
            return np.log2(expr.div(lib, axis=1) * 1e6 + 1), "counts -> log2(CPM+1)"
        return (np.log2(expr.clip(lower=0) + 1), "normalised values -> log2(x+1)") if finite.max() > 100 else \
            (expr, "already log scale")
    if finite.max() > 100:
        return np.log2(expr.clip(lower=1)), "array intensities -> log2"
    return expr, "array, already log2"


def read_table_any(path: Path) -> pd.DataFrame:
    sep = "," if path.name.lower().replace(".gz", "").endswith(".csv") else "\t"
    return pd.read_csv(path, sep=sep, index_col=0, low_memory=False)


def match_columns(cols, md: pd.DataFrame) -> tuple[dict, str]:
    """Match expression column names to sample ids using the metadata field with the most exact matches."""
    cols = [str(c) for c in cols]
    best, best_field = {}, ""
    cand = [md.index.to_series()] + [md[c] for c in md.columns
                                     if md[c].dtype == object or pd.api.types.is_string_dtype(md[c])]
    for s in cand:
        lut = {}
        for sid, v in zip(md.index, s.astype(str)):
            for token in {v.strip(), *[t.strip() for t in re.split(r"[|;,]", v)]}:
                if token:
                    lut.setdefault(token, sid)
        m = {c: lut[c] for c in cols if c in lut}
        if len(m) > len(best):
            best, best_field = m, (s.name if s.name else "id")
    return best, best_field


def read_sdrf(paths: list) -> pd.DataFrame:
    frames = [pd.read_csv(p, sep="\t", dtype=str) for p in paths]
    t = pd.concat(frames, ignore_index=True).drop_duplicates()
    md = pd.DataFrame(index=range(len(t)))
    for c in t.columns:
        m = re.match(r"(Characteristics|Factor Value|Comment)\s*\[(.+)\]", c, re.I)
        key = f"ch:{m.group(2).strip().lower()}" if m else c
        if key in md.columns:
            md[key] = md[key].fillna(t[c])
        else:
            md[key] = t[c].values
    sid = "Source Name" if "Source Name" in md.columns else md.columns[0]
    md = md.drop_duplicates(subset=[sid])
    md.index = md[sid].astype(str).values
    md.index.name = "sample"
    return md


def read_ae_processed(paths: list, md: pd.DataFrame) -> pd.DataFrame:
    """Processed data as one matrix or one file per hybridisation (zip archives are opened)."""
    tables = []
    for p in paths:
        if p.suffix == ".zip":
            with zipfile.ZipFile(p) as z:
                for n in z.namelist():
                    if n.endswith("/"):
                        continue
                    with z.open(n) as fh:
                        tables.append((Path(n).name, pd.read_csv(fh, sep="\t", index_col=0, low_memory=False)))
        else:
            tables.append((p.name, pd.read_csv(p, sep="\t", index_col=0, low_memory=False)))
    tables = [(n, t[~t.index.astype(str).str.startswith(("Reporter", "Hybridization", "!"))]) for n, t in tables]
    mats = [(n, t) for n, t in tables if t.shape[1] >= 10]
    if mats:
        n, t = max(mats, key=lambda x: x[1].shape[1])
        num = t.apply(pd.to_numeric, errors="coerce")
        keep = num.notna().mean() > 0.9
        annot = t.loc[:, ~keep]                                         # annotation columns (ids, symbols)
        num = num.loc[:, keep].dropna(how="all")
        if num.empty:
            raise RuntimeError(f"no numeric expression in {n}: columns {list(t.columns[:6])}, "
                               f"first values {t.iloc[:2, :3].values.tolist()}")
        num.attrs["annot"] = annot.loc[num.index]
        num.attrs["preview"] = f"{n}: index {list(t.index[:3])}, columns {list(t.columns[:4])}"
        return num
    # one file per hybridisation: first numeric column, named by the SDRF row that references the file
    fcol = next((c for c in md.columns if re.search(r"Derived Array Data File", c, re.I)), None)
    by_file = dict(zip(md[fcol].astype(str), md.index)) if fcol else {}
    cols = {}
    for n, t in tables:
        num = t.apply(pd.to_numeric, errors="coerce").dropna(axis=1, how="all")
        if num.shape[1] and n in by_file:
            cols[by_file[n]] = num.iloc[:, 0]
    if not cols:
        n, t = tables[0]
        raise RuntimeError(f"processed files could not be linked to SDRF samples; first file {n}: "
                           f"columns {list(t.columns[:5])}, index {list(t.index[:3])}")
    return pd.DataFrame(cols)


# =====================================================================================
# PRE-SPECIFIED METADATA RULES
# =====================================================================================
RX_AGE_KEY = re.compile(r"(^|[^a-z])age([^a-z]|$)|age_years|age \(y", re.I)
RX_SEX_KEY = re.compile(r"\b(sex|gender)\b", re.I)
RX_DEATH_KEY = re.compile(r"mortal|death|died|surviv|outcome|vital|deceased|alive", re.I)
RX_TIME_KEY = re.compile(r"time.?to|days?.?to|follow|survival.?time|time.?(death|event)", re.I)
RX_SEV_KEY = re.compile(r"sofa|apache|saps", re.I)
RX_PATIENT_KEY = re.compile(r"patient|subject|individual|donor|participant", re.I)
RX_TIMEPOINT_KEY = re.compile(r"time.?point|^ch:day|sampling day|^ch:time$|visit|hour", re.I)
RX_HEALTHY = re.compile(r"healthy|volunteer|control subject|normal control|\bhc\b", re.I)
RX_NONINF = re.compile(r"non.?infect|no.?infection|\bsirs\b|non.?septic|no.?sepsis|without infection|"
                       r"uninfected", re.I)
RX_SEPSIS = re.compile(r"seps|septic|pneumonia|^cap$|\bcap\b|peritonitis|bacter", re.I)
RX_NEGATED = re.compile(r"^(no|non|not|without|un)[\s\-_]?", re.I)
RX_SEPSIS_KEY = re.compile(r"seps|septic|infect|pneumonia|\bcap\b|peritonitis", re.I)


def parse_number(v) -> float:
    m = re.search(r"-?\d+(\.\d+)?", str(v))
    return float(m.group()) if m else np.nan


def parse_death(key: str, v: str) -> float:
    v = str(v).strip()
    if not v or v.lower() in ("na", "n/a", "nan", "unknown", "none"):
        return np.nan
    survival_key = re.search(r"surviv|alive", key, re.I) and not re.search(r"non.?surviv|time", key, re.I)
    if re.search(r"non.?surviv|died|dead|death|deceased|expired", v, re.I):
        return 1.0
    if re.search(r"\bsurviv|alive|discharged", v, re.I):
        return 0.0
    if re.fullmatch(r"(yes|true|1)", v, re.I):
        return 0.0 if survival_key else 1.0
    if re.fullmatch(r"(no|false|0)", v, re.I):
        return 1.0 if survival_key else 0.0
    return np.nan


def classify_group(fields: dict, default: str = "unclassified") -> str:
    """Rule on the characteristics of one sample (each value judged separately): sepsis (positive, non-negated
    value, or a sepsis/infection field set to yes/1/true) > non-infectious > healthy > cohort default.
    Free-text fields (title, description) are not used: they often describe the whole study."""
    vals = [(k, v.strip()) for k, v in fields.items() if v and v.strip() and v.strip().lower() not in ("nan", "na")]
    for k, v in vals:
        if RX_SEPSIS.search(v) and not RX_NEGATED.match(v) and not RX_NONINF.search(v) and not RX_HEALTHY.search(v):
            return "sepsis"
        if RX_SEPSIS_KEY.search(k) and re.fullmatch(r"(yes|true|1|y)", v, re.I) and not RX_NONINF.search(k):
            return "sepsis"
    if any(RX_NONINF.search(v) for _, v in vals):
        return "non-infectious ICU"
    if any(RX_HEALTHY.search(v) for _, v in vals):
        return "healthy"
    return default


def harmonise(acc: str, md: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Pre-specified rules: one row per sample with group, age, sex, death, time, severity, patient, time point."""
    ch = [c for c in md.columns if c.startswith("ch:")]
    out = pd.DataFrame(index=md.index)
    out["cohort"] = acc
    pick = lambda rx, exclude=None: next((c for c in ch if rx.search(c) and not (exclude and exclude.search(c))), None)
    k_age = pick(RX_AGE_KEY, re.compile(r"stage|page|percent|image|dosage", re.I))
    k_sex = pick(RX_SEX_KEY)
    k_time = pick(RX_TIME_KEY)
    k_death = pick(RX_DEATH_KEY, RX_TIME_KEY)
    if k_death is None:          # outcome stored inside a group field (e.g. 'sepsis non-survivor')
        k_death = next((c for c in ch if md[c].astype(str).str.contains(
            r"non.?surviv|surviv|died|alive|deceased", case=False).mean() >= 0.5), None)
    k_sev = [c for c in ch if RX_SEV_KEY.search(c)]
    k_pat = pick(RX_PATIENT_KEY)
    k_tp = pick(RX_TIMEPOINT_KEY, RX_TIME_KEY)
    out["age"] = md[k_age].map(parse_number) if k_age else np.nan
    out["age_field"] = k_age or ""
    out["sex"] = (md[k_sex].astype(str).str.lower().map(
        lambda v: "male" if re.match(r"m(ale)?$", v.strip()) else ("female" if re.match(r"f(emale)?$", v.strip())
                                                                    else np.nan)) if k_sex else np.nan)
    out["death"] = [parse_death(k_death, v) for v in md[k_death]] if k_death else np.nan
    out["death_field"] = k_death or ""
    out["time_days"] = md[k_time].map(parse_number) if k_time else np.nan
    for k in k_sev:
        out[f"sev_{re.sub(r'[^a-z0-9]+', '_', k[3:])}"] = md[k].map(parse_number)
    if k_pat and md[k_pat].nunique() < 0.3 * len(md):       # a group-type field, not a patient identifier
        k_pat = None
    out["patient"] = md[k_pat].astype(str) if k_pat else md.index.astype(str)
    out["patient_field"] = k_pat or "(sample = patient)"
    out["timepoint"] = md[k_tp].map(parse_number) if k_tp else np.nan
    default = "sepsis" if cfg.get("all_sepsis") else "unclassified"     # study design: every patient septic
    out["group"] = [classify_group({f: str(md.at[i, f]) for f in ch}, default) for i in md.index]
    if cfg["role"] == "healthy_reference":
        out["group"] = "healthy"
    out["outcome_definition"] = cfg["outcome"]
    out["tissue"] = cfg["tissue"]
    out["role"] = cfg["role"]
    return out


def first_sample(h: pd.DataFrame) -> pd.Series:
    """True for the earliest sample of each patient (ties: first in the series)."""
    order = h.assign(_o=np.arange(len(h)), _t=h["timepoint"].fillna(np.inf)).sort_values(["patient", "_t", "_o"])
    keep = ~order.duplicated("patient")
    return keep.reindex(h.index)


def col(md: pd.DataFrame, name: str) -> pd.Series:
    return md[name].astype(str).str.strip() if name in md.columns else pd.Series("", index=md.index)


def cohort_rules(acc: str, md: pd.DataFrame, h: pd.DataFrame) -> pd.DataFrame:
    """Amendment 1: cohort-specific rules (override the generic ones where the deposited fields require it)."""
    if acc == "GSE65682":
        endo = col(md, "ch:endotype_class").str.match(r"Mars\d")
        healthy = col(md, "ch:icu_acquired_infection").eq("healthy")
        nonInf = col(md, "ch:abdominal_sepsis_and_controls").eq("ctrl_GI") | col(md, "ch:pneumonia diagnoses").eq("no-cap")
        h["group"] = np.where(endo, "sepsis", np.where(healthy, "healthy",
                                                       np.where(nonInf, "non-infectious ICU", "other ICU")))
    elif acc in ("E-MTAB-5273", "E-MTAB-5274"):
        dis = col(md, "ch:disease").str.lower()
        h["group"] = np.where(dis.str.contains("peritonitis"), "sepsis",
                              np.where(dis.str.contains("pneumonia"), "excluded: CAP (also in E-MTAB-4421/4451)",
                                       np.where(dis.eq("normal"), "healthy", "unclassified")))
        h["patient"] = pd.Series(md.index.astype(str), index=md.index).str.replace(r"\.B(\.\d+)?$", "", regex=True)
        h["patient_field"] = "Source Name without .B.<day>"
        h["timepoint"] = pd.Series(md.index.astype(str), index=md.index).str.extract(r"\.B\.(\d+)$")[0].astype(float).fillna(1)
    elif acc == "GSE95233":
        surv = col(md, "ch:survival")
        h["group"] = np.where(surv.isin(["Survivor", "Non Survivor"]), "sepsis", "healthy")
        h["patient"] = col(md, "title").str.extract(r"-([A-Za-z]+_\d+)_D")[0].fillna(col(md, "title"))
        h["patient_field"] = "title identifier"
        h["timepoint"] = col(md, "ch:time point").str.extract(r"(\d+)")[0].astype(float)
        h["sex"] = col(md, "ch:gender").map({"M": "male", "H": "male", "F": "female"})
    elif acc == "GSE54514":
        h["patient"] = col(md, "ch:group_id")
        h["patient_field"] = "ch:group_id"
        h["timepoint"] = col(md, "ch:group_day").str.extract(r"D(\d+)")[0].astype(float)
    return h


def field_audit(acc: str, md: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for c in md.columns:
        s = md[c].astype(str)
        s = s[~s.isin(["", "nan"])]
        if s.empty:
            continue
        vc = s.value_counts()
        rows.append({"cohort": acc, "field": c, "n_non_missing": len(s), "n_unique": len(vc),
                     "examples": " || ".join(f"{k} (n={v})" for k, v in vc.head(6).items())})
    return pd.DataFrame(rows)


# =====================================================================================
# STEP 2/3 - LOAD, AUDIT, HARMONISE, EXPRESSION
# =====================================================================================
def load_geo_cohort(acc: str, cfg: dict, files, hg: HGNC) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    mats = sorted((C.RAW / acc).glob("*_series_matrix.txt.gz"))
    if not mats:
        raise RuntimeError("no series matrix")
    mds, exprs, notes = [], [], []
    for m in mats:
        md, ex = read_series_matrix(m)
        mds.append(md)
        if len(ex):
            gpl = md["platform_id"].iloc[0]
            if gpl in C.ILLUMINA_GPLS:
                ex.index = pd.Index([str(i) for i in ex.index])
                sym, how = map_features(ex.index, hg)
            else:
                sym = read_platform(download_gpl(gpl))
                sym = pd.Series(hg.map(sym.values).values, index=sym.index).dropna()
                how = f"{gpl} annotation"
            exprs.append(collapse(ex, sym))
            notes.append(f"{m.name}: {gpl}; {how}")
    md = pd.concat(mds)
    md = md[~md.index.duplicated()]
    if exprs:
        expr = pd.concat(exprs, axis=1, join="inner")
    else:                                      # RNA-seq: count table among the supplementary files
        sup = sorted((C.RAW / acc / "suppl").glob("*"))
        cands = [p for p in sup if re.search(r"count|expr|fpkm|tpm|cpm|norm|matrix", p.name, re.I)]
        if not cands:
            raise RuntimeError("RNA-seq series without a count table in 'suppl'")
        best = None
        for p in cands:
            t = read_table_any(p)
            m, fld = match_columns(t.columns, md)
            if best is None or len(m) > len(best[1]):
                best = (t, m, fld, p.name)
        t, m, fld, pname = best
        if len(m) < 0.5 * t.shape[1]:
            raise RuntimeError(f"count-table columns could not be matched to samples ({len(m)}/{t.shape[1]})")
        t = t[list(m)].rename(columns=m)
        t = t.apply(pd.to_numeric, errors="coerce").dropna(how="all")
        sym = hg.map(t.index)
        expr = collapse(t, sym.dropna())
        notes.append(f"{pname}: columns matched via '{fld}' ({len(m)})")
    if acc == "GSE33828":                      # age / sex for the healthy reference are in a supplementary file
        md = add_supplementary_age(md, acc)
    return md, expr, "; ".join(notes)


def load_healthy_supplementary(acc: str, hg: "HGNC") -> tuple[pd.DataFrame, pd.DataFrame, str]:
    """GSE33828: raw expression matrix + age/sex file from the GEO supplementary files (Amendment 1)."""
    sup = sorted((C.RAW / acc / "suppl").glob("*"))
    raw = next((p for p in sup if re.search(r"raw", p.name, re.I) and re.search(r"expression", p.name, re.I)), None)
    agef = next((p for p in sup if re.search(r"age", p.name, re.I) and not re.search(r"expression", p.name, re.I)), None)
    if raw is None or agef is None:
        raise RuntimeError(f"supplementary files missing (found: {[p.name for p in sup]})")
    t = pd.read_csv(raw, sep="\t", index_col=0, low_memory=False)
    num = t.apply(pd.to_numeric, errors="coerce")
    keep = num.notna().mean() > 0.9
    annot, num = t.loc[:, ~keep], num.loc[:, keep].dropna(how="all")
    a = pd.read_csv(agef, sep=None, engine="python", dtype=str)
    cols = set(map(str, num.columns))
    id_col = max(a.columns, key=lambda c: a[c].astype(str).isin(cols).sum())
    a = a[a[id_col].astype(str).isin(cols)].copy()
    if len(a) < 0.5 * num.shape[1]:
        raise RuntimeError(f"age file not matched to matrix columns: matrix {list(num.columns[:3])}, "
                           f"age file {a.columns.tolist()} {a.head(2).values.tolist()}")
    a.index = a[id_col].astype(str)
    md = pd.DataFrame(index=a.index)
    for c in a.columns:
        if c != id_col:
            md[f"ch:{c.strip().lower()}"] = a[c].values
    md["title"] = md.index
    num = num[[c for c in num.columns if str(c) in md.index]]
    num.columns = num.columns.astype(str)
    extra = {c: annot.loc[num.index, c] for c in annot.columns} if len(annot.columns) else None
    sym, how = map_features(num.index, hg, extra)
    num.index = pd.Index([str(i) for i in num.index])
    if len(sym) < 5000:
        raise RuntimeError(f"only {len(sym)} features mapped ({how}); matrix index {list(t.index[:3])}")
    x = num.clip(lower=1)
    x = np.log2(x) if np.nanmax(x.to_numpy()) > 100 else x
    x = quantile_normalise(x.fillna(x.median(axis=1).mean()))
    return md, collapse(x, sym), f"{raw.name}: log2 + quantile normalisation; genes via {how}; age from {agef.name}"


def add_supplementary_age(md: pd.DataFrame, acc: str) -> pd.DataFrame:
    for p in sorted((C.RAW / acc / "suppl").glob("*")):
        if not re.search(r"age|gender|sex", p.name, re.I):
            continue
        t = pd.read_csv(p, sep=None, engine="python", dtype=str)
        id_col = max(t.columns, key=lambda c: len(match_columns(t[c].astype(str), md)[0]))
        m, fld = match_columns(t[id_col].astype(str), md)
        if len(m) < 0.5 * len(t):
            LOG.warning(f"{acc}: supplementary age file {p.name} could not be matched to samples")
            continue
        t = t[t[id_col].astype(str).isin(m)].copy()
        t.index = t[id_col].astype(str).map(m)
        for c in t.columns:
            if c != id_col:
                md[f"ch:{c.strip().lower()}"] = t[c].reindex(md.index)
        LOG.info(f"{acc}: age/sex added from {p.name} (matched via '{fld}', {len(m)} samples)")
        break
    return md


def load_ae_cohort(acc: str, cfg: dict, files, hg: HGNC) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    sdrf = sorted((C.RAW / acc).glob("*sdrf.txt"))
    proc = [p for p in sorted((C.RAW / acc).glob("*")) if not p.name.endswith(("sdrf.txt", "idf.txt", ".part"))]
    if not sdrf or not proc:
        raise RuntimeError("SDRF or processed data missing")
    md = read_sdrf(sdrf)
    ex = read_ae_processed(proc, md)
    m, fld = match_columns(ex.columns, md)
    if len(m) < 0.5 * ex.shape[1]:
        raise RuntimeError(f"processed columns not matched to SDRF ({len(m)}/{ex.shape[1]})")
    annot, preview = ex.attrs.get("annot"), ex.attrs.get("preview", "")
    ex = ex[list(m)].rename(columns=m)
    extra = {c: annot[c] for c in annot.columns} if annot is not None and len(annot.columns) else None
    sym, how = map_features(ex.index, hg, extra)
    ex.index = pd.Index([str(i).replace(".0", "") if re.fullmatch(r"\d+\.0", str(i)) else str(i) for i in ex.index])
    if len(sym) < 5000:
        raise RuntimeError(f"only {len(sym)} features mapped to genes ({how}); file preview: {preview}")
    return md, collapse(ex, sym), f"processed columns matched via '{fld}' ({len(m)}); genes via {how}"


def coverage(expr_genes: set, genes: list) -> float:
    return float(np.mean([g in expr_genes for g in genes])) if genes else 0.0


def load_hallmark() -> list:
    p = C.RAW / "reference" / "MSigDB_Hallmark_2020.txt"
    if not p.exists():
        return []
    for line in p.read_text().splitlines():
        parts = line.split("\t")
        if parts[0].strip() == C.HALLMARK_TERM:
            return [g.split(",")[0].strip() for g in parts[2:] if g.strip()]
    return []


GAINS_CAP_IDS: set = set()
OVERLAP: dict = {}


def step23_main() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    hg = HGNC(download_hgnc())
    hallmark = load_hallmark()
    sets = dict(C.GENE_SETS)
    if hallmark:
        sets["inflammation_hallmark"] = hallmark
        (C.PREREG / "frozen_hallmark_inflammatory_response.txt").write_text("\n".join(hallmark))
    audits, clin, summ, cov = [], [], [], []
    for acc, cfg in C.COHORTS.items():
        try:
            if acc == "GSE33828":
                md, expr, note = load_healthy_supplementary(acc, hg)
            else:
                loader = load_geo_cohort if cfg["source"] == "geo" else load_ae_cohort
                md, expr, note = loader(acc, cfg, None, hg)
            LOG.info(f"{acc}: metadata {md.shape[0]} samples x {md.shape[1]} fields; expression {expr.shape}")
        except Exception as exc:
            LOG.warning(f"{acc}: not available for analysis ({exc}). Manual: {manual_instructions(acc, cfg)}")
            sup = sorted(p.name for p in (C.RAW / acc / "suppl").glob("*")) if (C.RAW / acc / "suppl").exists() else []
            if sup:
                LOG.warning(f"{acc}: supplementary files present: {sup}")
            for m in sorted((C.RAW / acc).glob("*_series_matrix.txt.gz")):     # metadata still audited
                audits.append(field_audit(acc, read_series_matrix(m)[0]))
            summ.append({"cohort": acc, "role": cfg["role"], "status": f"not available: {exc}"})
            continue
        audits.append(field_audit(acc, md))
        pd.concat(audits, ignore_index=True).to_csv(C.TAB / "S1_metadata_field_audit.tsv", sep="\t", index=False)
        (C.TAB / "metadata").mkdir(exist_ok=True)
        md.to_csv(C.TAB / "metadata" / f"metadata_{acc}.tsv", sep="\t")
        try:
            h = cohort_rules(acc, md, harmonise(acc, md, cfg))
            h = h.loc[h.index.intersection(expr.columns)]
            expr = expr[h.index]
            if expr.size == 0 or not np.isfinite(expr.to_numpy(dtype=float)).any():
                raise RuntimeError(f"no expression values linked to samples (expression {expr.shape}; "
                                   f"metadata {md.shape})")
            expr, scale = to_log2(expr, cfg["data"])
        except Exception as exc:
            LOG.warning(f"{acc}: metadata/expression could not be combined ({exc}); field audit written")
            summ.append({"cohort": acc, "role": cfg["role"], "status": f"not combined: {exc}"})
            continue
        h["first_sample"] = first_sample(h)
        elig = (h["group"].isin(["sepsis"]) if cfg["role"] != "healthy_reference" else h["group"].eq("healthy"))
        elig &= h["first_sample"] & (h["age"] >= C.MIN_AGE)
        if cfg["role"] != "healthy_reference":
            elig &= h["death"].isin([0.0, 1.0])
        h["eligible"] = elig
        if acc in ("E-MTAB-4421", "E-MTAB-4451"):
            GAINS_CAP_IDS.update(h["patient"].astype(str))
        if acc in ("E-MTAB-5273", "E-MTAB-5274"):            # Amendment 1: independent patients only
            ov = h["eligible"] & h["patient"].astype(str).isin(GAINS_CAP_IDS)
            OVERLAP[acc] = int(ov.sum())
            h.loc[ov, "eligible"] = False
            h.loc[ov, "group"] = "excluded: patient also in E-MTAB-4421/4451"
        clin.append(h)
        expr.astype("float32").to_csv(C.PROC / f"expr_{acc}.tsv.gz", sep="\t", float_format="%.4f")
        genes = set(expr.index)
        for k, v in sets.items():
            cov.append({"cohort": acc, "gene_set": k, "coverage": coverage(genes, v), "n_genes": len(v)})
        e = h[h["eligible"]]
        summ.append({"cohort": acc, "role": cfg["role"], "status": "loaded", "data": cfg["data"],
                     "tissue": cfg["tissue"], "outcome_definition": cfg["outcome"], "samples": len(h),
                     "groups": "; ".join(f"{g} {n}" for g, n in h["group"].value_counts().items()),
                     "patients": h["patient"].nunique(), "eligible": int(len(e)),
                     "events": int(e["death"].sum()) if cfg["role"] != "healthy_reference" else 0,
                     "time_to_event": bool(e["time_days"].notna().mean() > 0.9) if len(e) else False,
                     "age_median": float(e["age"].median()) if len(e) else np.nan,
                     "age_iqr": f"{e['age'].quantile(.25):.0f}-{e['age'].quantile(.75):.0f}" if len(e) else "",
                     "female_pct": float((e["sex"] == "female").mean() * 100) if len(e) else np.nan,
                     "severity_fields": ", ".join(c for c in h.columns if c.startswith("sev_")),
                     "age_field": h["age_field"].iloc[0], "death_field": h["death_field"].iloc[0],
                     "patient_field": h["patient_field"].iloc[0], "genes": len(genes), "scale": scale,
                     "note": note})
        LOG.info(f"{acc}: {len(h)} samples, {h['patient'].nunique()} patients, eligible {len(e)}, "
                 f"genes {len(genes)}; {scale}; groups: {dict(h['group'].value_counts())}")
    audit = pd.concat(audits, ignore_index=True) if audits else pd.DataFrame()
    audit.to_csv(C.TAB / "S1_metadata_field_audit.tsv", sep="\t", index=False)
    summ = pd.DataFrame(summ)
    summ.to_csv(C.TAB / "S2_cohort_summary.tsv", sep="\t", index=False)
    clin = pd.concat(clin) if clin else pd.DataFrame()
    clin.rename_axis("sample").to_csv(C.TAB / "S3_harmonised_clinical.tsv", sep="\t")
    cov = pd.DataFrame(cov)
    cov.to_csv(C.TAB / "S4_gene_set_coverage.tsv", sep="\t", index=False)
    return summ, clin, cov


# =====================================================================================
# STEP 4 - GATES AND FIGURE
# =====================================================================================
def step4_main(summ: pd.DataFrame, clin: pd.DataFrame, cov: pd.DataFrame) -> None:
    S = summ.set_index("cohort")
    loaded = S[S["status"] == "loaded"]
    rows = []
    d = "GSE65682"
    if d in loaded.index:
        r = loaded.loc[d]
        rows.append({"gate": f"discovery {d}: eligible sepsis patients / deaths",
                     "value": f"{int(r.eligible)} / {int(r.events)}",
                     "criterion": f">= {C.DISCOVERY_MIN_N} / >= {C.DISCOVERY_MIN_EVENTS}",
                     "passed": bool(r.eligible >= C.DISCOVERY_MIN_N and r.events >= C.DISCOVERY_MIN_EVENTS),
                     "severity": "critical"})
        rows.append({"gate": f"discovery {d}: time-to-event available", "value": bool(r.time_to_event),
                     "criterion": "reported (Cox if True, logistic otherwise)", "passed": True, "severity": "warning"})
        noninf = clin[(clin["cohort"] == d) & (clin["group"] == "non-infectious ICU")]
        rows.append({"gate": f"{d}: non-infectious ICU patients with death label (secondary S2)",
                     "value": int(noninf["death"].notna().sum()), "criterion": ">= 30 (else S2 not testable)",
                     "passed": bool(noninf["death"].notna().sum() >= 30), "severity": "warning"})
    else:
        rows.append({"gate": f"discovery {d} loaded", "value": "no", "criterion": "loaded", "passed": False,
                     "severity": "critical"})
    val = loaded[loaded["role"] == "validation"]
    usable = val[(val["eligible"] >= C.VALIDATION_MIN_N) & (val["events"] >= C.VALIDATION_MIN_EVENTS)]
    for acc, r in val.iterrows():
        ok = bool(r.eligible >= C.VALIDATION_MIN_N and r.events >= C.VALIDATION_MIN_EVENTS)
        rows.append({"gate": f"validation {acc}: eligible / deaths ({r.outcome_definition})",
                     "value": f"{int(r.eligible)} / {int(r.events)}",
                     "criterion": f">= {C.VALIDATION_MIN_N} / >= {C.VALIDATION_MIN_EVENTS} (else excluded)",
                     "passed": ok, "severity": "warning"})
    for acc in S.index[S["status"] != "loaded"]:
        rows.append({"gate": f"{acc} available", "value": S.loc[acc, "status"][:80], "criterion": "loaded",
                     "passed": False, "severity": "critical" if S.loc[acc, "role"] != "validation" else "warning"})
    for acc, n in OVERLAP.items():
        rows.append({"gate": f"{acc}: faecal-peritonitis patients removed because also in E-MTAB-4421/4451",
                     "value": n, "criterion": "reported (Amendment 1: independent cohorts)", "passed": True,
                     "severity": "warning"})
    rows.append({"gate": "validation cohorts usable / total deaths", "value": f"{len(usable)} / {int(usable['events'].sum())}",
                 "criterion": f">= {C.VALIDATION_MIN_COHORTS} / >= {C.VALIDATION_MIN_TOTAL_EVENTS}",
                 "passed": bool(len(usable) >= C.VALIDATION_MIN_COHORTS
                                and usable["events"].sum() >= C.VALIDATION_MIN_TOTAL_EVENTS), "severity": "critical"})
    if "GSE33828" in loaded.index:
        r = loaded.loc["GSE33828"]
        rows.append({"gate": "healthy reference GSE33828: donors with age", "value": int(r.eligible),
                     "criterion": f">= {C.HEALTHY_MIN_N}", "passed": bool(r.eligible >= C.HEALTHY_MIN_N),
                     "severity": "critical"})
    for acc, r in loaded.iterrows():
        rows.append({"gate": f"{acc}: genes mapped to HGNC symbols", "value": int(r.genes),
                     "criterion": f">= {C.MIN_GENES}", "passed": bool(r.genes >= C.MIN_GENES),
                     "severity": "critical" if acc in ("GSE65682", "GSE33828") else "warning"})
    low = cov[cov["coverage"] < C.MIN_SET_COVERAGE]
    rows.append({"gate": "pre-specified gene sets with coverage below threshold (cohort: set)",
                 "value": "; ".join(f"{a}: {b}" for a, b in zip(low["cohort"], low["gene_set"])) or "none",
                 "criterion": f"coverage >= {int(C.MIN_SET_COVERAGE * 100)}% (else score not computed there)",
                 "passed": low.empty, "severity": "warning"})
    rows.append({"gate": "MSigDB Hallmark inflammatory response gene set available",
                 "value": int((cov["gene_set"] == "inflammation_hallmark").any()), "criterion": "available",
                 "passed": bool((cov["gene_set"] == "inflammation_hallmark").any()), "severity": "warning"})
    record_gates("stage1", rows)
    plot_fig1(summ, clin, cov)


def plot_fig1(summ: pd.DataFrame, clin: pd.DataFrame, cov: pd.DataFrame) -> None:
    set_style()
    S = summ[summ["status"] == "loaded"].copy()
    order = ["discovery", "validation", "healthy_reference"]
    S["o"] = S["role"].map({r: i for i, r in enumerate(order)})
    S = S.sort_values(["o", "eligible"], ascending=[True, False])
    fig = plt.figure(figsize=(FIG_W, 78 * MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 1.05], wspace=0.55)
    # a: patients and deaths
    ax = fig.add_subplot(gs[0, 0])
    y = np.arange(len(S))[::-1]
    role_col = {"discovery": OKABE[0], "validation": OKABE[5], "healthy_reference": "#999999"}
    ax.barh(y, S["eligible"], color=[role_col[r] for r in S["role"]], height=0.7)
    ax.barh(y, S["events"], color=OKABE[1], height=0.7)
    for yy, (_, r) in zip(y, S.iterrows()):
        txt = f"{int(r.eligible)}" + (f" ({int(r.events)} deaths)" if r.role != "healthy_reference" else "")
        ax.text(r.eligible + S["eligible"].max() * 0.02, yy, txt, va="center", fontsize=5)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{c}\n{d}, {t}" for c, d, t in zip(S["cohort"], S["data"], S["tissue"])], fontsize=5)
    ax.set_xlim(0, S["eligible"].max() * 1.45)
    ax.set_xlabel("Patients / donors analysed (one sample each)")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (OKABE[0], OKABE[5], "#999999", OKABE[1])]
    ax.legend(handles, ["discovery", "validation", "healthy reference", "deaths"], loc="upper center",
              bbox_to_anchor=(0.45, -0.2), ncol=4, frameon=False, fontsize=5, handlelength=1, columnspacing=0.8)
    panel_label(ax, "a", x=-0.55)
    # b: age distributions of the analysed patients (no outcome information shown)
    ax = fig.add_subplot(gs[0, 1])
    data = [clin[(clin["cohort"] == c) & clin["eligible"]]["age"].dropna().values for c in S["cohort"]]
    keep = [i for i, v in enumerate(data) if len(v)]
    kw = dict(positions=[y[i] for i in keep], widths=0.6, patch_artist=True, showfliers=False,
              medianprops={"color": "black", "lw": 0.8}, boxprops={"lw": 0.5}, whiskerprops={"lw": 0.5},
              capprops={"lw": 0.5})
    try:
        bp = ax.boxplot([data[i] for i in keep], orientation="horizontal", **kw)     # matplotlib >= 3.10
    except TypeError:
        bp = ax.boxplot([data[i] for i in keep], vert=False, **kw)
    for patch, i in zip(bp["boxes"], keep):
        patch.set_facecolor(role_col[S["role"].iloc[i]])
        patch.set_alpha(0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_ylim(-0.7, len(S) - 0.3)
    ax.set_xlabel("Age (years)")
    ax.set_title("Age of analysed patients", fontsize=6)
    panel_label(ax, "b", x=-0.12)
    # c: gene-set coverage
    ax = fig.add_subplot(gs[0, 2])
    M = cov.pivot(index="cohort", columns="gene_set", values="coverage").reindex(S["cohort"])
    im = ax.imshow(M.values, aspect="auto", cmap="Blues", vmin=0, vmax=1, interpolation="none")
    ax.set_yticks(range(len(M)))
    ax.set_yticklabels(M.index, fontsize=5)
    ax.set_xticks(range(M.shape[1]))
    nice = {"b_cell": "B cells", "classical_monocyte": "Classical monocytes", "immature_neutrophil":
            "Immature neutrophils", "inflammation_hallmark": "Inflammation (Hallmark)", "mature_neutrophil":
            "Mature neutrophils", "mhc_class_ii": "MHC class II", "ms1_monocyte": "MS1 monocytes",
            "nk_cell": "NK cells", "t_cell": "T cells"}
    ax.set_xticklabels([nice.get(c, c) for c in M.columns], rotation=55, ha="right", fontsize=5)
    for (i, j), v in np.ndenumerate(M.values):          # only incomplete coverage is printed
        if np.isfinite(v) and v < 0.995:
            ax.text(j, i, f"{v * 100:.0f}", ha="center", va="center", fontsize=5,
                    color="white" if v > 0.6 else "black")
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    cb.set_label("Genes measured (%)", fontsize=5)
    cb.ax.tick_params(labelsize=5)
    cb.set_ticks([0, 0.5, 1])
    cb.set_ticklabels(["0", "50", "100"])
    ax.set_title("Coverage of pre-specified gene sets\n(numbers shown if < 100%)", fontsize=6)
    panel_label(ax, "c", x=-0.45)
    save_figure(fig, "Fig1_cohorts")


# =====================================================================================
# RUN
# =====================================================================================
def run_stage1() -> None:
    np.random.seed(C.SEED)
    LOG.info(f"Project folder: {C.PROJECT_DIR}")
    write_preregistration()
    write_amendment()
    step1_main()
    summ, clin, cov = step23_main()
    step4_main(summ, clin, cov)
    (C.LOGS / "stage1_run_info.json").write_text(json.dumps({
        "finished_utc": utc_now(), "python": sys.version.split()[0], "pandas": pd.__version__,
        "numpy": np.__version__}, indent=1))
    LOG.info(gate_summary())


run_stage1()
