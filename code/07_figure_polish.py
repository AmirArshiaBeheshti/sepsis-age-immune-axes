# =====================================================================================
# SEPSIS PROJECT - 07 FIGURE POLISH: redraws Fig 4 (panel a new), Fig 5 (panel b new) and Fig 7 (panels a and f new)
#   ONLY from the saved result tables (S19, S21, S22b-S26, S34-S36 and the cached single-cell object). No analysis is
#   re-run and no result can change; the other panels are drawn by exactly the same code as before.
#   The previous figure files are copied to results/figures/before_polish/ before they are replaced.
#
# HOW TO RUN: paste the whole file into ONE Jupyter cell in the same working directory as 01-06. No downloads.
#   Run time: under a minute.
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


# #####################################################################################
# PART B - ANALYSIS STAGE 1: immune-state scores, prognosis meta-analysis (H1) and age-effect mediation (H2)
# #####################################################################################
import warnings
import statsmodels.api as sm
from scipy import stats
from scipy.special import expit

LOG = get_logger("sepsis_analysis1")
GATES_FILE = C.TAB / "S0_analysis1_quality_gates.tsv"
_GATES.clear()


class A1:
    AGE_STEP = 10.0                       # exposure contrast: 10 years older
    N_BOOT = 1000
    SCORES = ["immature_neutrophil", "mature_neutrophil", "ms1_monocyte", "classical_monocyte", "t_cell", "b_cell",
              "nk_cell", "mhc_class_ii", "inflammation_hallmark"]
    PRIMARY_MEDIATORS = ["ml_index", "mhc_class_ii", "inflammation_hallmark"]
    NO_RISK_SCALE = ("E-MTAB-4451",)
    MIN_EVENTS_INTERACTION = 40      # outcome-dependent sampling: excluded from risk-difference mediation
    ALPHA = 0.05


AMENDMENT_3 = """# Amendment 3 to Pre-registration v1 - revised primary hypotheses (written before ANY outcome analysis)

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
"""


def write_amendment_3() -> None:
    p, meta = C.PREREG / "AMENDMENT_3.md", C.PREREG / "AMENDMENT_3.json"
    h = hashlib.sha256(AMENDMENT_3.encode()).hexdigest()
    if p.exists():
        old = json.loads(meta.read_text())
        if hashlib.sha256(p.read_text().encode()).hexdigest() != old["sha256"] or h != old["sha256"]:
            raise SystemExit("AMENDMENT_3.md differs from its registered version - do not edit it.")
        LOG.info(f"Amendment 3 present (registered {old['registered_utc']}, unchanged)")
        return
    if (C.TAB / "S10_immune_scores.tsv").exists():
        raise SystemExit("Scores already exist but Amendment 3 is missing - the amendment must precede the analysis.")
    p.write_text(AMENDMENT_3)
    meta.write_text(json.dumps({"registered_utc": utc_now(), "sha256": h, "file": p.name}, indent=1))
    LOG.info(f"Amendment 3 written (SHA-256 {h[:16]}...) - before any outcome analysis")


# -------------------------------------------------------------------------------------
# data and scores
# -------------------------------------------------------------------------------------
def analysed_cohorts() -> tuple[str, list]:
    s = pd.read_csv(C.TAB / "S2_cohort_summary.tsv", sep="\t")
    s = s[s["status"] == "loaded"]
    disc = s.loc[s["role"] == "discovery", "cohort"].iloc[0]
    val = s[(s["role"] == "validation") & (s["eligible"] >= C.VALIDATION_MIN_N) &
            (s["events"] >= C.VALIDATION_MIN_EVENTS)]["cohort"].tolist()
    return disc, val


def platform_of(cohort: str) -> str:
    s = pd.read_csv(C.TAB / "S2_cohort_summary.tsv", sep="\t").set_index("cohort")
    note = str(s.loc[cohort, "note"])
    if C.COHORTS[cohort]["data"] == "rnaseq":
        return "RNA-seq"
    return "Affymetrix" if re.search(r"GPL13667|GPL570", note) else "Illumina"


def hallmark_genes() -> list:
    p = C.PREREG / "frozen_hallmark_inflammatory_response.txt"
    return [g for g in p.read_text().split() if g] if p.exists() else []


def cohort_scores(cohort: str, clin: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    h = clin[(clin["cohort"] == cohort) & clin["eligible"].astype(bool)].copy()
    E = pd.read_csv(C.PROC / f"expr_{cohort}.tsv.gz", sep="\t", index_col=0)
    E = E[h.index.intersection(E.columns)]
    h = h.loc[E.columns]
    sd = E.std(axis=1, ddof=0)
    E = E[sd > 0]
    Z = E.sub(E.mean(axis=1), axis=0).div(E.std(axis=1, ddof=0), axis=0)
    sets = dict(C.GENE_SETS)
    hm = hallmark_genes()
    if hm:
        sets["inflammation_hallmark"] = hm
    cov = {}
    for k in A1.SCORES:
        genes = [g for g in sets.get(k, []) if g in Z.index]
        cov[k] = len(genes) / max(len(sets.get(k, [])), 1)
        if cov[k] >= C.MIN_SET_COVERAGE:
            s = Z.loc[genes].mean(axis=0)
            h[k] = ((s - s.mean()) / s.std(ddof=0)).values
        else:
            h[k] = np.nan
    ml = h["immature_neutrophil"] - h["t_cell"]
    h["ml_index"] = (ml - ml.mean()) / ml.std(ddof=0)
    h["age10"] = h["age"] / A1.AGE_STEP
    h["male"] = (h["sex"] == "male").astype(float).where(h["sex"].isin(["male", "female"]))
    if h["male"].notna().mean() < 0.8 or h["male"].nunique() < 2:      # sex not usable: not adjusted in this cohort
        LOG.warning(f"{cohort}: sex missing or constant - models in this cohort are not adjusted for sex")
        h["male"] = 0.0
    h["platform"] = platform_of(cohort)
    return h, cov


# -------------------------------------------------------------------------------------
# models
# -------------------------------------------------------------------------------------
def logit_fit(y, X):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = sm.Logit(y, sm.add_constant(X, has_constant="add")).fit(disp=0, method="newton", maxiter=200)
    if not r.mle_retvals.get("converged", True):
        raise RuntimeError("logistic model did not converge")
    return r


def assoc(d: pd.DataFrame, pred: str, adjust=("age10", "male"), model: str = "logit") -> dict:
    cols = [pred] + [a for a in adjust if a != pred and not (a == "male" and d["male"].nunique() < 2)]
    dd = d.dropna(subset=cols + ["death"] + (["time_days"] if model == "cox" else []))
    X = dd[cols].astype(float)
    try:
        if model == "cox":
            t = dd["time_days"].clip(lower=0.5).astype(float)
            r = sm.PHReg(t, X, status=dd["death"].astype(int), ties="breslow").fit(disp=0)
            b, se = float(r.params[0]), float(r.bse[0])
        else:
            r = logit_fit(dd["death"].astype(int), X)
            b, se = float(r.params[pred]), float(r.bse[pred])
    except Exception as exc:                       # separation / non-convergence
        LOG.warning(f"model failed ({pred}, {model}, n={len(dd)}): {exc}")
        return {"n": len(dd), "events": int(dd["death"].sum()), "logratio": np.nan, "se": np.nan}
    return {"n": len(dd), "events": int(dd["death"].sum()), "logratio": b, "se": se}


def dl_meta(est: np.ndarray, se: np.ndarray) -> dict:
    """DerSimonian-Laird random-effects meta-analysis."""
    ok = np.isfinite(est) & np.isfinite(se) & (se > 0)
    est, se = np.asarray(est)[ok], np.asarray(se)[ok]
    k = len(est)
    if k == 0:
        return {"k": 0, "est": np.nan, "se": np.nan, "P": np.nan, "I2": np.nan, "tau2": np.nan}
    w = 1 / se ** 2
    fe = np.sum(w * est) / np.sum(w)
    q = np.sum(w * (est - fe) ** 2)
    tau2 = max(0.0, (q - (k - 1)) / (np.sum(w) - np.sum(w ** 2) / np.sum(w))) if k > 1 else 0.0
    wr = 1 / (se ** 2 + tau2)
    m = np.sum(wr * est) / np.sum(wr)
    s = np.sqrt(1 / np.sum(wr))
    i2 = max(0.0, (q - (k - 1)) / q) if k > 1 and q > 0 else 0.0
    return {"k": k, "est": m, "se": s, "P": 2 * stats.norm.sf(abs(m / s)), "I2": i2, "tau2": tau2}


def g_formula(d: pd.DataFrame, meds: list) -> dict:
    """Natural direct / indirect / total effect (risk difference) of +10 years of age, mediators jointly,
    exposure-mediator interaction (VanderWeele regression-based approach with g-computation)."""
    y = d["death"].astype(int).to_numpy()
    a = d["age10"].to_numpy(float)
    use_sex = d["male"].nunique() > 1
    sex = d["male"].to_numpy(float).reshape(len(d), -1) if use_sex else np.zeros((len(d), 0))
    M = d[meds].to_numpy(float)
    # mediator models: M_k ~ age + sex
    Xm = np.column_stack([np.ones(len(a)), a, sex])          # sex column omitted if not usable
    beta = np.linalg.lstsq(Xm, M, rcond=None)[0]          # 3 x K
    shift = beta[1]                                       # per 10 years
    # outcome model: death ~ age + M + age x M + sex
    inter = y.sum() >= A1.MIN_EVENTS_INTERACTION
    design = (lambda age_, Mx: np.column_stack([np.ones(len(age_)), age_, Mx, age_[:, None] * Mx, sex])) if inter \
        else (lambda age_, Mx: np.column_stack([np.ones(len(age_)), age_, Mx, sex]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = sm.Logit(y, design(a, M)).fit(disp=0, method="newton", maxiter=200)
    if not r.mle_retvals.get("converged", True) or not np.all(np.isfinite(r.params)):
        raise RuntimeError("outcome model did not converge")
    g = r.params

    def risk(age_, Mx):
        return float(np.mean(expit(design(age_, Mx) @ g)))
    a1, M1 = a + 1.0, M + shift
    y11, y10, y00 = risk(a1, M1), risk(a1, M), risk(a, M)
    te, nde, nie = y11 - y00, y10 - y00, y11 - y10
    return {"TE": te, "NDE": nde, "NIE": nie, "PM": nie / te if abs(te) > 1e-9 else np.nan}


def mediation(d: pd.DataFrame, meds: list, n_boot: int, seed: int) -> dict:
    cols = meds + ["age10", "male", "death"]
    d = d.dropna(subset=cols)
    est = g_formula(d, meds)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(d), len(d))
        try:
            boots.append(g_formula(d.iloc[idx], meds))
        except Exception:
            continue
    B = pd.DataFrame(boots)
    out = {"n": len(d), "events": int(d["death"].sum()), "boots_ok": len(B),
           "interaction": bool(d["death"].sum() >= A1.MIN_EVENTS_INTERACTION)}
    for k in ("TE", "NDE", "NIE", "PM"):
        out[k] = est[k]
        out[f"{k}_se"] = float(B[k].std(ddof=1)) if len(B) > 10 else np.nan
        out[f"{k}_lo"], out[f"{k}_hi"] = (np.nanpercentile(B[k], [2.5, 97.5]) if len(B) > 10 else (np.nan, np.nan))
    return out


# -------------------------------------------------------------------------------------
# figures
# -------------------------------------------------------------------------------------
def log_axis(ax) -> None:
    ax.set_xscale("log")
    lo, hi = ax.get_xlim()
    grid = (0.125, 0.25, 0.5, 1, 2, 4, 8) if hi / lo > 3 else (0.5, 0.67, 0.8, 1, 1.25, 1.5, 2)
    ticks = [t for t in grid if lo <= t <= hi]
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())


def forest(ax, rows: pd.DataFrame, title: str, xlabel: str) -> None:
    y = np.arange(len(rows))[::-1]
    for yy, (_, r) in zip(y, rows.iterrows()):
        pooled = str(r["label"]).startswith(("Validation", "All"))
        col_ = OKABE[1] if pooled else (OKABE[0] if r.get("role") == "discovery" else "#555555")
        lo, hi = np.exp(r["est"] - 1.96 * r["se"]), np.exp(r["est"] + 1.96 * r["se"])
        ax.plot([lo, hi], [yy, yy], color=col_, lw=0.8)
        ax.scatter([np.exp(r["est"])], [yy], marker="D" if pooled else "s", s=18 if pooled else 10, color=col_,
                   zorder=3)
        ax.text(1.02, yy, f"{np.exp(r['est']):.2f} ({lo:.2f}-{hi:.2f})", transform=ax.get_yaxis_transform(),
                va="center", fontsize=5)
    ax.axvline(1, color="#777777", lw=0.5, ls="--")
    log_axis(ax)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{l}" for l in rows["label"]], fontsize=5)
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=6)


def plot_fig2(h1: pd.DataFrame, comp: pd.DataFrame, agesc: pd.DataFrame) -> None:
    set_style()
    fig = plt.figure(figsize=(FIG_W, 78 * MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.8, 0.8], wspace=1.6, left=0.12, right=0.98, top=0.84,
                          bottom=0.2)
    ax = fig.add_subplot(gs[0, 0])
    forest(ax, h1, "H1: ML index and death\n(adjusted for age and sex)", "OR per SD (HR for GSE65682)")
    panel_label(ax, "a", x=-0.62, y=1.15)
    ax = fig.add_subplot(gs[0, 1])
    c = comp.sort_values("est")
    y = np.arange(len(c))
    ax.errorbar(np.exp(c["est"]), y, xerr=[np.exp(c["est"]) - np.exp(c["est"] - 1.96 * c["se"]),
                                            np.exp(c["est"] + 1.96 * c["se"]) - np.exp(c["est"])],
                fmt="o", ms=3, color=OKABE[0], elinewidth=0.7, capsize=0)
    ax.axvline(1, color="#777777", lw=0.5, ls="--")
    log_axis(ax)
    ax.set_yticks(y)
    ax.set_yticklabels([NICE.get(s, s) for s in c["score"]], fontsize=5)
    ax.set_xlabel("Pooled OR per SD (all cohorts)")
    ax.set_title("Each immune score and death\n(random effects, age and sex adjusted)", fontsize=6)
    panel_label(ax, "b", x=-0.62, y=1.15)
    ax = fig.add_subplot(gs[0, 2])
    a = agesc.set_index("score").reindex(c["score"])
    ax.errorbar(a["est"], y, xerr=1.96 * a["se"], fmt="o", ms=3, color=OKABE[2], elinewidth=0.7, capsize=0)
    ax.axvline(0, color="#777777", lw=0.5, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_xlabel("SD change per 10 years of age")
    ax.set_title("Age and immune state in sepsis\n(pooled, sex adjusted)", fontsize=6)
    panel_label(ax, "c", x=-0.15, y=1.15)
    save_figure(fig, "Fig2_prognosis_and_age")


def plot_fig3(med: pd.DataFrame, pooled: dict, comp_med: pd.DataFrame) -> None:
    set_style()
    fig = plt.figure(figsize=(FIG_W, 72 * MM))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.75, left=0.12, right=0.97, top=0.85, bottom=0.2)
    ax = fig.add_subplot(gs[0, 0])
    rows = med.copy()
    y = np.arange(len(rows))[::-1]
    for off, k, col_ in ((0.18, "TE", "#555555"), (0.0, "NDE", OKABE[5]), (-0.18, "NIE", OKABE[1])):
        ax.errorbar(rows[k] * 100, y + off, xerr=[(rows[k] - rows[f"{k}_lo"]) * 100, (rows[f"{k}_hi"] - rows[k]) * 100],
                    fmt="o", ms=2.5, color=col_, elinewidth=0.6, capsize=0,
                    label={"TE": "total", "NDE": "direct", "NIE": "indirect (via immune state)"}[k])
    ax.axvline(0, color="#777777", lw=0.5, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(rows["label"], fontsize=5)
    ax.set_xlabel("Risk difference for 10 years older (percentage points)")
    ax.set_title("H2: decomposition of the age effect on death", fontsize=6)
    ax.legend(frameon=False, fontsize=5, loc="upper center", bbox_to_anchor=(0.45, -0.2), ncol=3)
    panel_label(ax, "a", x=-0.45, y=1.06)
    ax = fig.add_subplot(gs[0, 1])
    c = comp_med.sort_values("PM_est")
    y = np.arange(len(c))
    ax.errorbar(c["PM_est"] * 100, y, xerr=[(c["PM_est"] - c["PM_lo"]) * 100, (c["PM_hi"] - c["PM_est"]) * 100],
                fmt="o", ms=3, color=OKABE[1], elinewidth=0.7, capsize=0)
    ax.axvline(0, color="#777777", lw=0.5, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([NICE.get(s, s) for s in c["mediator"]], fontsize=5)
    ax.set_xlabel("Proportion of the age effect mediated (%)")
    t = (f"Primary (joint): {pooled['PM'] * 100:.0f}% ({pooled['PM_lo'] * 100:.0f} to {pooled['PM_hi'] * 100:.0f})"
         if np.isfinite(pooled.get("PM", np.nan)) else "Primary (joint): not estimable")
    ax.set_title(f"Single mediators (pooled validation)\n{t}", fontsize=6)
    panel_label(ax, "b", x=-0.6, y=1.06)
    save_figure(fig, "Fig3_age_mediation")


NICE = {"ml_index": "ML index", "immature_neutrophil": "Immature neutrophils",
        "mature_neutrophil": "Mature neutrophils", "ms1_monocyte": "MS1 monocytes",
        "classical_monocyte": "Classical monocytes", "t_cell": "T cells", "b_cell": "B cells", "nk_cell": "NK cells",
        "mhc_class_ii": "MHC class II", "inflammation_hallmark": "Inflammation (Hallmark)"}


# -------------------------------------------------------------------------------------
# run
# -------------------------------------------------------------------------------------
def pooled_ratio(nie: dict, te: dict, n: int = 200000, seed: int = 0) -> tuple:
    if not (np.isfinite(nie.get("est", np.nan)) and np.isfinite(te.get("est", np.nan)) and
            np.isfinite(nie.get("se", np.nan)) and np.isfinite(te.get("se", np.nan))):
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(seed)
    a = rng.normal(nie["est"], nie["se"], n)
    b = rng.normal(te["est"], te["se"], n)
    r = a / b
    return nie["est"] / te["est"], *np.nanpercentile(r[np.abs(b) > 1e-6], [2.5, 97.5])


# numbers produced by the data-preparation run on the real data (S2_cohort_summary, 2 Nov 2026); a mismatch means
# that the folder holds other / older data -> warning (the analysis still uses what is in the folder)
EXPECTED = {"GSE65682": (478, 114), "E-MTAB-4421": (265, 56), "E-MTAB-4451": (106, 52), "E-MTAB-5273": (94, 14),
            "GSE95233": (51, 17), "GSE54514": (36, 10), "GSE185263": (345, 52)}


def to_bool(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


def preflight() -> None:
    """Stops early, with a clear message, if anything this stage needs is missing or inconsistent."""
    need = {C.TAB / "S2_cohort_summary.tsv": "data-preparation stage (sepsis_stage1_jupyter.py)",
            C.TAB / "S3_harmonised_clinical.tsv": "data-preparation stage (sepsis_stage1_jupyter.py)",
            C.TAB / "S0_stage1_quality_gates.tsv": "data-preparation stage (sepsis_stage1_jupyter.py)",
            C.PREREG / "PREREG_v1.md": "data-preparation stage", C.PREREG / "AMENDMENT_1.md": "data-preparation stage"}
    miss = [f"{p.name} (made by the {w})" for p, w in need.items() if not p.exists()]
    if miss:
        raise SystemExit(f"Project folder {C.PROJECT_DIR} is missing: {'; '.join(miss)}. Run this cell in the same "
                         f"working directory as before, or set SEPSIS_PROJECT_DIR to that folder.")
    g = pd.read_csv(C.TAB / "S0_stage1_quality_gates.tsv", sep="\t")
    if (g["status"] == "FAIL").any():
        raise SystemExit("The data-preparation quality gates contain FAIL rows - fix those first:\n"
                         + g.loc[g["status"] == "FAIL", "gate"].to_string(index=False))
    if not (C.PREREG / "frozen_hallmark_inflammatory_response.txt").exists():
        LOG.warning("Frozen Hallmark gene set missing: the inflammation score and H2 (which uses it) will be incomplete")
    disc, val = analysed_cohorts()
    for c in [disc] + val:
        f = C.PROC / f"expr_{c}.tsv.gz"
        if not f.exists():
            raise SystemExit(f"Expression file {f.name} is missing - re-run the data-preparation stage.")
    clin = pd.read_csv(C.TAB / "S3_harmonised_clinical.tsv", sep="\t", index_col=0, low_memory=False)
    for col_ in ("cohort", "eligible", "age", "sex", "death", "time_days", "group"):
        if col_ not in clin.columns:
            raise SystemExit(f"S3_harmonised_clinical.tsv lacks column '{col_}' - re-run the data-preparation stage.")
    el = clin[to_bool(clin["eligible"])]
    for c, (n, e) in EXPECTED.items():
        sub = el[el["cohort"] == c]
        if c in [disc] + val and (len(sub), int(sub["death"].sum())) != (n, e):
            LOG.warning(f"{c}: {len(sub)} patients / {int(sub['death'].sum())} deaths in this folder; expected {n} / {e} "
                        f"from the verified run - check that this is the intended project folder")
    import statsmodels
    LOG.info(f"Preflight passed: discovery {disc}, {len(val)} validation cohorts; "
             f"python {sys.version.split()[0]}, pandas {pd.__version__}, numpy {np.__version__}, "
             f"statsmodels {statsmodels.__version__}")


def data_hash(d: pd.DataFrame, cols: list) -> str:
    return hashlib.sha256(pd.util.hash_pandas_object(d[cols].round(10), index=True).values.tobytes()).hexdigest()[:16]


def cached_mediation(d: pd.DataFrame, meds: list, n_boot: int, seed: int, tag: str) -> dict:
    """Mediation with on-disk cache: a re-run (e.g. after an interruption) reuses finished cohorts."""
    cols = meds + ["age10", "male", "death"]
    key = f"{tag}_{'+'.join(meds)}_{n_boot}_{seed}_{data_hash(d.dropna(subset=cols), cols)}"
    p = C.PROJECT_DIR / "results" / "cache" / "mediation" / f"{key}.json"
    if p.exists():
        return json.loads(p.read_text())
    r = mediation(d, meds, n_boot, seed)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps({k: (float(v) if isinstance(v, (np.floating, float)) else v) for k, v in r.items()}))
    return r


def safe_plot(fn, *args) -> None:
    try:
        fn(*args)
    except Exception as exc:                     # a figure problem must not lose the computed results
        LOG.warning(f"figure {fn.__name__} failed ({exc}); all tables were written")


def run_analysis1() -> None:
    t0 = time.time()
    warnings.filterwarnings("ignore", category=FutureWarning)
    np.random.seed(C.SEED)
    LOG.info(f"Project folder: {C.PROJECT_DIR}")
    preflight()
    write_preregistration()
    write_amendment()
    write_amendment_3()
    clin = pd.read_csv(C.TAB / "S3_harmonised_clinical.tsv", sep="\t", index_col=0, low_memory=False)
    clin["eligible"] = to_bool(clin["eligible"])
    disc, val = analysed_cohorts()
    LOG.info(f"Discovery {disc}; validation {', '.join(val)}")
    data, gates = {}, []
    for c in [disc] + val:
        try:
            d, cov = cohort_scores(c, clin)
        except Exception as exc:
            if c == disc:
                raise SystemExit(f"Discovery cohort {c} could not be scored: {exc}")
            LOG.warning(f"{c}: scores could not be computed ({exc}) - cohort excluded")
            gates.append({"gate": f"{c}: scores computed", "value": str(exc)[:80], "criterion": "computed",
                          "passed": False, "severity": "warning"})
            continue
        data[c] = d
        low = [k for k, v in cov.items() if v < C.MIN_SET_COVERAGE]
        LOG.info(f"{c}: {len(d)} patients, {int(d['death'].sum())} deaths; platform {d['platform'].iloc[0]}"
                 + (f"; scores not computed (coverage): {low}" if low else ""))
        gates.append({"gate": f"{c}: complete cases (age, sex, ML index, death)",
                      "value": int(d.dropna(subset=["age10", "male", "ml_index", "death"]).shape[0]),
                      "criterion": f">= {C.VALIDATION_MIN_N}", "severity": "critical" if c == disc else "warning",
                      "passed": bool(d.dropna(subset=["age10", "male", "ml_index", "death"]).shape[0] >= C.VALIDATION_MIN_N)})
    val = [c for c in val if c in data]
    cohorts = [disc] + val
    allsc = pd.concat(data.values())
    allsc.to_csv(C.TAB / "S10_immune_scores.tsv", sep="\t")

    # ---------------- H1 ----------------
    rows = []
    for c in cohorts:
        model = "cox" if (c == disc and data[c]["time_days"].notna().mean() > 0.9) else "logit"
        r = assoc(data[c], "ml_index", model=model)
        rows.append({"label": c, "role": "discovery" if c == disc else "validation", "model": model, **r})
        if c == disc:
            rows.append({"label": f"{c} (logistic, sensitivity)", "role": "sensitivity", "model": "logit",
                         **assoc(data[c], "ml_index")})
    H1 = pd.DataFrame(rows)
    vm = dl_meta(H1.loc[H1["role"] == "validation", "logratio"].to_numpy(),
                 H1.loc[H1["role"] == "validation", "se"].to_numpy())
    am = dl_meta(H1.loc[H1["role"].isin(["validation"]) | (H1["label"] == f"{disc} (logistic, sensitivity)"),
                        "logratio"].to_numpy(),
                 H1.loc[H1["role"].isin(["validation"]) | (H1["label"] == f"{disc} (logistic, sensitivity)"),
                        "se"].to_numpy())
    H1 = pd.concat([H1, pd.DataFrame([
        {"label": f"Validation pooled (k={vm['k']}, I2={vm['I2'] * 100:.0f}%)", "role": "pooled", "model": "RE",
         "logratio": vm["est"], "se": vm["se"], "P": vm["P"]},
        {"label": f"All cohorts pooled (k={am['k']}, I2={am['I2'] * 100:.0f}%)", "role": "pooled", "model": "RE",
         "logratio": am["est"], "se": am["se"], "P": am["P"]}])], ignore_index=True)
    H1["P"] = H1["P"].fillna(pd.Series(2 * stats.norm.sf(np.abs(H1["logratio"] / H1["se"])), index=H1.index))
    H1["ratio"], H1["lo"], H1["hi"] = (np.exp(H1["logratio"]), np.exp(H1["logratio"] - 1.96 * H1["se"]),
                                       np.exp(H1["logratio"] + 1.96 * H1["se"]))
    H1.to_csv(C.TAB / "S11_H1_ml_index_mortality.tsv", sep="\t", index=False)
    disc_b = H1.loc[H1["label"] == disc, "logratio"].iloc[0]
    LOG.info(f"H1 discovery {disc}: ratio {np.exp(disc_b):.2f}; validation pooled OR {np.exp(vm['est']):.2f} "
             f"({np.exp(vm['est'] - 1.96 * vm['se']):.2f}-{np.exp(vm['est'] + 1.96 * vm['se']):.2f}), P = {vm['P']:.2g}")

    # sensitivity: leave-one-out, leave-4451-out, platform strata, severity-adjusted
    sens = []
    vrows = H1[H1["role"] == "validation"]
    for c in vrows["label"]:
        r_ = vrows[vrows["label"] != c]
        m_ = dl_meta(r_["logratio"].to_numpy(), r_["se"].to_numpy())
        sens.append({"analysis": f"validation without {c}", "OR": np.exp(m_["est"]), "P": m_["P"], "k": m_["k"]})
    for plat in ("Affymetrix", "Illumina", "RNA-seq"):
        cs = [c for c in val if data[c]["platform"].iloc[0] == plat]
        if cs:
            r_ = vrows[vrows["label"].isin(cs)]
            m_ = dl_meta(r_["logratio"].to_numpy(), r_["se"].to_numpy())
            sens.append({"analysis": f"validation, {plat} only", "OR": np.exp(m_["est"]), "P": m_["P"], "k": m_["k"]})
    for c in cohorts:
        sevc = [x for x in data[c].columns if x.startswith("sev_") and data[c][x].notna().mean() > 0.7]
        for sv in sevc:
            d_ = data[c].copy()
            d_["sev_z"] = (d_[sv] - d_[sv].mean()) / d_[sv].std()
            r = assoc(d_, "ml_index", adjust=("age10", "male", "sev_z"))
            sens.append({"analysis": f"{c}: adjusted for {sv[4:]}", "OR": np.exp(r["logratio"]),
                         "P": 2 * stats.norm.sf(abs(r["logratio"] / r["se"])) if r["se"] == r["se"] else np.nan,
                         "k": 1})
    # interaction age x ML index (pooled over all cohorts)
    ib = []
    for c in cohorts:
        d_ = data[c].dropna(subset=["ml_index", "age10", "male", "death"]).copy()
        d_["age_x_ml"] = d_["age10"] * d_["ml_index"]
        r = assoc(d_, "age_x_ml", adjust=("ml_index", "age10", "male"))
        ib.append((r["logratio"], r["se"]))
    im = dl_meta(np.array([x[0] for x in ib]), np.array([x[1] for x in ib]))
    sens.append({"analysis": "age (per 10 y) x ML index interaction, all cohorts", "OR": np.exp(im["est"]),
                 "P": im["P"], "k": im["k"]})
    pd.DataFrame(sens).to_csv(C.TAB / "S12_H1_sensitivity.tsv", sep="\t", index=False)

    # component scores (secondary) and age -> score (descriptive)
    comp, agesc = [], []
    for s in ["ml_index"] + A1.SCORES:
        est, se, ae, as_ = [], [], [], []
        for c in cohorts:
            if data[c][s].notna().sum() < C.VALIDATION_MIN_N:
                continue
            r = assoc(data[c], s)
            est.append(r["logratio"])
            se.append(r["se"])
            d_ = data[c].dropna(subset=[s, "age10", "male"])
            xc = ["age10", "male"] if d_["male"].nunique() > 1 else ["age10"]
            ols = sm.OLS(d_[s], sm.add_constant(d_[xc])).fit()
            ae.append(ols.params["age10"])
            as_.append(ols.bse["age10"])
        m_ = dl_meta(np.array(est), np.array(se))
        comp.append({"score": s, **m_})
        m2 = dl_meta(np.array(ae), np.array(as_))
        agesc.append({"score": s, **m2})
    comp, agesc = pd.DataFrame(comp), pd.DataFrame(agesc)
    comp.to_csv(C.TAB / "S13_scores_mortality_meta.tsv", sep="\t", index=False)
    agesc.to_csv(C.TAB / "S14_age_scores_meta.tsv", sep="\t", index=False)

    # ---------------- H2 ----------------
    LOG.info(f"H2: mediation with {A1.N_BOOT} bootstraps per cohort (mediators: {', '.join(A1.PRIMARY_MEDIATORS)})")
    med = []
    for i, c in enumerate(cohorts):
        if c in A1.NO_RISK_SCALE:
            LOG.info(f"  {c}: excluded from risk-scale mediation (outcome-dependent sampling)")
            continue
        if any(data[c][m].isna().all() for m in A1.PRIMARY_MEDIATORS):
            LOG.warning(f"  {c}: a primary mediator is missing - excluded")
            continue
        r = cached_mediation(data[c], A1.PRIMARY_MEDIATORS, A1.N_BOOT, C.SEED + i, c)
        med.append({"label": c, "role": "discovery" if c == disc else "validation", **r})
        LOG.info(f"  {c}: TE {r['TE'] * 100:.1f} pp, NIE {r['NIE'] * 100:.1f} pp "
                 f"({r['NIE_lo'] * 100:.1f} to {r['NIE_hi'] * 100:.1f}); PM {r['PM'] * 100:.0f}%")
    med = pd.DataFrame(med)
    v = med[med["role"] == "validation"]
    nie_m, te_m = dl_meta(v["NIE"].to_numpy(), v["NIE_se"].to_numpy()), dl_meta(v["TE"].to_numpy(), v["TE_se"].to_numpy())
    nde_m = dl_meta(v["NDE"].to_numpy(), v["NDE_se"].to_numpy())
    pm, pm_lo, pm_hi = pooled_ratio(nie_m, te_m) if te_m["k"] else (np.nan, np.nan, np.nan)
    pooled = {"PM": pm, "PM_lo": pm_lo, "PM_hi": pm_hi}
    med = pd.concat([med, pd.DataFrame([{
        "label": f"Validation pooled (k={nie_m['k']})", "role": "pooled",
        "TE": te_m["est"], "TE_lo": te_m["est"] - 1.96 * te_m["se"], "TE_hi": te_m["est"] + 1.96 * te_m["se"],
        "NDE": nde_m["est"], "NDE_lo": nde_m["est"] - 1.96 * nde_m["se"], "NDE_hi": nde_m["est"] + 1.96 * nde_m["se"],
        "NIE": nie_m["est"], "NIE_lo": nie_m["est"] - 1.96 * nie_m["se"], "NIE_hi": nie_m["est"] + 1.96 * nie_m["se"],
        "NIE_P": nie_m["P"], "PM": pm, "PM_lo": pm_lo, "PM_hi": pm_hi, "I2_NIE": nie_m["I2"]}])], ignore_index=True)
    med.to_csv(C.TAB / "S15_H2_mediation.tsv", sep="\t", index=False)
    LOG.info(f"H2 validation pooled: NIE {nie_m['est'] * 100:.2f} pp (P = {nie_m['P']:.2g}); "
             f"proportion mediated {pm * 100:.0f}% ({pm_lo * 100:.0f} to {pm_hi * 100:.0f})")
    # single mediators (secondary)
    cm = []
    for j, s in enumerate(["ml_index"] + A1.SCORES):
        LOG.info(f"  single mediator {j + 1}/{len(A1.SCORES) + 1}: {s}")
        r_ = []
        for i, c in enumerate(val):
            if c in A1.NO_RISK_SCALE or data[c][s].isna().all():
                continue
            r_.append(cached_mediation(data[c], [s], A1.N_BOOT, C.SEED + 100 + i, c))
        r_ = pd.DataFrame(r_)
        if r_.empty:
            continue
        n_, t_ = dl_meta(r_["NIE"].to_numpy(), r_["NIE_se"].to_numpy()), dl_meta(r_["TE"].to_numpy(), r_["TE_se"].to_numpy())
        p_, lo_, hi_ = pooled_ratio(n_, t_)
        cm.append({"mediator": s, "NIE": n_["est"], "NIE_P": n_["P"], "PM_est": p_, "PM_lo": lo_, "PM_hi": hi_,
                   "k": n_["k"]})
    cm = pd.DataFrame(cm)
    cm.to_csv(C.TAB / "S16_single_mediators.tsv", sep="\t", index=False)

    # ---------------- Holm and gates ----------------
    d_nie = med.loc[med["label"] == disc, "NIE"]
    p1 = vm["P"] if (np.sign(vm["est"]) == np.sign(disc_b) and vm["est"] > 0) else 1.0
    p2 = nie_m["P"] if (len(d_nie) and np.isfinite(nie_m["est"]) and nie_m["est"] > 0
                        and np.sign(d_nie.iloc[0]) == 1) else 1.0
    p1 = p1 if np.isfinite(p1) else 1.0
    holm = sorted([("H1", p1), ("H2", p2)], key=lambda x: x[1])
    adj, prev = {}, 0.0
    for j, (hname, p) in enumerate(holm):
        prev = max(prev, min(1.0, p * (2 - j)))
        adj[hname] = prev
    res = pd.DataFrame([
        {"hypothesis": "H1 ML index -> death (age, sex adjusted)", "discovery": f"{np.exp(disc_b):.2f}",
         "validation_pooled": f"OR {np.exp(vm['est']):.2f} ({np.exp(vm['est'] - 1.96 * vm['se']):.2f}-"
                              f"{np.exp(vm['est'] + 1.96 * vm['se']):.2f})", "P": vm["P"], "Holm_P": adj["H1"],
         "supported": adj["H1"] < A1.ALPHA},
        {"hypothesis": "H2 age effect mediated by immune state (NIE, risk difference)",
         "discovery": f"{d_nie.iloc[0] * 100:.2f} pp" if len(d_nie) else "n/a",
         "validation_pooled": f"NIE {nie_m['est'] * 100:.2f} pp; PM {pm * 100:.0f}% ({pm_lo * 100:.0f} to {pm_hi * 100:.0f})",
         "P": nie_m["P"], "Holm_P": adj["H2"], "supported": adj["H2"] < A1.ALPHA}])
    res.to_csv(C.TAB / "S17_primary_results.tsv", sep="\t", index=False)
    LOG.info("PRIMARY RESULTS (pre-registered, Holm over 2):\n" + res.to_string(index=False))
    for _, r in res.iterrows():
        gates.append({"gate": f"RESULT (reported, not a gate): {r['hypothesis']}",
                      "value": f"{r['validation_pooled']}; Holm P = {r['Holm_P']:.2g}; supported = {r['supported']}",
                      "criterion": "reported", "passed": True, "severity": "warning"})
    bmin = med.loc[med["role"] != "pooled", "boots_ok"].min() if "boots_ok" in med.columns else np.nan
    gates.append({"gate": "mediation bootstraps converged (min over cohorts)",
                  "value": int(bmin) if np.isfinite(bmin) else "none",
                  "criterion": f">= {int(0.95 * A1.N_BOOT)}",
                  "passed": bool(np.isfinite(bmin) and bmin >= 0.95 * A1.N_BOOT), "severity": "warning"})
    record_gates("analysis1", gates)
    h1_plot = H1[H1["role"].isin(["discovery", "validation", "pooled"])].rename(columns={"logratio": "est"})
    safe_plot(plot_fig2, h1_plot, comp, agesc)
    safe_plot(plot_fig3, med[["label", "role", "TE", "TE_lo", "TE_hi", "NDE", "NDE_lo", "NDE_hi", "NIE", "NIE_lo",
                              "NIE_hi"]], pooled, cm)
    LOG.info(f"Run time {(time.time() - t0) / 60:.1f} min")
    LOG.info(gate_summary().replace("STAGE 1 QUALITY GATES", "ANALYSIS STAGE 1 QUALITY GATES").replace(
        "Stage 2 can start after the field audit has been checked.", "Analysis stage 1 complete.").replace(
        "do not start Stage 2", "check the failed gates"))


# #####################################################################################
# PART C - ANALYSIS STAGE 2: gene-level meta-analysis, pathway enrichment (pre-ranked GSEA) and bulk
#          ligand-receptor axes
# #####################################################################################
LOG = get_logger("sepsis_analysis2")
GATES_FILE = C.TAB / "S0_analysis2_quality_gates.tsv"
_GATES.clear()


class A2:
    ENRICHR = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={lib}"
    REACTOME_LIBS = ("Reactome_Pathways_2024", "Reactome_2022")
    SET_MIN, SET_MAX = 15, 500
    N_PERM = 10000
    FDR = 0.05
    AGE_P = 0.05
    SIG_GENE_P = 0.05
    MIN_SIG_GENES = 20
    # pre-specified ligand-receptor axes (Amendment 3); a side with several genes = mean of their z
    LR_PAIRS = {
        "HLA class II - CD4": (["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1"], ["CD4"]),
        "MIF - CD74": (["MIF"], ["CD74"]),
        "CD80 - CD28": (["CD80"], ["CD28"]),
        "CD86 - CD28": (["CD86"], ["CD28"]),
        "CD86 - CTLA4": (["CD86"], ["CTLA4"]),
        "CD274 - PDCD1": (["CD274"], ["PDCD1"]),
        "IL10 - IL10RA/IL10RB": (["IL10"], ["IL10RA", "IL10RB"]),
        "S100A8/S100A9 - TLR4": (["S100A8", "S100A9"], ["TLR4"]),
        "ARG1 (neutrophil suppression)": (["ARG1"], []),
    }


AMENDMENT_4 = """# Amendment 4 - specification of analysis stage 2 (pathways, ligand-receptor axes in bulk)

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
"""


def write_amendment_4() -> None:
    p, meta = C.PREREG / "AMENDMENT_4.md", C.PREREG / "AMENDMENT_4.json"
    h = hashlib.sha256(AMENDMENT_4.encode()).hexdigest()
    if p.exists():
        old = json.loads(meta.read_text())
        if hashlib.sha256(p.read_text().encode()).hexdigest() != old["sha256"] or h != old["sha256"]:
            raise SystemExit("AMENDMENT_4.md differs from its registered version - do not edit it.")
        LOG.info(f"Amendment 4 present (registered {old['registered_utc']}, unchanged)")
        return
    if (C.TAB / "S18_gene_meta.tsv.gz").exists():
        raise SystemExit("Stage-2 results exist but Amendment 4 is missing - the amendment must precede the analysis.")
    p.write_text(AMENDMENT_4)
    meta.write_text(json.dumps({"registered_utc": utc_now(), "sha256": h, "file": p.name}, indent=1))
    LOG.info(f"Amendment 4 written (SHA-256 {h[:16]}...) - before stage 2")


# -------------------------------------------------------------------------------------
# gene sets
# -------------------------------------------------------------------------------------
def read_enrichr_library(path: Path) -> dict:
    out = {}
    for line in path.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) > 2:
            out[parts[0].strip()] = sorted({g.split(",")[0].strip() for g in parts[2:] if g.strip()})
    return out


def gene_sets() -> dict:
    ref = C.RAW / "reference"
    hall = read_enrichr_library(ref / "MSigDB_Hallmark_2020.txt")
    react = {}
    for lib in A2.REACTOME_LIBS:
        p = ref / f"{lib}.txt"
        try:
            fetch(A2.ENRICHR.format(lib=lib), p)
            react = read_enrichr_library(p)
            if len(react) > 100:
                LOG.info(f"Reactome library: {lib} ({len(react)} sets)")
                break
        except Exception as exc:
            LOG.warning(f"{lib} not downloaded ({exc})")
    return {"Hallmark": hall, "Reactome": react}


# -------------------------------------------------------------------------------------
# gene-level statistics (vectorised over genes)
# -------------------------------------------------------------------------------------
def batched_logit(y: np.ndarray, G: np.ndarray, Z: np.ndarray, iters: int = 30) -> tuple[np.ndarray, np.ndarray]:
    """For every gene g: logit P(y) = b0 + b1 * G[:, g] + Z c. Newton-Raphson for all genes at once.
    Returns (b1, se(b1)); NaN where the fit did not converge (separation)."""
    n, ng = G.shape
    X = np.concatenate([np.ones((ng, n, 1)), G.T[:, :, None], np.broadcast_to(Z, (ng, n, Z.shape[1]))], axis=2)
    p = X.shape[2]
    b = np.zeros((ng, p))
    b[:, 0] = np.log(y.mean() / (1 - y.mean()))
    conv = np.zeros(ng, dtype=bool)
    for _ in range(iters):
        eta = np.einsum("gnp,gp->gn", X, b)
        mu = expit(eta)
        w = np.clip(mu * (1 - mu), 1e-10, None)
        grad = np.einsum("gnp,gn->gp", X, y[None, :] - mu)
        H = np.einsum("gnp,gn,gnq->gpq", X, w, X)
        try:
            step = np.linalg.solve(H, grad[:, :, None])[:, :, 0]
        except np.linalg.LinAlgError:
            step = np.einsum("gpq,gq->gp", np.linalg.pinv(H), grad)
        b += step
        conv = np.abs(step).max(axis=1) < 1e-8
        if conv.all():
            break
    eta = np.einsum("gnp,gp->gn", X, b)
    mu = expit(eta)
    H = np.einsum("gnp,gn,gnq->gpq", X, np.clip(mu * (1 - mu), 1e-10, None), X)
    se = np.sqrt(np.abs(np.linalg.pinv(H)[:, 1, 1]))
    bad = ~conv | ~np.isfinite(b[:, 1]) | (np.abs(b[:, 1]) > 10)
    return np.where(bad, np.nan, b[:, 1]), np.where(bad, np.nan, se)


def gene_stats(cohort: str, d: pd.DataFrame, genes: list) -> pd.DataFrame:
    E = pd.read_csv(C.PROC / f"expr_{cohort}.tsv.gz", sep="\t", index_col=0)
    d = d.dropna(subset=["age10", "male", "death"])
    E = E.loc[genes, d.index]
    Zg = E.sub(E.mean(axis=1), axis=0).div(E.std(axis=1, ddof=0).replace(0, np.nan), axis=0).T   # samples x genes
    ok = Zg.notna().all(axis=0).to_numpy()
    covs = ["age10"] + (["male"] if d["male"].nunique() > 1 else [])
    Zc = d[covs].to_numpy(float)
    y = d["death"].to_numpy(float)
    b = np.full(len(genes), np.nan)
    s = np.full(len(genes), np.nan)
    idx = np.where(ok)[0]
    for i0 in range(0, len(idx), 1000):                       # chunks keep memory small
        j = idx[i0:i0 + 1000]
        b[j], s[j] = batched_logit(y, Zg.iloc[:, j].to_numpy(float), Zc)
    # gene ~ age + sex (OLS), all genes at once
    Xa = np.column_stack([np.ones(len(d)), Zc])
    Y = Zg.fillna(0).to_numpy(float)
    XtXi = np.linalg.pinv(Xa.T @ Xa)
    B = XtXi @ Xa.T @ Y
    res = Y - Xa @ B
    s2 = (res ** 2).sum(axis=0) / (len(d) - Xa.shape[1])
    age_b, age_se = B[1], np.sqrt(s2 * XtXi[1, 1])
    return pd.DataFrame({"gene": genes, f"death_b_{cohort}": b, f"death_se_{cohort}": s,
                         f"age_b_{cohort}": np.where(ok, age_b, np.nan), f"age_se_{cohort}": np.where(ok, age_se, np.nan)})


def meta_vectorised(B: np.ndarray, S: np.ndarray) -> pd.DataFrame:
    """DerSimonian-Laird per row (genes x cohorts)."""
    ok = np.isfinite(B) & np.isfinite(S) & (S > 0)
    W = np.where(ok, 1 / np.where(ok, S, 1) ** 2, 0.0)
    Bz = np.where(ok, B, 0.0)
    k = ok.sum(axis=1)
    fe = (W * Bz).sum(axis=1) / np.where(W.sum(axis=1) > 0, W.sum(axis=1), np.nan)
    Q = (W * (Bz - fe[:, None]) ** 2 * ok).sum(axis=1)
    c = W.sum(axis=1) - (W ** 2).sum(axis=1) / np.where(W.sum(axis=1) > 0, W.sum(axis=1), np.nan)
    tau2 = np.where(k > 1, np.maximum(0, (Q - (k - 1)) / np.where(c > 0, c, np.nan)), 0.0)
    Wr = np.where(ok, 1 / (np.where(ok, S, 1) ** 2 + tau2[:, None]), 0.0)
    est = (Wr * Bz).sum(axis=1) / np.where(Wr.sum(axis=1) > 0, Wr.sum(axis=1), np.nan)
    se = np.sqrt(1 / np.where(Wr.sum(axis=1) > 0, Wr.sum(axis=1), np.nan))
    z = est / se
    return pd.DataFrame({"est": est, "se": se, "z": z, "P": 2 * stats.norm.sf(np.abs(z)), "k": k,
                         "I2": np.where((k > 1) & (Q > 0), np.maximum(0, (Q - (k - 1)) / np.where(Q > 0, Q, 1)), 0)})


# -------------------------------------------------------------------------------------
# pre-ranked GSEA (gene-set permutation), implemented in NumPy
# -------------------------------------------------------------------------------------
class PrerankGSEA:
    def __init__(self, ranking: pd.Series, n_perm: int, kmax: int, seed: int):
        r = ranking.dropna().sort_values(ascending=False)
        self.genes = r.index.to_numpy()
        self.pos = {g: i for i, g in enumerate(self.genes)}
        self.w = np.abs(r.to_numpy(float))
        self.n = len(r)
        rng = np.random.default_rng(seed)
        # random gene orders; the first k columns of a row form a random k-subset (shared null for all sizes)
        self.perm = np.vstack([np.argsort(rng.random((500, self.n)), axis=1)[:, :kmax]
                               for _ in range(int(np.ceil(n_perm / 500)))])[:n_perm].astype(np.int32)

    def _es(self, idx: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """ES for one or many sets given as row-wise hit positions; also returns the argmax position."""
        idx = np.sort(np.atleast_2d(idx), axis=1)
        k = idx.shape[1]
        w = self.w[idx]
        nr = w.sum(axis=1, keepdims=True)
        after = np.cumsum(w, axis=1) / nr - (idx - np.arange(k)) / (self.n - k)
        before = after - w / nr
        mx, mn = after.max(axis=1), before.min(axis=1)
        es = np.where(mx >= -mn, mx, mn)
        peak = np.where(mx >= -mn, idx[np.arange(len(idx)), after.argmax(axis=1)],
                        idx[np.arange(len(idx)), before.argmin(axis=1)])
        return es, peak

    def run(self, sets: dict) -> pd.DataFrame:
        null_cache, rows = {}, []
        for name, genes in sets.items():
            hit = np.array(sorted({self.pos[g] for g in genes if g in self.pos}))
            k = len(hit)
            if k < A2.SET_MIN or k > A2.SET_MAX:
                continue
            if k not in null_cache:
                null_cache[k] = self._es(self.perm[:, :k])[0]
            null = null_cache[k]
            es, peak = self._es(hit)
            es, peak = float(es[0]), int(peak[0])
            same = null[null >= 0] if es >= 0 else null[null < 0]
            nes = es / abs(same.mean()) if len(same) else np.nan
            p = ((same >= es).sum() + 1) / (len(same) + 1) if es >= 0 else ((same <= es).sum() + 1) / (len(same) + 1)
            le = [self.genes[i] for i in hit if (i <= peak if es >= 0 else i >= peak)]
            rows.append({"pathway": name, "size": k, "ES": es, "NES": nes, "P": p, "leading_edge": ";".join(le)})
        out = pd.DataFrame(rows)
        if len(out):
            out["FDR"] = np.maximum(self._qvalues(out, null_cache), 1.0 / len(self.perm))   # resolution floor
        return out

    @staticmethod
    def _qvalues(out: pd.DataFrame, null_cache: dict) -> np.ndarray:
        """GSEA FDR q-value (Subramanian et al. 2005): for NES* > 0, q = [mean over sets of the fraction of the
        set's positive null NES >= NES*] / [fraction of observed positive NES >= NES*]; mirrored for NES* < 0."""
        nulls = {}
        for k, null in null_cache.items():
            pos, neg = null[null >= 0], null[null < 0]
            nulls[k] = (np.sort(pos / pos.mean()) if len(pos) else np.array([]),
                        np.sort(-neg / neg.mean()) if len(neg) else np.array([]))   # negative NES, ascending
        sizes = out["size"].to_numpy()
        nes = out["NES"].to_numpy(float)
        q = np.ones(len(out))
        for sign in (1, -1):
            sel = np.where(np.sign(nes) == sign)[0] if sign > 0 else np.where(nes < 0)[0]
            if not len(sel):
                continue
            obs = nes[sel]
            for i, x in zip(sel, obs):
                fr = []
                for k in sizes[np.sign(nes) == sign] if sign > 0 else sizes[nes < 0]:
                    arr = nulls[k][0] if sign > 0 else nulls[k][1]
                    if len(arr):
                        fr.append((len(arr) - np.searchsorted(arr, x, side="left")) / len(arr) if sign > 0
                                  else np.searchsorted(arr, x, side="right") / len(arr))
                num = np.mean(fr) if fr else 1.0
                den = (obs >= x).mean() if sign > 0 else (obs <= x).mean()
                q[i] = min(1.0, num / den) if den > 0 else 1.0
            order = sel[np.argsort(-np.abs(nes[sel]))]                       # most extreme first
            q[order] = np.minimum.accumulate(q[order][::-1])[::-1]           # monotone in |NES|
        return q


# -------------------------------------------------------------------------------------
# ligand-receptor axes in bulk
# -------------------------------------------------------------------------------------
def lr_scores(cohort: str, d: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    E = pd.read_csv(C.PROC / f"expr_{cohort}.tsv.gz", sep="\t", index_col=0)
    E = E[d.index.intersection(E.columns)]
    Z = E.sub(E.mean(axis=1), axis=0).div(E.std(axis=1, ddof=0).replace(0, np.nan), axis=0)
    out, cov = pd.DataFrame(index=E.columns), {}
    for pair, (lig, rec) in A2.LR_PAIRS.items():
        lg, rg = [g for g in lig if g in Z.index], [g for g in rec if g in Z.index]
        cov[pair] = f"{len(lg)}/{len(lig)} + {len(rg)}/{len(rec)}"
        if not lg or (rec and not rg):
            continue
        sides = [Z.loc[lg].mean(axis=0)] + ([Z.loc[rg].mean(axis=0)] if rg else [])
        s = sum(sides) / len(sides)
        out[pair] = ((s - s.mean()) / s.std(ddof=0)).values
    return out, cov


# -------------------------------------------------------------------------------------
# figure
# -------------------------------------------------------------------------------------
def plot_fig4(hall_d: pd.DataFrame, hall_a: pd.DataFrame, react_d: pd.DataFrame, lr: pd.DataFrame) -> None:
    set_style()
    fig = plt.figure(figsize=(FIG_W, 125 * MM))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1], height_ratios=[1, 0.85], wspace=1.45, hspace=0.6,
                          left=0.2, right=0.95, top=0.93, bottom=0.08)
    # a: Hallmark NES for death and age
    ax = fig.add_subplot(gs[0, 0])
    h = hall_d.merge(hall_a[["pathway", "NES", "FDR"]], on="pathway", suffixes=("", "_age"))
    h = h[(h["FDR"] < A2.FDR) | (h["FDR_age"] < A2.FDR)]
    h = pd.concat([h.nsmallest(10, "NES"), h.nlargest(14, "NES")]).drop_duplicates("pathway").sort_values("NES")
    y = np.arange(len(h))
    ax.barh(y, h["NES"], color=np.where(h["NES"] > 0, OKABE[1], OKABE[0]), height=0.65, alpha=0.9)
    ax.scatter(h["NES_age"], y, marker="D", s=10, color="black", zorder=3, label="age NES")
    ax.axvline(0, color="#777777", lw=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([(p[0].upper() + p[1:]).replace("_", " ")[:40] for p in h["pathway"]], fontsize=5)
    for t_, f_ in zip(ax.get_yticklabels(), h["FDR"]):    # black = death FDR < 0.05; grey = age only
        t_.set_color("black" if f_ < A2.FDR else "#888888")
    ax.set_xlabel("NES (bars: death; diamonds: age)\nlabels: black = death FDR < 0.05, grey = age only")
    ax.set_title("Hallmark pathways\n(FDR < 0.05 for death or age)", fontsize=6)
    panel_label(ax, "a", x=-0.9, y=1.03)
    # b: top Reactome death pathways
    ax = fig.add_subplot(gs[0, 1])
    r = react_d[react_d["FDR"] < A2.FDR].copy()
    r = pd.concat([r.nlargest(10, "NES"), r.nsmallest(10, "NES")]).drop_duplicates("pathway").sort_values("NES")
    y = np.arange(len(r))
    sc_ = ax.scatter(r["NES"], y, s=np.clip(r["size"] / 4, 6, 40), c=-np.log10(r["FDR"].clip(lower=1e-10)),
                     cmap="viridis", linewidths=0)
    ax.axvline(0, color="#777777", lw=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([re.sub(r"\s*R-HSA-\d+", "", p)[:46] for p in r["pathway"]], fontsize=4.8)
    ax.set_xlabel("NES (death)")
    ax.set_title("Reactome pathways and death\n(FDR < 0.05; top 10 per direction)", fontsize=6)
    cb = fig.colorbar(sc_, ax=ax, fraction=0.05, pad=0.03)
    cb.set_label("-log10 FDR", fontsize=5)
    cb.ax.tick_params(labelsize=5)
    panel_label(ax, "b", x=-1.35, y=1.03)
    # c, d: ligand-receptor axes
    l_ = lr.sort_values("death_OR")
    y = np.arange(len(l_))
    ax = fig.add_subplot(gs[1, 0])
    for off, col_, est, lo, hi, lab in ((0.12, OKABE[0], "death_OR", "death_lo", "death_hi", "age, sex adjusted"),
                                        (-0.12, OKABE[1], "death_OR_adjML", "adjML_lo", "adjML_hi", "+ ML index")):
        ax.errorbar(l_[est], y + off, xerr=[l_[est] - l_[lo], l_[hi] - l_[est]], fmt="o", ms=2.5, color=col_,
                    elinewidth=0.6, capsize=0, label=lab)
    ax.axvline(1, color="#777777", lw=0.5, ls="--")
    log_axis(ax)
    ax.set_yticks(y)
    ax.set_yticklabels(l_["axis"], fontsize=5)
    ax.set_xlabel("Pooled OR per SD (all cohorts)")
    ax.set_title("Ligand-receptor axes and death", fontsize=6)
    ax.legend(frameon=False, fontsize=5, loc="upper center", bbox_to_anchor=(0.45, -0.22), ncol=2)
    panel_label(ax, "c", x=-0.9, y=1.03)
    ax = fig.add_subplot(gs[1, 1])
    ax.errorbar(l_["age_b"], y, xerr=1.96 * l_["age_se"], fmt="o", ms=2.5, color=OKABE[2], elinewidth=0.6, capsize=0)
    ax.axvline(0, color="#777777", lw=0.5, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_xlabel("SD change per 10 years of age")
    ax.set_title("Ligand-receptor axes and age", fontsize=6)
    panel_label(ax, "d", x=-0.12, y=1.03)
    save_figure(fig, "Fig4_pathways_ligand_receptor")


# -------------------------------------------------------------------------------------
# run
# -------------------------------------------------------------------------------------
def run_analysis2() -> None:
    t0 = time.time()
    warnings.filterwarnings("ignore", category=FutureWarning)
    LOG.info(f"Project folder: {C.PROJECT_DIR}")
    preflight()
    if not (C.TAB / "S17_primary_results.tsv").exists():
        raise SystemExit("Analysis stage 1 has not finished (S17_primary_results.tsv missing) - run 02_analysis_H1_H2.py.")
    write_preregistration()
    write_amendment()
    write_amendment_3()
    write_amendment_4()
    clin = pd.read_csv(C.TAB / "S3_harmonised_clinical.tsv", sep="\t", index_col=0, low_memory=False)
    clin["eligible"] = to_bool(clin["eligible"])
    disc, val = analysed_cohorts()
    cohorts = [disc] + val
    data = {c: cohort_scores(c, clin)[0] for c in cohorts}
    gates = []
    # ---------------- gene-level meta-analysis ----------------
    common = None
    for c in cohorts:
        g = pd.read_csv(C.PROC / f"expr_{c}.tsv.gz", sep="\t", index_col=0, usecols=[0]).index
        common = set(g) if common is None else common & set(g)
    common = sorted(common)
    LOG.info(f"Gene-level models: {len(common)} genes measured in all {len(cohorts)} cohorts")
    gates.append({"gate": "genes measured in all analysed cohorts", "value": len(common), "criterion": ">= 8000",
                  "passed": len(common) >= 8000, "severity": "critical"})
    tabs = []
    for c in cohorts:
        t_ = time.time()
        tabs.append(gene_stats(c, data[c], common).set_index("gene"))
        conv = tabs[-1][f"death_b_{c}"].notna().mean()
        LOG.info(f"  {c}: gene-level models done ({conv:.1%} converged, {time.time() - t_:.0f} s)")
        gates.append({"gate": f"{c}: gene-level logistic models converged", "value": f"{conv:.1%}",
                      "criterion": ">= 95%", "passed": bool(conv >= 0.95), "severity": "warning"})
    G = pd.concat(tabs, axis=1)
    md = meta_vectorised(G[[f"death_b_{c}" for c in cohorts]].to_numpy(), G[[f"death_se_{c}" for c in cohorts]].to_numpy())
    ma = meta_vectorised(G[[f"age_b_{c}" for c in cohorts]].to_numpy(), G[[f"age_se_{c}" for c in cohorts]].to_numpy())
    from statsmodels.stats.multitest import multipletests
    meta = pd.DataFrame({"gene": G.index, "death_logOR": md["est"].values, "death_se": md["se"].values,
                         "death_z": md["z"].values, "death_P": md["P"].values, "death_k": md["k"].values,
                         "death_I2": md["I2"].values, "age_b": ma["est"].values, "age_z": ma["z"].values,
                         "age_P": ma["P"].values}).set_index("gene")
    meta["death_FDR"] = np.nan
    okp = meta["death_P"].notna()
    meta.loc[okp, "death_FDR"] = multipletests(meta.loc[okp, "death_P"], method="fdr_bh")[1]
    meta.to_csv(C.TAB / "S18_gene_meta.tsv.gz", sep="\t")
    LOG.info(f"Genes associated with death (meta FDR < 0.05): {int((meta['death_FDR'] < 0.05).sum())} "
             f"(higher in non-survivors {int(((meta['death_FDR'] < 0.05) & (meta['death_z'] > 0)).sum())}, "
             f"lower {int(((meta['death_FDR'] < 0.05) & (meta['death_z'] < 0)).sum())})")
    # ---------------- GSEA ----------------
    sets = gene_sets()
    for k_, v_ in sets.items():
        gates.append({"gate": f"{k_} gene sets loaded", "value": len(v_), "criterion": ">= 40" if k_ == "Hallmark" else ">= 500",
                      "passed": len(v_) >= (40 if k_ == "Hallmark" else 500), "severity": "warning"})
    res = {}
    for rank_name, col_ in (("death", "death_z"), ("age", "age_z")):
        eng = PrerankGSEA(meta[col_], A2.N_PERM, A2.SET_MAX, C.SEED)
        for coll, s in sets.items():
            r = eng.run(s)
            r.insert(0, "collection", coll)
            r.insert(0, "ranking", rank_name)
            res[(rank_name, coll)] = r
            LOG.info(f"GSEA {rank_name} / {coll}: {len(r)} sets tested, {int((r['FDR'] < A2.FDR).sum())} at FDR < {A2.FDR}")
    gsea = pd.concat(res.values(), ignore_index=True)
    gsea.to_csv(C.TAB / "S19_gsea.tsv", sep="\t", index=False)
    # ---------------- reversal signature ----------------
    sig_rows = []
    for coll in sets:
        d_, a_ = res[("death", coll)], res[("age", coll)]
        m = d_.merge(a_[["pathway", "NES", "P"]], on="pathway", suffixes=("", "_age"))
        m["concordant"] = (m["FDR"] < A2.FDR) & (np.sign(m["NES"]) == np.sign(m["NES_age"])) & (m["P_age"] < A2.AGE_P)
        m["death_only"] = m["FDR"] < A2.FDR
        sig_rows.append(m.assign(collection=coll))
    pw = pd.concat(sig_rows, ignore_index=True)

    def signature(mask: pd.Series) -> pd.DataFrame:
        genes = {}
        for _, r in pw[mask].iterrows():
            for g in str(r["leading_edge"]).split(";"):
                if g and g in meta.index and meta.loc[g, "death_P"] < A2.SIG_GENE_P and \
                        np.sign(meta.loc[g, "death_z"]) == np.sign(r["NES"]):
                    genes.setdefault(g, []).append(r["pathway"])
        s = pd.DataFrame({"gene": list(genes), "n_pathways": [len(v) for v in genes.values()],
                          "pathways": [" | ".join(v[:5]) for v in genes.values()]})
        if s.empty:
            return s
        s["death_z"] = meta.loc[s["gene"], "death_z"].values
        s["age_z"] = meta.loc[s["gene"], "age_z"].values
        s["direction"] = np.where(s["death_z"] > 0, "UP", "DOWN")
        return s.sort_values("death_z", key=np.abs, ascending=False)

    pre = signature(pw["concordant"])
    n_up, n_dn = int((pre.get("direction", pd.Series(dtype=str)) == "UP").sum()), int((pre.get("direction", pd.Series(dtype=str)) == "DOWN").sum())
    usable = n_up >= A2.MIN_SIG_GENES and n_dn >= A2.MIN_SIG_GENES
    pre.to_csv(C.TAB / "S20_reversal_signature_prespecified.tsv", sep="\t", index=False)
    LOG.info(f"Pre-specified (death-and-age concordant) signature: UP {n_up}, DOWN {n_dn} -> "
             f"{'usable' if usable else 'NOT usable (< 20 genes in a direction)'}")
    gates.append({"gate": "pre-specified reversal signature (death-and-age concordant pathways)",
                  "value": f"UP {n_up} / DOWN {n_dn}", "criterion": f">= {A2.MIN_SIG_GENES} per direction",
                  "passed": bool(usable), "severity": "warning"})
    if not usable:
        exp = signature(pw["death_only"])
        exp.to_csv(C.TAB / "S20b_reversal_signature_death_only_EXPLORATORY.tsv", sep="\t", index=False)
        LOG.info(f"EXPLORATORY death-only signature: UP {int((exp['direction'] == 'UP').sum()) if len(exp) else 0}, "
                 f"DOWN {int((exp['direction'] == 'DOWN').sum()) if len(exp) else 0}")
    pw.drop(columns=["leading_edge"]).to_csv(C.TAB / "S19b_pathway_death_age_concordance.tsv", sep="\t", index=False)
    # ---------------- ligand-receptor axes ----------------
    lrd, covs = {}, {}
    for c in cohorts:
        s_, cv = lr_scores(c, data[c])
        lrd[c] = data[c].join(s_)
        covs[c] = cv
    pd.DataFrame(covs).to_csv(C.TAB / "S21b_lr_gene_coverage.tsv", sep="\t")
    rows = []
    for pair in A2.LR_PAIRS:
        e1, e2, e3, a1, a2 = [], [], [], [], []
        for c in cohorts:
            d_ = lrd[c]
            if pair not in d_.columns or d_[pair].notna().sum() < C.VALIDATION_MIN_N:
                continue
            r = assoc(d_, pair)
            e1.append(r["logratio"])
            e2.append(r["se"])
            r2 = assoc(d_, pair, adjust=("ml_index", "age10", "male"))
            e3.append((r2["logratio"], r2["se"]))
            dd = d_.dropna(subset=[pair, "age10", "male"])
            xc = ["age10", "male"] if dd["male"].nunique() > 1 else ["age10"]
            o = sm.OLS(dd[pair], sm.add_constant(dd[xc])).fit()
            a1.append(o.params["age10"])
            a2.append(o.bse["age10"])
        if len(e1) < 5:                             # axis measurable in fewer than 5 cohorts
            LOG.warning(f"{pair}: measurable in only {len(e1)} cohorts - not analysed")
            continue
        m1 = dl_meta(np.array(e1), np.array(e2))
        m3 = dl_meta(np.array([x[0] for x in e3]), np.array([x[1] for x in e3]))
        ma_ = dl_meta(np.array(a1), np.array(a2))
        rows.append({"axis": pair, "k": m1["k"], "death_OR": np.exp(m1["est"]),
                     "death_lo": np.exp(m1["est"] - 1.96 * m1["se"]), "death_hi": np.exp(m1["est"] + 1.96 * m1["se"]),
                     "death_P": m1["P"], "death_I2": m1["I2"], "death_OR_adjML": np.exp(m3["est"]),
                     "adjML_lo": np.exp(m3["est"] - 1.96 * m3["se"]), "adjML_hi": np.exp(m3["est"] + 1.96 * m3["se"]),
                     "adjML_P": m3["P"], "age_b": ma_["est"], "age_se": ma_["se"], "age_P": ma_["P"]})
    lr = pd.DataFrame(rows)
    gates.append({"gate": "ligand-receptor axes measurable (genes present) in >= 5 cohorts", "value": len(lr),
                  "criterion": f">= 7 of {len(A2.LR_PAIRS)}", "passed": len(lr) >= 7, "severity": "warning"})
    if len(lr):
        for c_ in ("death", "adjML", "age"):
            lr[f"{c_}_FDR"] = multipletests(lr[f"{c_}_P"].fillna(1.0), method="fdr_bh")[1]
        lr.to_csv(C.TAB / "S21_ligand_receptor_bulk.tsv", sep="\t", index=False)
        LOG.info("Ligand-receptor axes (pooled over all cohorts):\n" + lr[["axis", "k", "death_OR", "death_FDR",
                                                                          "death_OR_adjML", "adjML_FDR", "age_b",
                                                                          "age_FDR"]].round(4).to_string(index=False))
    gates.append({"gate": "RESULT (reported): ligand-receptor axes associated with death (FDR < 0.05)",
                  "value": (", ".join(lr.loc[lr["death_FDR"] < 0.05, "axis"]) or "none") if len(lr) else "none",
                  "criterion": "reported", "passed": True, "severity": "warning"})
    gates.append({"gate": "RESULT (reported): pathways with death FDR < 0.05 (Hallmark / Reactome)",
                  "value": f"{int((res[('death', 'Hallmark')]['FDR'] < 0.05).sum())} / "
                           f"{int((res[('death', 'Reactome')]['FDR'] < 0.05).sum()) if len(res[('death', 'Reactome')]) else 0}",
                  "criterion": "reported", "passed": True, "severity": "warning"})
    record_gates("analysis2", gates)
    safe_plot(plot_fig4, res[("death", "Hallmark")], res[("age", "Hallmark")], res[("death", "Reactome")], lr)
    LOG.info(f"Run time {(time.time() - t0) / 60:.1f} min")
    LOG.info(gate_summary().replace("STAGE 1 QUALITY GATES", "ANALYSIS STAGE 2 QUALITY GATES").replace(
        "Stage 2 can start after the field audit has been checked.", "Analysis stage 2 complete.").replace(
        "do not start Stage 2", "check the failed gates"))


# #####################################################################################
# PART D - ANALYSIS STAGE 3: drug repurposing in human PBMCs (OP3 / NeurIPS 2023 official differential expression)
# #####################################################################################
from statsmodels.stats.multitest import multipletests

LOG = get_logger("sepsis_analysis3")
GATES_FILE = C.TAB / "S0_analysis3_quality_gates.tsv"
_GATES.clear()


class A3:
    # official OP3 limma results (pseudobulk per compound x cell type vs DMSO), as released for the Kaggle
    # competition; paths from the OP3 repository (scripts/create_resources/neurips-2023-kaggle.sh); public
    OP3_BASE = "https://openproblems-bio.s3.amazonaws.com/public/neurips-2023-competition/2023-09-14_kaggle_upload/"
    OP3_FILES = ("2023-09-12_de_by_cell_type_train.h5ad", "2023-09-12_de_by_cell_type_test.h5ad")
    STAT_LAYERS = ("sign_log10_pval",)          # sign(logFC) x -log10(P); in X when no layer of that name exists
    CLIP = 4.0                                  # = P 1e-4, the clipping used by the OP3 benchmark
    MIN_DIR_GENES = 20
    FULL_COVERAGE = 0.50
    N_NULL = 10000
    FDR = 0.05
    TAU = 0.10
    NOMINAL = -np.log10(0.05)
    MIN_MHC_GENES_CHANGED = 4                   # MHC-II readout: >= 4 genes individually changed (P < 0.05)
    PRIMARY_CT = r"myeloid|monocyt"
    MHC = list(C.GENE_SETS["mhc_class_ii"])
    ALARMIN = ["S100A8", "S100A9"]
    GLUCOCORTICOIDS = (r"dexamethasone|hydrocortisone|cortisol|prednisolone|prednisone|methylprednisolone|"
                       r"budesonide|mometasone|fluticasone|betamethasone|triamcinolone|clobetasol|beclomethasone|"
                       r"halcinonide|fluocinolone|fluocinonide|desonide|flunisolide|ciclesonide|amcinonide|"
                       r"halobetasol|desoximetasone|fluorometholone|loteprednol|difluprednate|deflazacort")
    GR_TARGETS = ["TSC22D3", "FKBP5", "ZBTB16", "PER1"]
    ERYTHROID = set("RHAG GYPB GYPE GYPA SPTA1 HBQ1 HBD HBB HBA1 HBA2 HBM AHSP CA1 KLF1 GATA1 ALAS2 EPB42 ANK1 "
                    "SLC4A1 TMCC2 XK TRIM10 DMTN BPGM TSPO2 GLRX5 TFRC HMBS RBM38 MPP1 TFDP2 E2F2".split())
    NEUTROPHIL = set("ELANE MMP8 ARG1 ARG2 VNN1 MAOA CEACAM8 LCN2 LTF CAMP OLFM4 DEFA4 MPO CTSG AZU1 BPI".split())
    ENRICHR = "https://maayanlab.cloud/Enrichr"
    LIGAND_LIBS = ("Ligand_Perturbations_from_GEO_up", "Ligand_Perturbations_from_GEO_down")
    LIGAND_CONTROLS = {"IFN-gamma": r"\bIFN-?G(AMMA)?\b|\bIFNG\b|interferon[ -]gamma",
                       "GM-CSF": r"\bGM-?CSF\b|\bCSF2\b"}
    DGIDB = "https://dgidb.org/api/graphql"
    LR_TARGETS = {"HLA-DRA": "increase", "HLA-DRB1": "increase", "HLA-DPA1": "increase", "HLA-DPB1": "increase",
                  "CD4": "increase", "CD74": "increase", "MIF": "increase", "CD86": "increase", "CD28": "increase",
                  "CTLA4": "decrease", "ARG1": "decrease", "S100A8": "decrease", "S100A9": "decrease",
                  "TLR4": "decrease"}
    CONTROL_RECEPTORS = {"CSF2RA": "increase", "CSF2RB": "increase", "IFNGR1": "increase", "IFNGR2": "increase"}
    CONTROL_DRUGS = r"sargramostim|molgramostim|interferon gamma|interferon-gamma|ifn-gamma"


AMENDMENT_5 = """# Amendment 5 - specification of analysis stage 3 (drug repurposing); written before any perturbation data
# were downloaded or opened

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
"""


AMENDMENT_6_HEAD = """# Amendment 6 - change of the stage-3 data source; supersedes Amendment 5
# written before any drug-signature result (reversal, MHC class II, candidates) was computed

Amendment 5 specified the OP3 raw single-cell counts (a file of several GB) with our own pseudobulk model. Because
of the download size, the official OP3 differential-expression files (same experiment; limma pseudobulk models of
the OP3 authors) are used instead. No reversal, MHC class II or candidate result had been computed under Amendment 5;
the script checked this before writing this amendment (audit in AMENDMENT_6.json). The specification below
replaces Amendment 5 in full.

"""


def _sha(t: str) -> str:
    return hashlib.sha256(t.encode()).hexdigest()


def _register(name: str, text: str, extra: dict | None = None) -> None:
    p, meta = C.PREREG / f"{name}.md", C.PREREG / f"{name}.json"
    p.write_text(text)
    meta.write_text(json.dumps({"registered_utc": utc_now(), "sha256": _sha(text), "file": p.name, **(extra or {})},
                               indent=1))
    LOG.info(f"{name.replace('_', ' ').title()} written (SHA-256 {_sha(text)[:16]}...) - before the analyses it specifies")


def _verify(name: str, text: str | None = None) -> dict:
    p, meta = C.PREREG / f"{name}.md", C.PREREG / f"{name}.json"
    old = json.loads(meta.read_text())
    if _sha(p.read_text()) != old["sha256"] or (text is not None and _sha(text) != old["sha256"]):
        raise SystemExit(f"{name}.md differs from its registered version - do not edit it (restore the registered file from the repository).")
    return old


def write_amendment_5() -> None:
    """Amendment 5 (this design) - or, if an earlier Amendment 5 (raw-count design) is registered, Amendment 6."""
    p5 = C.PREREG / "AMENDMENT_5.md"
    results = [C.TAB / "S23_op3_compound_results.tsv", C.TAB / "S24_op3_candidates.tsv"]
    if not p5.exists():
        if any(r.exists() for r in results):
            raise SystemExit("Stage-3 results exist but Amendment 5 is missing - the amendment must precede the analysis.")
        _register("AMENDMENT_5", AMENDMENT_5)
        return
    old5 = _verify("AMENDMENT_5")
    if old5["sha256"] == _sha(AMENDMENT_5):
        LOG.info(f"Amendment 5 present (registered {old5['registered_utc']}, unchanged)")
        return
    if "Pseudobulk = sum of raw counts" not in p5.read_text():
        raise SystemExit("AMENDMENT_5.md belongs to another version of this stage - send it to me before going on.")
    LOG.info(f"Amendment 5 (raw-count design) registered {old5['registered_utc']} - superseded by Amendment 6")
    text6 = AMENDMENT_6_HEAD + AMENDMENT_5.replace(
        "# Amendment 5 - specification of analysis stage 3 (drug repurposing); written before any perturbation data\n"
        "# were downloaded or opened\n", "## Specification of analysis stage 3 (drug repurposing)\n")
    if (C.PREREG / "AMENDMENT_6.md").exists():
        old6 = _verify("AMENDMENT_6", text6)
        LOG.info(f"Amendment 6 present (registered {old6['registered_utc']}, unchanged)")
        return
    raw = C.RAW / "op3" / "sc_counts_reannotated_with_counts.h5ad"
    part = raw.with_name(raw.name + ".part")
    audit = {"amendment_5_registered_utc": old5["registered_utc"],
             "drug_results_present": [r.name for r in results if r.exists()],
             "raw_count_file_present": raw.exists(), "raw_count_partial_download_bytes":
                 part.stat().st_size if part.exists() else 0,
             "raw_pseudobulk_cache_present": (C.PROC / "op3_pseudobulk.npz").exists()}
    if audit["drug_results_present"]:
        raise SystemExit("Drug results of the Amendment 5 design already exist - check them before continuing; the design must "
                         "not be changed silently.")
    _register("AMENDMENT_6", text6, {"audit": audit})
    LOG.info(f"Audit before Amendment 6: {audit}")
    for f_ in (part, raw):
        if f_.exists():
            LOG.info(f"Not needed any more (can be deleted to free {f_.stat().st_size / 1e9:.1f} GB): {f_}")


# -------------------------------------------------------------------------------------
# OP3 differential-expression profiles
# -------------------------------------------------------------------------------------
def h5_ok(p: Path) -> bool:
    try:
        import h5py
        with h5py.File(p, "r") as f:
            return "obs" in f
    except Exception:
        return False


def download_resumable(url: str, dest: Path, attempts: int = 60) -> Path:
    """Download with HTTP Range: after a broken connection the next attempt continues where the last one stopped.
    The partial file (<name>.part) is kept between runs, so re-running the cell also continues."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    total = None
    for k in range(attempts):
        have = part.stat().st_size if part.exists() else 0
        if total is not None and have >= total:
            break
        try:
            with HTTP.get(url, headers={"Range": f"bytes={have}-"} if have else {}, stream=True,
                          timeout=(30, 120)) as r:
                if r.status_code == 416:                       # nothing left to fetch
                    break
                if r.status_code == 404:
                    raise SystemExit(f"Not found: {url}")
                r.raise_for_status()
                if r.status_code == 206 and "Content-Range" in r.headers:
                    total = int(r.headers["Content-Range"].split("/")[-1])
                    mode = "ab"
                else:                                          # server sent the whole file
                    total = int(r.headers.get("Content-Length", 0)) or None
                    mode, have = "wb", 0
                if k == 0 or have == 0:
                    LOG.info(f"Downloading {dest.name}" + (f" ({total / 1e6:.0f} MB)" if total else ""))
                done, nxt = have, have + 50e6
                with open(part, mode) as fh:
                    for chunk in r.iter_content(1 << 20):
                        fh.write(chunk)
                        done += len(chunk)
                        if done >= nxt:
                            LOG.info(f"  {dest.name}: {done / 1e6:.0f}" + (f" / {total / 1e6:.0f} MB" if total else " MB"))
                            nxt += 50e6
            if total is None:
                break
        except requests.RequestException as exc:
            got = part.stat().st_size if part.exists() else 0
            LOG.warning(f"  connection interrupted at {got / 1e6:.0f} MB ({type(exc).__name__}); resuming "
                        f"(attempt {k + 2}/{attempts})")
            time.sleep(min(30, 3 * (k + 1)))
    got = part.stat().st_size if part.exists() else 0
    if total is None or got < total:
        raise SystemExit(f"{dest.name}: download incomplete ({got / 1e6:.0f} MB). Run the cell again - it continues "
                         f"from {got / 1e6:.0f} MB. Or download by hand (resumable):\n  wget -c -O {dest} \"{url}\"")
    part.replace(dest)
    return dest


def op3_files() -> list:
    d = C.RAW / "op3"
    out = []
    for name in A3.OP3_FILES:
        dest = d / name
        if dest.exists() and not h5_ok(dest):
            LOG.warning(f"{dest.name} is incomplete or not an h5ad file - it is moved aside and downloaded again")
            dest.replace(dest.with_name(dest.name + ".broken"))
        if not dest.exists():
            download_resumable(A3.OP3_BASE + name, dest)
            if not h5_ok(dest):
                dest.replace(dest.with_name(dest.name + ".broken"))
                raise SystemExit(f"{name} was downloaded but is not a readable h5ad file - run the cell again.")
            LOG.info(f"OP3 file ready: {dest.name} ({dest.stat().st_size / 1e6:.0f} MB)")
        out.append(dest)
    return out


def op3_profiles(files: list) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame | None, str]:
    """Rows = compound x cell type, columns = genes. Returns statistic, metadata, logFC (if present), name."""
    import anndata as ad
    stats_, metas, lfcs, used = [], [], [], set()
    for p in files:
        a = ad.read_h5ad(p)
        lay = next((l_ for l_ in A3.STAT_LAYERS if l_ in a.layers), None)
        X = a.layers[lay] if lay else a.X
        if X is None:
            raise SystemExit(f"{p.name}: no sign_log10_pval layer and no X matrix")
        used.add(lay or "X (sign_log10_pval)")
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        stats_.append(pd.DataFrame(X.astype(float), columns=a.var_names.astype(str)))
        lfcs.append(pd.DataFrame(np.asarray(a.layers["logFC"].toarray() if hasattr(a.layers["logFC"], "toarray")
                                            else a.layers["logFC"], float), columns=a.var_names.astype(str))
                    if "logFC" in a.layers else None)
        for col_ in ("cell_type", "sm_name"):
            if col_ not in a.obs.columns:
                raise SystemExit(f"{p.name}: obs column '{col_}' missing (columns: {list(a.obs.columns)})")
        m = a.obs[["cell_type", "sm_name"]].astype(str).reset_index(drop=True)
        m["split"] = a.obs["split"].astype(str).values if "split" in a.obs.columns else ""
        m["file"] = p.name
        metas.append(m)
    T = pd.concat(stats_, ignore_index=True)
    M = pd.concat(metas, ignore_index=True)
    L = pd.concat(lfcs, ignore_index=True) if all(x is not None for x in lfcs) else None
    dup = M.duplicated(subset=["cell_type", "sm_name"]).to_numpy()
    if dup.any():
        LOG.warning(f"{int(dup.sum())} duplicated compound x cell-type profiles removed")
    keep = ~dup
    T, M = T[keep].reset_index(drop=True), M[keep].reset_index(drop=True)
    L = L[keep].reset_index(drop=True) if L is not None else None
    if T.isna().all(axis=0).any():
        T = T.loc[:, ~T.isna().all(axis=0)]
    return T, M, L, ", ".join(sorted(used))


# -------------------------------------------------------------------------------------
# connectivity and gene-set readouts against random gene sets
# -------------------------------------------------------------------------------------
def random_sets(pool: np.ndarray, k: int, n: int, rng) -> np.ndarray:
    if k <= len(pool):
        return pool[np.argpartition(rng.random((n, len(pool))), k - 1, axis=1)[:, :k]]
    return rng.choice(pool, (n, k), replace=True)


def null_scores(M: np.ndarray, du: np.ndarray, dd: np.ndarray | None) -> np.ndarray:
    null = np.empty((M.shape[0], du.shape[0]))
    for a in range(0, du.shape[0], 250):
        b = min(a + 250, du.shape[0])
        null[:, a:b] = M[:, du[a:b]].mean(axis=2) - (M[:, dd[a:b]].mean(axis=2) if dd is not None else 0)
    return null


def summarise(obs: np.ndarray, null: np.ndarray) -> pd.DataFrame:
    mu, sd = null.mean(axis=1), null.std(axis=1, ddof=1)
    z = (obs - mu) / sd
    pe = (1 + (np.abs(null - mu[:, None]) >= np.abs(obs - mu)[:, None]).sum(axis=1)) / (null.shape[1] + 1)
    return pd.DataFrame({"score": obs, "z": z, "P": 2 * stats.norm.sf(np.abs(z)), "P_empirical": pe})


def connectivity(T: pd.DataFrame, up: list, dn: list, seed: int) -> pd.DataFrame:
    genes = list(T.columns)
    pos = {g: i for i, g in enumerate(genes)}
    u, d = [pos[g] for g in up if g in pos], [pos[g] for g in dn if g in pos]
    pool = np.array(sorted(set(range(len(genes))) - set(u) - set(d)))
    rng = np.random.default_rng(seed)
    du, dd = random_sets(pool, len(u), A3.N_NULL, rng), random_sets(pool, len(d), A3.N_NULL, rng)
    M = np.clip(np.nan_to_num(T.to_numpy(float)), -A3.CLIP, A3.CLIP)
    out = summarise(M[:, u].mean(axis=1) - M[:, d].mean(axis=1), null_scores(M, du, dd))
    return out.rename(columns={"score": "R", "z": "R_z", "P": "R_P", "P_empirical": "R_P_empirical"})


def set_shift(T: pd.DataFrame, genes: list, seed: int, prefix: str) -> pd.DataFrame:
    cols = list(T.columns)
    pos = {g: i for i, g in enumerate(cols)}
    g = [pos[x] for x in genes if x in pos]
    M = np.clip(np.nan_to_num(T.to_numpy(float)), -A3.CLIP, A3.CLIP)
    if len(g) < 2:
        return pd.DataFrame(index=range(len(T)))
    pool = np.array(sorted(set(range(len(cols))) - set(g)))
    du = random_sets(pool, len(g), A3.N_NULL, np.random.default_rng(seed))
    out = summarise(M[:, g].mean(axis=1), null_scores(M, du, None))
    raw = np.nan_to_num(T.to_numpy(float))[:, g]
    out["n_up"] = (raw > A3.NOMINAL).sum(axis=1)
    out["n_down"] = (raw < -A3.NOMINAL).sum(axis=1)
    out["n_genes"] = len(g)
    return out.add_prefix(prefix + "_")


# -------------------------------------------------------------------------------------
# Enrichr ligand perturbations (IFN-gamma, GM-CSF) and DGIdb target layer
# -------------------------------------------------------------------------------------
def enrichr(genes: list, library: str) -> list:
    r = HTTP.post(A3.ENRICHR + "/addList", files={"list": (None, "\n".join(genes)), "description": (None, "sepsis")},
                  timeout=120)
    r.raise_for_status()
    uid = r.json()["userListId"]
    r = http_get(f"{A3.ENRICHR}/enrich?userListId={uid}&backgroundType={library}", timeout=180)
    return r.json().get(library, [])


def ligand_controls(up: list, dn: list) -> pd.DataFrame:
    rows = []
    for lib in A3.LIGAND_LIBS:
        lig_dir = "up" if lib.endswith("_up") else "down"
        for sig_dir, genes in (("UP", up), ("DOWN", dn)):
            cache = C.RAW / "enrichr" / f"{lib}_{sig_dir}.json"
            if cache.exists():
                res = json.loads(cache.read_text())
            else:
                res = enrichr(genes, lib)
                cache.parent.mkdir(parents=True, exist_ok=True)
                cache.write_text(json.dumps(res))
            for row in res:
                term = str(row[1])
                for ctrl, rx in A3.LIGAND_CONTROLS.items():
                    if re.search(rx, term, re.I):
                        rev = (sig_dir == "DOWN" and lig_dir == "up") or (sig_dir == "UP" and lig_dir == "down")
                        rows.append({"control": ctrl, "library": lib, "signature_genes": sig_dir, "term": term,
                                     "P": float(row[2]), "adj_P": float(row[6]), "overlap": ";".join(row[5]),
                                     "relation": "reversing" if rev else "mimicking"})
    return pd.DataFrame(rows, columns=["control", "library", "signature_genes", "term", "P", "adj_P", "overlap",
                                       "relation"])


def dgidb(genes: list) -> pd.DataFrame:
    q = ("query($n:[String!]){genes(names:$n){nodes{name interactions{interactionScore "
         "interactionTypes{type directionality} drug{name conceptId approved}}}}}")
    cache = C.RAW / "dgidb" / ("targets_" + hashlib.sha256(",".join(genes).encode()).hexdigest()[:10] + ".json")
    if cache.exists():
        res = json.loads(cache.read_text())
    else:
        r = HTTP.post(A3.DGIDB, json={"query": q, "variables": {"n": genes}}, timeout=120)
        r.raise_for_status()
        res = r.json()
        if res.get("errors") or not res.get("data"):
            raise RuntimeError(str(res.get("errors"))[:300])
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(res))
    want_all = {**A3.LR_TARGETS, **A3.CONTROL_RECEPTORS}
    rows = []
    for nd in res["data"]["genes"]["nodes"]:
        for it in nd.get("interactions") or []:
            types = [t.get("type") for t in it.get("interactionTypes") or [] if t.get("type")]
            dirs = [str(t.get("directionality")).upper() for t in it.get("interactionTypes") or []
                    if t.get("directionality")]
            want = want_all.get(nd["name"], "")
            inh = "INHIBITORY" in dirs or any(re.search(r"inhibit|antagon|block|antibod|negative", t, re.I) for t in types)
            act = "ACTIVATING" in dirs or any(re.search(r"agonist|activat|induc|stimul|positive", t, re.I) for t in types)
            rows.append({"gene": nd["name"], "role": "positive control" if nd["name"] in A3.CONTROL_RECEPTORS
                         else "death-associated axis", "desired": want, "drug": it["drug"]["name"],
                         "approved": bool(it["drug"].get("approved")), "interaction_types": ";".join(types),
                         "directionality": ";".join(dirs),
                         "matches_desired_direction": bool((want == "decrease" and inh) or (want == "increase" and act)),
                         "score": it.get("interactionScore")})
    return pd.DataFrame(rows, columns=["gene", "role", "desired", "drug", "approved", "interaction_types",
                                       "directionality", "matches_desired_direction", "score"])




# -------------------------------------------------------------------------------------
# figure 5
# -------------------------------------------------------------------------------------
def short_ct(c: str) -> str:
    return (c.replace("T regulatory cells", "Treg").replace("T cells CD4+", "CD4 T").replace("T cells CD8+", "CD8 T")
            .replace(" cells", ""))


def plot_fig5(genes_tab: pd.DataFrame, res: pd.DataFrame, cand: pd.DataFrame, pc: pd.DataFrame, dg: pd.DataFrame,
              primary: str) -> None:
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    set_style()
    fig = plt.figure(figsize=(FIG_W, 140 * MM))
    outer = fig.add_gridspec(2, 1, height_ratios=[1, 1.1], hspace=0.55, left=0.08, right=0.98, top=0.95, bottom=0.08)
    top_ = outer[0].subgridspec(1, 2, width_ratios=[1, 1.9], wspace=0.35)
    bot_ = outer[1].subgridspec(1, 3, width_ratios=[1.05, 1, 1], wspace=1.25)
    # a: which signature genes are testable in PBMCs
    ax = fig.add_subplot(top_[0, 0])
    cats = [("in OP3 (testable)", OKABE[2]), ("erythroid (absent)", "#E8A0A0"), ("neutrophil (absent)", "#F2C879"),
            ("other (absent)", "#CCCCCC")]
    for i, dirn in enumerate(("UP", "DOWN")):
        s = genes_tab[genes_tab["direction"] == dirn]
        n = len(s)
        vals = [s["in_OP3"].sum(), ((~s["in_OP3"]) & (s["lineage"] == "erythroid")).sum(),
                ((~s["in_OP3"]) & (s["lineage"] == "neutrophil")).sum(), ((~s["in_OP3"]) & (s["lineage"] == "")).sum()]
        left = 0
        for v, (_, col) in zip(vals, cats):
            ax.barh(1 - i, v / n, left=left, color=col, height=0.55)
            left += v / n
        ax.text(1.02, 1 - i, f"{vals[0]}/{n}", va="center", fontsize=5)
    ax.axvline(A3.FULL_COVERAGE, color="#555555", lw=0.5, ls="--")
    ax.set_yticks([1, 0])
    ax.set_yticklabels(["UP in non-survivors", "DOWN in non-survivors"])
    ax.set_xlim(0, 1.18)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xlabel("Fraction of signature genes")
    ax.legend(handles=[Patch(color=c_, label=l_) for l_, c_ in cats], frameon=False, fontsize=4.5, ncol=2,
              loc="upper center", bbox_to_anchor=(0.45, -0.25))
    ax.set_title("Testable part of the death signature\n(OP3 PBMC gene universe)", fontsize=6)
    panel_label(ax, "a", x=-0.45)
    # b: primary cell type, reversal vs MHC-II
    ax = fig.add_subplot(top_[0, 1])
    p = res[res["cell_type"] == primary].copy()
    if len(p) and "mhc_z" in p:
        has_r = p["R_z"].notna().any()
        xv = p["R_z"] if has_r else pd.Series(0.0, index=p.index)
        tier = p["sm_name"].map(cand.set_index("compound")["tier"]) if len(cand) else pd.Series(np.nan, index=p.index)
        colr = np.array(["#9E9E9E"] * len(p), dtype=object)
        colr[(tier == "Tier 2").to_numpy()] = OKABE[4]
        colr[(tier == "Tier 1").to_numpy()] = OKABE[2]
        colr[(tier == "Caution").to_numpy()] = OKABE[1]
        ax.scatter(xv, p["mhc_z"], s=9, c=list(colr), edgecolor="white", linewidth=0.3, zorder=3)
        ax.axhline(0, color="#777777", lw=0.5)
        ax.axvline(0, color="#777777", lw=0.5)
        lab = p[tier.isin(["Tier 1", "Tier 2", "Caution"]).to_numpy() |
                p["sm_name"].str.contains(A3.GLUCOCORTICOIDS, case=False, regex=True).to_numpy()]
        order = (lab["mhc_z"].abs() + xv.loc[lab.index].abs()).sort_values(ascending=False).index
        for i_ in order[:14]:
            ax.annotate(str(p.loc[i_, "sm_name"])[:22], (xv.loc[i_], p.loc[i_, "mhc_z"]), fontsize=4.5,
                        xytext=(3, 2), textcoords="offset points")
        ax.set_xlabel("Signature reversal z (< 0: reverses)" if has_r else "Signature reversal not testable (gate)")
        ax.set_ylabel("MHC class II shift z (> 0: restores)")
        hs = [Line2D([], [], marker="o", ls="", color=c_, markersize=3, label=l_) for c_, l_ in
              ((OKABE[2], "Tier 1"), (OKABE[4], "Tier 2"), (OKABE[1], "Caution"), ("#9E9E9E", "Other"))]
        ax.legend(handles=hs, frameon=False, fontsize=4.5, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    ax.set_title(f"{primary}: {len(p)} compounds (OP3, 24 h)", fontsize=6)
    panel_label(ax, "b", x=-0.12)
    # c: candidates x cell types
    ax = fig.add_subplot(bot_[0, 0])
    cts = sorted(res["cell_type"].unique())
    if res["R_z"].notna().any():
        top = (cand.loc[cand["tier"] != "Caution", "compound"].head(15).tolist() if len(cand) else []) or \
            res[res["cell_type"] == primary].nsmallest(15, "R_z")["sm_name"].tolist()
        Mz = res.pivot_table(index="sm_name", columns="cell_type", values="R_z").reindex(index=top, columns=cts)
        Fz = res.pivot_table(index="sm_name", columns="cell_type", values="R_FDR").reindex(index=top, columns=cts)
        lim = np.nanmax(np.abs(Mz.to_numpy())) if np.isfinite(Mz.to_numpy()).any() else 1
        im = ax.imshow(Mz.to_numpy(), cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto", interpolation="none")
        for (i, j), v in np.ndenumerate(Fz.to_numpy()):
            if np.isfinite(v) and v < A3.FDR:
                ax.text(j, i, "*", ha="center", va="center", fontsize=6)
        ax.set_xticks(range(len(cts)))
        ax.set_xticklabels([short_ct(c) for c in cts], rotation=45, ha="right")
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels([t[:24] for t in top], fontsize=4.5)
        cb = fig.colorbar(im, ax=ax, fraction=0.06, pad=0.04, shrink=0.8)
        cb.set_label("Reversal z", fontsize=5)
        cb.ax.tick_params(labelsize=4.5)
    else:
        ax.text(0.5, 0.5, "Reversal not testable\n(expression gate)", ha="center", va="center", fontsize=6)
        ax.axis("off")
    ax.set_title("Reversal across PBMC cell types\n(* FDR < 0.05)", fontsize=6)
    panel_label(ax, "c", x=-0.75)
    # d: positive controls
    ax = fig.add_subplot(bot_[0, 1])
    if len(pc):
        yy = np.arange(len(pc))[::-1]
        ax.barh(yy, pc["value"], color=[OKABE[2] if ok else "#BBBBBB" for ok in pc["passed"]], height=0.6)
        ax.set_yticks(yy)
        ax.set_yticklabels([s[:34] for s in pc["label"]], fontsize=4.5)
        ax.axvline(0, color="#777777", lw=0.5)
    else:
        ax.text(0.5, 0.5, "No control testable", ha="center", va="center", fontsize=6)
    ax.set_xlabel("Effect (n genes, z or -log10 adj. P; S25)")
    ax.set_title("Positive controls\n(green: as expected)", fontsize=6)
    panel_label(ax, "d", x=-1.05)
    # e: DGIdb
    ax = fig.add_subplot(bot_[0, 2])
    if len(dg):
        t = dg[dg["approved"]].groupby("gene").agg(n=("drug", "nunique"), match=("matches_desired_direction", "sum"))
        order = [g for g in list(A3.LR_TARGETS) + list(A3.CONTROL_RECEPTORS) if g in t.index]
        t = t.reindex(order)
        yy = np.arange(len(t))[::-1]
        ax.barh(yy, t["n"], color="#CCCCCC", height=0.6, label="approved drugs")
        ax.barh(yy, t["match"], color=OKABE[0], height=0.6, label="desired direction")
        ax.set_yticks(yy)
        ax.set_yticklabels([f"{g} ({(A3.LR_TARGETS | A3.CONTROL_RECEPTORS)[g][:3]}.)" for g in t.index], fontsize=4.5)
        ax.legend(frameon=False, fontsize=4.5, loc="upper center", bbox_to_anchor=(0.4, -0.16), ncol=2)
    else:
        ax.text(0.5, 0.5, "DGIdb not reachable", ha="center", va="center", fontsize=6)
    ax.set_xlabel("Approved drugs (DGIdb)")
    ax.set_title("Druggable axis genes", fontsize=6)
    panel_label(ax, "e", x=-0.75)
    save_figure(fig, "Fig5_drug_repurposing")


# -------------------------------------------------------------------------------------
# run
# -------------------------------------------------------------------------------------
def preflight3() -> None:
    need = {C.TAB / "S20_reversal_signature_prespecified.tsv": "03_analysis_pathways_LR.py",
            C.TAB / "S21_ligand_receptor_bulk.tsv": "03_analysis_pathways_LR.py",
            C.PREREG / "AMENDMENT_4.md": "03_analysis_pathways_LR.py",
            C.TAB / "S17_primary_results.tsv": "02_analysis_H1_H2.py"}
    miss = [f"{p.name} (made by {w})" for p, w in need.items() if not p.exists()]
    if miss:
        raise SystemExit(f"Project folder {C.PROJECT_DIR} is missing: {'; '.join(miss)}.")
    try:
        import anndata
        import h5py
        LOG.info(f"anndata {getattr(anndata, '__version__', '?')}, h5py {h5py.__version__}")
    except ImportError as exc:
        raise SystemExit(f"Package missing ({exc}); install with: pip install anndata h5py")


def run_analysis3() -> None:
    t0 = time.time()
    warnings.filterwarnings("ignore", category=FutureWarning)
    LOG.info(f"Project folder: {C.PROJECT_DIR}")
    preflight()
    preflight3()
    write_preregistration()
    write_amendment()
    write_amendment_3()
    write_amendment_4()
    write_amendment_5()
    sig = pd.read_csv(C.TAB / "S20_reversal_signature_prespecified.tsv", sep="\t")
    up = sig.loc[sig["direction"] == "UP", "gene"].tolist()
    dn = sig.loc[sig["direction"] == "DOWN", "gene"].tolist()
    if min(len(up), len(dn)) < A3.MIN_DIR_GENES:
        raise SystemExit(f"S20 signature has UP {len(up)} / DOWN {len(dn)} genes (< {A3.MIN_DIR_GENES}) - re-run 03.")
    LOG.info(f"Reversal signature: UP {len(up)}, DOWN {len(dn)} genes "
             f"(UP erythroid {sum(g in A3.ERYTHROID for g in up)}, neutrophil {sum(g in A3.NEUTROPHIL for g in up)})")
    gates = []
    # ---------------- OP3 profiles ----------------
    T, M, LFC, stat = op3_profiles(op3_files())
    dmso = M["sm_name"].str.contains(r"^(?:dimethyl sulfoxide|dmso)$", case=False, regex=True)
    if dmso.any():
        T, M = T[~dmso.values].reset_index(drop=True), M[~dmso.values].reset_index(drop=True)
        LFC = LFC[~dmso.values].reset_index(drop=True) if LFC is not None else None
    cts = sorted(M["cell_type"].unique())
    LOG.info(f"OP3: {len(M)} compound x cell-type profiles, {M['sm_name'].nunique()} compounds, {T.shape[1]} genes; "
             f"cell types {cts}; statistic {stat}")
    prim = [c for c in cts if re.search(A3.PRIMARY_CT, c, re.I)]
    primary = prim[0] if prim else None
    gates.append({"gate": "OP3 profiles: compounds / primary cell type", "value": f"{M['sm_name'].nunique()} / "
                  f"{primary or 'none'}", "criterion": ">= 100 compounds; myeloid cells present",
                  "passed": bool(M["sm_name"].nunique() >= 100 and primary is not None), "severity": "critical"})
    n_prim = int((M["cell_type"] == primary).sum()) if primary else 0
    gates.append({"gate": f"OP3: compounds profiled in {primary}", "value": n_prim, "criterion": ">= 100",
                  "passed": n_prim >= 100, "severity": "critical"})
    if not gates[-1]["passed"] or not gates[-2]["passed"]:
        record_gates("analysis3", gates)
        raise SystemExit("OP3 files do not have the expected content (both the train and the test file are needed).")
    # ---------------- expression gate ----------------
    universe = set(T.columns)
    gtab = pd.DataFrame({"gene": up + dn, "direction": ["UP"] * len(up) + ["DOWN"] * len(dn)})
    gtab["in_OP3"] = gtab["gene"].isin(universe)
    gtab["lineage"] = np.where(gtab["gene"].isin(A3.ERYTHROID), "erythroid",
                               np.where(gtab["gene"].isin(A3.NEUTROPHIL), "neutrophil", ""))
    gtab.to_csv(C.TAB / "S22b_op3_signature_genes.tsv", sep="\t", index=False)
    cov = gtab.groupby("direction").agg(n_genes=("gene", "size"), in_OP3=("in_OP3", "sum")).reset_index()
    cov["fraction"] = cov["in_OP3"] / cov["n_genes"]
    cov["coverage"] = np.where(cov["fraction"] >= A3.FULL_COVERAGE, "full", "partial")
    miss = gtab[~gtab["in_OP3"]].groupby(["direction", "lineage"]).size().rename("n_absent").reset_index()
    cov.to_csv(C.TAB / "S22_op3_signature_coverage.tsv", sep="\t", index=False)
    LOG.info("Signature genes in the OP3 PBMC gene universe:\n" + cov.to_string(index=False) +
             "\nAbsent genes by lineage:\n" + miss.replace("", "other").to_string(index=False))
    cu, cd = (cov.set_index("direction").loc[d_, "in_OP3"] for d_ in ("UP", "DOWN"))
    gate_ok = bool(min(cu, cd) >= A3.MIN_DIR_GENES)
    fr = cov.set_index("direction")["fraction"]
    gates.append({"gate": "EXPRESSION GATE: signature genes in the OP3 PBMC gene universe (UP / DOWN)",
                  "value": f"{int(cu)} ({fr['UP']:.0%}) / {int(cd)} ({fr['DOWN']:.0%})",
                  "criterion": f">= {A3.MIN_DIR_GENES} per direction (>= 50% = full coverage)", "passed": gate_ok,
                  "severity": "warning"})
    # ---------------- readouts ----------------
    parts = [M]
    if gate_ok:
        parts.append(connectivity(T, up, dn, C.SEED))
    else:
        LOG.warning("Fewer than 20 signature genes in a direction - reversal not tested (pre-specified); "
                    "the MHC class II readout and the target layer are still reported")
    parts.append(set_shift(T, A3.MHC, C.SEED + 1, "mhc"))
    parts.append(set_shift(T, A3.ALARMIN, C.SEED + 2, "alarmin"))
    res = pd.concat(parts, axis=1)
    for c_ in ("R", "R_z", "R_P", "R_P_empirical", "mhc_z", "mhc_P", "mhc_n_up", "mhc_n_down", "alarmin_z", "alarmin_P"):
        if c_ not in res:
            res[c_] = np.nan
    for c_ in ("R", "mhc", "alarmin"):
        pcol, fcol = (f"{c_}_P", f"{c_}_FDR") if c_ != "R" else ("R_P", "R_FDR")
        res[fcol] = np.nan
        for _, idx in res.groupby("cell_type").groups.items():
            ok = res.loc[idx, pcol].notna()
            if ok.any():
                res.loc[ok[ok].index, fcol] = multipletests(res.loc[ok[ok].index, pcol], method="fdr_bh")[1]
    res["R_tau"] = res.groupby("cell_type")["R"].rank(pct=True)
    res["reverses"] = (res["R"] < 0) & (res["R_FDR"] < A3.FDR) & (res["R_tau"] <= A3.TAU)
    res["mimics"] = (res["R"] > 0) & (res["R_FDR"] < A3.FDR) & (res["R_tau"] >= 1 - A3.TAU)
    res["mhc_up"] = (res["mhc_z"] > 0) & (res["mhc_FDR"] < A3.FDR) & (res["mhc_n_up"] >= A3.MIN_MHC_GENES_CHANGED)
    res["mhc_down"] = (res["mhc_z"] < 0) & (res["mhc_FDR"] < A3.FDR) & (res["mhc_n_down"] >= A3.MIN_MHC_GENES_CHANGED)
    res["alarmin_down"] = (res["alarmin_z"] < 0) & (res["alarmin_FDR"] < A3.FDR)
    res.to_csv(C.TAB / "S23_op3_compound_results.tsv", sep="\t", index=False)
    LOG.info("Per cell type: compounds reversing / mimicking the signature, raising / lowering MHC class II:\n" +
             res.groupby("cell_type")[["reverses", "mimics", "mhc_up", "mhc_down"]].sum().to_string())
    # ---------------- candidates ----------------
    pr = res[res["cell_type"] == primary].set_index("sm_name")
    n_rev = res[res["reverses"]].groupby("sm_name")["cell_type"].apply(lambda v: ", ".join(sorted(v)))
    rows = []
    for c in sorted(res["sm_name"].unique()):
        p_ = pr.loc[c] if c in pr.index else None
        rev_p = bool(p_ is not None and p_["reverses"])
        mhc_p = bool(p_ is not None and p_["mhc_up"])
        n_ct = len(n_rev[c].split(", ")) if c in n_rev.index else 0
        if rev_p and mhc_p:
            tier = "Tier 1"
        elif rev_p or mhc_p or n_ct >= 2:
            tier = "Tier 2"
        elif p_ is not None and (p_["mimics"] or p_["mhc_down"]):
            tier = "Caution"
        else:
            continue
        g_ = (lambda k: p_[k] if p_ is not None else np.nan)
        rows.append({"compound": c, "tier": tier, "reversal_primary": rev_p, "mhc_up_primary": mhc_p,
                     "reversed_in": n_rev.get(c, ""), "n_cell_types_reversed": n_ct, "R_z_primary": g_("R_z"),
                     "R_FDR_primary": g_("R_FDR"), "mhc_z_primary": g_("mhc_z"), "mhc_FDR_primary": g_("mhc_FDR"),
                     "mhc_genes_up_primary": g_("mhc_n_up"),
                     "alarmin_down_primary": bool(p_ is not None and p_["alarmin_down"]),
                     "mimics_primary": bool(p_ is not None and p_["mimics"]),
                     "mhc_down_primary": bool(p_ is not None and p_["mhc_down"]),
                     "glucocorticoid": bool(re.search(A3.GLUCOCORTICOIDS, c, re.I))})
    cand = pd.DataFrame(rows, columns=["compound", "tier", "reversal_primary", "mhc_up_primary", "reversed_in",
                                       "n_cell_types_reversed", "R_z_primary", "R_FDR_primary", "mhc_z_primary",
                                       "mhc_FDR_primary", "mhc_genes_up_primary", "alarmin_down_primary",
                                       "mimics_primary", "mhc_down_primary", "glucocorticoid"])
    cand["_o"] = cand["tier"].map({"Tier 1": 0, "Tier 2": 1, "Caution": 2})
    cand = cand.sort_values(["_o", "R_z_primary", "mhc_z_primary"], ascending=[True, True, False]).drop(columns="_o")
    cand.to_csv(C.TAB / "S24_op3_candidates.tsv", sep="\t", index=False)
    LOG.info(f"Candidates ({primary}): Tier 1 {int((cand['tier'] == 'Tier 1').sum())}, Tier 2 "
             f"{int((cand['tier'] == 'Tier 2').sum())}, caution {int((cand['tier'] == 'Caution').sum())}\n"
             + cand.head(30).round(3).to_string(index=False))
    # ---------------- positive controls ----------------
    pcs = []
    gcs = sorted(c for c in res["sm_name"].unique() if re.search(A3.GLUCOCORTICOIDS, c, re.I))
    LOG.info(f"Glucocorticoids in OP3: {gcs or 'none'}")
    gr_pass, gc_mhc_pass = False, False
    Tn = np.nan_to_num(T.to_numpy(float))
    for g_ in gcs:
        for i_ in M.index[M["sm_name"] == g_]:
            s_ = pd.Series(Tn[i_], index=T.columns)
            fdr = pd.Series(multipletests(10 ** (-np.abs(s_.to_numpy())), method="fdr_bh")[1], index=T.columns)
            hits = [x for x in A3.GR_TARGETS if x in s_.index and s_[x] > 0 and fdr[x] < A3.FDR]
            gr_pass = gr_pass or len(hits) >= 2
            ct = M.loc[i_, "cell_type"]
            pcs.append({"control": "(a) GR targets induced", "compound": g_, "cell_type": ct,
                        "label": f"{g_[:18]} GR targets ({short_ct(ct)})", "value": len(hits),
                        "detail": ", ".join(f"{x}: {s_[x]:+.1f}" for x in A3.GR_TARGETS if x in s_.index),
                        "passed": len(hits) >= 2})
        if g_ in pr.index and pd.notna(pr.loc[g_, "mhc_z"]):
            ok_ = bool(pr.loc[g_, "mhc_z"] < 0 and pr.loc[g_, "mhc_P"] < 0.05)
            gc_mhc_pass = gc_mhc_pass or ok_
            pcs.append({"control": "(b) GC lower myeloid MHC-II", "compound": g_, "cell_type": primary,
                        "label": f"{g_[:18]} MHC-II z (myeloid)", "value": float(pr.loc[g_, "mhc_z"]),
                        "detail": f"P = {pr.loc[g_, 'mhc_P']:.2g}", "passed": ok_})
    gates.append({"gate": "POSITIVE CONTROL (a): glucocorticoids induce GR targets (TSC22D3, FKBP5, ZBTB16, PER1)",
                  "value": (f"{'passed' if gr_pass else 'not passed'} ({', '.join(gcs)})" if gcs else
                            "no glucocorticoid in OP3 (not testable)"),
                  "criterion": ">= 2 targets FDR < 0.05 in >= 1 cell type", "passed": bool(gr_pass), "severity": "warning"})
    gates.append({"gate": f"POSITIVE CONTROL (b): glucocorticoids lower MHC class II in {primary}",
                  "value": "passed" if gc_mhc_pass else ("not passed" if gcs else "not testable"),
                  "criterion": "z < 0, P < 0.05", "passed": bool(gc_mhc_pass), "severity": "warning"})
    try:
        lig = ligand_controls(up, dn)
        lig.to_csv(C.TAB / "S25b_positive_control_cytokines_enrichr.tsv", sep="\t", index=False)
        for ctrl in A3.LIGAND_CONTROLS:
            s_ = lig[(lig["control"] == ctrl) & (lig["relation"] == "reversing")]
            best = s_["adj_P"].min() if len(s_) else np.nan
            okc = bool(len(s_) and best < 0.05)
            if len(s_):
                pcs.append({"control": "(c) cytokine reverses (Enrichr)", "compound": ctrl, "cell_type": "",
                            "label": f"{ctrl} -log10 adj.P (Enrichr)", "value": float(-np.log10(best)),
                            "detail": s_.nsmallest(1, "adj_P")["term"].iloc[0], "passed": okc})
            gates.append({"gate": f"POSITIVE CONTROL (c): {ctrl} signatures reverse the death signature (Enrichr)",
                          "value": f"{len(s_)} terms; best adj. P = {best:.2g}" if len(s_) else "no term found",
                          "criterion": "adj. P < 0.05", "passed": okc, "severity": "warning"})
    except Exception as exc:
        LOG.warning(f"Enrichr ligand libraries not reachable ({exc}); control (c) skipped")
    # ---------------- DGIdb target layer ----------------
    lr = pd.read_csv(C.TAB / "S21_ligand_receptor_bulk.tsv", sep="\t")
    axes = lr.loc[lr["death_FDR"] < 0.05, "axis"].tolist()
    genes = sorted({g for a_ in axes for side in A2.LR_PAIRS.get(a_, ([], [])) for g in side} & set(A3.LR_TARGETS))
    genes += sorted(A3.CONTROL_RECEPTORS)
    try:
        dg = dgidb(genes)
        dg.to_csv(C.TAB / "S26_dgidb_targets.tsv", sep="\t", index=False)
        ctl = dg[(dg["role"] == "positive control") & dg["drug"].str.contains(A3.CONTROL_DRUGS, case=False, regex=True)]
        gates.append({"gate": "POSITIVE CONTROL (d): DGIdb lists GM-CSF / IFN-gamma drugs for their receptors",
                      "value": ", ".join(sorted(ctl["drug"].unique())) or "none", "criterion": ">= 1",
                      "passed": bool(len(ctl)), "severity": "warning"})
        LOG.info("DGIdb approved drugs per gene (in the desired direction):\n" + dg[dg["approved"]].groupby("gene").agg(
            drugs=("drug", "nunique"), desired=("matches_desired_direction", "sum")).to_string())
    except Exception as exc:
        LOG.warning(f"DGIdb not reachable ({exc}); target layer skipped")
        dg = pd.DataFrame(columns=["gene", "role", "drug", "approved", "matches_desired_direction"])
    pc = pd.DataFrame(pcs, columns=["control", "compound", "cell_type", "label", "value", "detail", "passed"])
    pc.to_csv(C.TAB / "S25_positive_controls.tsv", sep="\t", index=False)
    gates.append({"gate": f"RESULT (reported): candidates in {primary} (Tier 1 / Tier 2 / caution)",
                  "value": f"{int((cand['tier'] == 'Tier 1').sum())} / {int((cand['tier'] == 'Tier 2').sum())} / "
                           f"{int((cand['tier'] == 'Caution').sum())}: " + ", ".join(cand["compound"].head(8)),
                  "criterion": "reported", "passed": True, "severity": "warning"})
    record_gates("analysis3", gates)
    pc_plot = pc[pc["control"] != "(a) GR targets induced"]
    pc_a = pc[pc["control"] == "(a) GR targets induced"].sort_values("value", ascending=False).drop_duplicates("compound")
    safe_plot(plot_fig5, gtab, res, cand, pd.concat([pc_a, pc_plot]).head(14), dg, primary)
    LOG.info(f"Run time {(time.time() - t0) / 60:.1f} min")
    LOG.info(gate_summary().replace("STAGE 1 QUALITY GATES", "ANALYSIS STAGE 3 QUALITY GATES").replace(
        "Stage 2 can start after the field audit has been checked.", "Analysis stage 3 complete.").replace(
        "do not start Stage 2", "check the failed gates"))


# #####################################################################################
# PART E - ANALYSIS STAGE 4 (closing analyses): severity adjustment, incremental discrimination, comparison with
#          an established mortality signature, and the in-vivo glucocorticoid-response score
# #####################################################################################
LOG = get_logger("sepsis_analysis4")
GATES_FILE = C.TAB / "S0_analysis4_quality_gates.tsv"
_GATES.clear()


class A4:
    GC_STAT = 3.0               # OP3 statistic sign(logFC) x -log10(P): |value| >= 3 (P <= 1e-3) in EVERY glucocorticoid
    MAX_GC_GENES = 50           # per direction, ranked by the mean statistic
    MIN_GC_GENES = 10           # gate per direction
    MIN_SEV_COMPLETE = 0.70     # severity score available in >= 70% of a cohort's patients
    N_BOOT = 500
    # Sweeney et al. (Nat Commun 2018; prospective 11-gene version, Crit Care Med / PMC8219678): higher = higher risk
    SWEENEY_UP = {"DEFA4": [], "CD163": [], "RGS1": [], "PER1": [], "HIF1A": [], "SELENOP": ["SEPP1"],
                  "IFTAP": ["C11orf74"], "CIT": []}
    SWEENEY_DN = {"LY86": [], "TST": [], "KCNJ2": []}
    MIN_SWEENEY = (5, 2)        # genes needed (UP, DOWN) per cohort


AMENDMENT_7 = """# Amendment 7 - specification of the closing analyses (stage 05); written before any of them was run

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
"""


def write_amendment_7() -> None:
    p = C.PREREG / "AMENDMENT_7.md"
    if p.exists():
        old = _verify("AMENDMENT_7", AMENDMENT_7)
        LOG.info(f"Amendment 7 present (registered {old['registered_utc']}, unchanged)")
        return
    if (C.TAB / "S28_gc_score_results.tsv").exists():
        raise SystemExit("Stage-05 results exist but Amendment 7 is missing - the amendment must precede the analysis.")
    _register("AMENDMENT_7", AMENDMENT_7)


# -------------------------------------------------------------------------------------
# helpers
# -------------------------------------------------------------------------------------
def cohort_z(cohort: str, idx) -> pd.DataFrame:
    E = pd.read_csv(C.PROC / f"expr_{cohort}.tsv.gz", sep="\t", index_col=0)
    E = E[[i for i in idx if i in E.columns]]
    sd = E.std(axis=1, ddof=0)
    E = E[sd > 0]
    return E.sub(E.mean(axis=1), axis=0).div(E.std(axis=1, ddof=0), axis=0)


def set_score(Z: pd.DataFrame, up: list, dn: list) -> tuple[pd.Series, int, int]:
    u, d = [g for g in up if g in Z.index], [g for g in dn if g in Z.index]
    s = (Z.loc[u].mean(axis=0) if u else 0) - (Z.loc[d].mean(axis=0) if d else 0)
    s = pd.Series(s, index=Z.columns, dtype=float)
    return (s - s.mean()) / s.std(ddof=0), len(u), len(d)


def resolve(Zindex, spec: dict) -> list:
    out = []
    for g, aliases in spec.items():
        hit = next((x for x in [g] + aliases if x in Zindex), None)
        if hit:
            out.append(hit)
    return out


def auc(y: np.ndarray, s: np.ndarray) -> float:
    return float(stats.mannwhitneyu(s[y == 1], s[y == 0]).statistic / ((y == 1).sum() * (y == 0).sum()))


def delong(y: np.ndarray, s1: np.ndarray, s2: np.ndarray) -> dict:
    """Paired DeLong test (Sun & Xu 2014) for two correlated AUCs."""
    pos, neg = y == 1, y == 0
    m, n = pos.sum(), neg.sum()

    def comps(s):
        x, z = s[pos], s[neg]
        v10 = np.array([(np.sum(xi > z) + 0.5 * np.sum(xi == z)) / n for xi in x])
        v01 = np.array([(np.sum(x > zj) + 0.5 * np.sum(x == zj)) / m for zj in z])
        return v10.mean(), v10, v01
    a1, v10a, v01a = comps(s1)
    a2, v10b, v01b = comps(s2)
    s10 = np.cov(np.vstack([v10a, v10b]))
    s01 = np.cov(np.vstack([v01a, v01b]))
    S = s10 / m + s01 / n
    var = S[0, 0] + S[1, 1] - 2 * S[0, 1]
    se = float(np.sqrt(max(var, 1e-12)))
    return {"auc1": a1, "auc2": a2, "diff": a2 - a1, "se": se, "P": float(2 * stats.norm.sf(abs((a2 - a1) / se))),
            "se1": float(np.sqrt(S[0, 0])), "se2": float(np.sqrt(S[1, 1]))}


def design(d: pd.DataFrame, cols: list) -> np.ndarray:
    return sm.add_constant(d[cols].astype(float), has_constant="add").to_numpy()


def ols_age(d: pd.DataFrame, s: str) -> dict:
    dd = d.dropna(subset=[s, "age10", "male"])
    xc = ["age10", "male"] if dd["male"].nunique() > 1 else ["age10"]
    r = sm.OLS(dd[s].astype(float), sm.add_constant(dd[xc].astype(float))).fit()
    return {"n": len(dd), "b": float(r.params["age10"]), "se": float(r.bse["age10"])}


def spearman_fz(x: pd.Series, y: pd.Series) -> dict:
    ok = x.notna() & y.notna()
    rho = stats.spearmanr(x[ok], y[ok]).statistic
    n = int(ok.sum())
    return {"rho": float(rho), "z": float(np.arctanh(np.clip(rho, -0.999, 0.999))), "se": 1 / np.sqrt(max(n - 3, 1)),
            "n": n}


def optimism_auc(d: pd.DataFrame, base: list, full: list, n_boot: int, seed: int) -> dict:
    """Within-cohort AUC of base vs full logistic model, Harrell bootstrap optimism correction."""
    dd = d.dropna(subset=list(dict.fromkeys(base + full + ["death"])))
    y = dd["death"].astype(int).to_numpy()

    def fit_auc(train, test):
        out = []
        for cols in (base, full):
            r = logit_fit(train["death"].astype(int), train[cols].astype(float))
            out.append(auc(test["death"].astype(int).to_numpy(), design(test, cols) @ r.params.to_numpy()))
        return out
    app = fit_auc(dd, dd)
    rng = np.random.default_rng(seed)
    opt = []
    for _ in range(n_boot):
        b = dd.iloc[rng.integers(0, len(dd), len(dd))]
        if b["death"].nunique() < 2:
            continue
        try:
            ab, ao = fit_auc(b, b), fit_auc(b, dd)
            opt.append([ab[0] - ao[0], ab[1] - ao[1]])
        except Exception:
            continue
    o = np.mean(opt, axis=0) if opt else np.array([np.nan, np.nan])
    return {"n": len(dd), "events": int(y.sum()), "auc_base": app[0] - o[0], "auc_full": app[1] - o[1],
            "delta": (app[1] - o[1]) - (app[0] - o[0]), "boots_ok": len(opt)}


def gc_signature(common: set, s20: set) -> pd.DataFrame:
    T, M, _, _ = op3_profiles(op3_files())
    my = M["cell_type"].str.contains(A3.PRIMARY_CT, case=False, regex=True)
    gc = M["sm_name"].str.contains(A3.GLUCOCORTICOIDS, case=False, regex=True)
    rows = M.index[my & gc]
    names = M.loc[rows, "sm_name"].tolist()
    if len(rows) < 2:
        raise SystemExit(f"Fewer than 2 glucocorticoid profiles in OP3 myeloid cells ({names}) - GC score not defined.")
    G = T.loc[rows].T
    G.columns = names
    excl = s20 | {g for v in C.GENE_SETS.values() for g in v}
    G = G[G.index.isin(common) & ~G.index.isin(excl)]
    mn = G.mean(axis=1)
    up = G[(G >= A4.GC_STAT).all(axis=1)].index
    dn = G[(G <= -A4.GC_STAT).all(axis=1)].index
    up = mn.loc[up].sort_values(ascending=False).head(A4.MAX_GC_GENES)
    dn = mn.loc[dn].sort_values().head(A4.MAX_GC_GENES)
    out = pd.concat([pd.DataFrame({"gene": up.index, "direction": "UP", "mean_statistic": up.values}),
                     pd.DataFrame({"gene": dn.index, "direction": "DOWN", "mean_statistic": dn.values})])
    for n_ in names:
        out[n_] = G.loc[out["gene"], n_].values
    return out


def pooled_row(label: str, est: list, se: list, exp_: bool = True) -> dict:
    m = dl_meta(np.array(est, float), np.array(se, float))
    f = np.exp if exp_ else (lambda v: v)
    return {"label": label, "est": m["est"], "se": m["se"], "value": f(m["est"]), "lo": f(m["est"] - 1.96 * m["se"]),
            "hi": f(m["est"] + 1.96 * m["se"]), "P": m["P"], "k": m["k"], "I2": m["I2"]}


# -------------------------------------------------------------------------------------
# figure 6
# -------------------------------------------------------------------------------------
def lin_forest(ax, rows: pd.DataFrame, title: str, xlabel: str, ref: float, log: bool = False) -> None:
    y = np.arange(len(rows))[::-1]
    for yy, (_, r) in zip(y, rows.iterrows()):
        pooled = str(r["label"]).startswith(("Pooled", "All", "Validation"))
        col_ = OKABE[1] if pooled else (OKABE[0] if r.get("role") == "discovery" else "#555555")
        ax.plot([r["lo"], r["hi"]], [yy, yy], color=col_, lw=0.8)
        ax.scatter([r["value"]], [yy], marker="D" if pooled else "s", s=18 if pooled else 10, color=col_, zorder=3)
        ax.text(1.03, yy, f"{r['value']:.2f} ({r['lo']:.2f} to {r['hi']:.2f})", transform=ax.get_yaxis_transform(),
                va="center", fontsize=4.8)
    ax.axvline(ref, color="#777777", lw=0.5, ls="--")
    if log:
        log_axis(ax)
    else:
        ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(4))
    ax.set_yticks(y)
    ax.set_yticklabels(rows["label"], fontsize=4.8)
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=6)


def plot_fig6(sev: pd.DataFrame, dauc: pd.DataFrame, gcd: pd.DataFrame, gca: pd.DataFrame, cor: pd.DataFrame) -> None:
    set_style()
    fig = plt.figure(figsize=(FIG_W, 200 * MM))
    gs = fig.add_gridspec(3, 2, hspace=0.75, wspace=1.45, left=0.17, right=0.86, top=0.95, bottom=0.06)
    panels = [(gs[0, 0], sev, "ML index: unadjusted vs\nseverity-adjusted", "Odds ratio per SD", 1, True, "a"),
              (gs[0, 1], dauc, "Added value of the ML index\n(external, models from discovery)", "Delta AUC", 0, False,
               "b"),
              (gs[1, 0], gcd, "Glucocorticoid-response score\nand 28-day death", "Odds ratio per SD", 1, True, "c"),
              (gs[1, 1], gca, "Glucocorticoid-response score\nand age", "Change in score (SD) per 10 years", 0,
               False, "d"),
              (gs[2, 0], cor, "Glucocorticoid-response score:\ncorrelation with MHC-II and ML index", "Spearman rho", 0, False, "e")]
    for spec, rows, title, xl, ref, lg, lab in panels:
        ax = fig.add_subplot(spec)
        if rows is not None and len(rows):
            lin_forest(ax, rows.reset_index(drop=True), title, xl, ref, lg)
        else:
            ax.text(0.5, 0.5, "not available", ha="center", va="center", fontsize=6)
            ax.set_title(title, fontsize=6)
        panel_label(ax, lab, x=-0.75)
    save_figure(fig, "Fig6_closing_analyses")


# -------------------------------------------------------------------------------------
# run
# -------------------------------------------------------------------------------------
def run_analysis4() -> None:
    t0 = time.time()
    warnings.filterwarnings("ignore", category=FutureWarning)
    LOG.info(f"Project folder: {C.PROJECT_DIR}")
    preflight()
    preflight3()
    if not (C.TAB / "S23_op3_compound_results.tsv").exists():
        raise SystemExit("Stage 4 (04_drug_repurposing.py) has not finished - S23_op3_compound_results.tsv missing.")
    write_preregistration()
    write_amendment()
    write_amendment_3()
    write_amendment_4()
    write_amendment_5()
    write_amendment_7()
    clin = pd.read_csv(C.TAB / "S3_harmonised_clinical.tsv", sep="\t", index_col=0, low_memory=False)
    clin["eligible"] = to_bool(clin["eligible"])
    disc, val = analysed_cohorts()
    cohorts = [disc] + val
    data = {}
    for c in cohorts:
        try:
            data[c] = cohort_scores(c, clin)[0]
        except Exception as exc:
            LOG.warning(f"{c}: scores not computed ({exc}) - excluded from stage 05")
    cohorts = [c for c in cohorts if c in data]
    gates = []
    common = None
    for c in cohorts:
        g = set(pd.read_csv(C.PROC / f"expr_{c}.tsv.gz", sep="\t", index_col=0, usecols=[0]).index)
        common = g if common is None else common & g
    s20 = set(pd.read_csv(C.TAB / "S20_reversal_signature_prespecified.tsv", sep="\t")["gene"])
    # ---------------- GC signature (from OP3) ----------------
    try:
        gsig = gc_signature(common, s20)
    except (SystemExit, Exception) as exc:
        LOG.warning(f"GC-response signature not available ({exc}) - analysis (D) skipped, (A)-(C) still run")
        gsig = pd.DataFrame(columns=["gene", "direction", "mean_statistic"])
    gsig.to_csv(C.TAB / "S27_gc_response_signature.tsv", sep="\t", index=False)
    gup = gsig.loc[gsig["direction"] == "UP", "gene"].tolist()
    gdn = gsig.loc[gsig["direction"] == "DOWN", "gene"].tolist()
    LOG.info(f"Glucocorticoid-response signature (OP3 myeloid, all GCs concordant, |stat| >= {A4.GC_STAT}): "
             f"UP {len(gup)}, DOWN {len(gdn)}\n  UP: {', '.join(gup[:25])}\n  DOWN: {', '.join(gdn[:25])}")
    gc_ok = min(len(gup), len(gdn)) >= A4.MIN_GC_GENES
    gates.append({"gate": "GC-response signature genes (UP / DOWN)", "value": f"{len(gup)} / {len(gdn)}",
                  "criterion": f">= {A4.MIN_GC_GENES} per direction", "passed": gc_ok, "severity": "warning"})
    # ---------------- per-cohort scores ----------------
    sw_cov = {}
    for c in cohorts:
        d = data[c]
        Z = cohort_z(c, d.index)
        if gc_ok:
            d["gc_score"], _, _ = set_score(Z, gup, gdn)
        su, sd_ = resolve(Z.index, A4.SWEENEY_UP), resolve(Z.index, A4.SWEENEY_DN)
        sw_cov[c] = (len(su), len(sd_))
        d["sweeney"] = set_score(Z, su, sd_)[0] if (len(su) >= A4.MIN_SWEENEY[0] and len(sd_) >= A4.MIN_SWEENEY[1]) \
            else np.nan
        data[c] = d
    LOG.info("Sweeney genes measured (UP / DOWN): " + "; ".join(f"{c} {u}/{d_}" for c, (u, d_) in sw_cov.items()))
    pd.concat(data.values())[["cohort", "ml_index", "mhc_class_ii", "sweeney"] +
                             (["gc_score"] if gc_ok else [])].to_csv(C.TAB / "S28b_patient_scores_stage05.tsv", sep="\t")
    role = {c: ("discovery" if c == disc else "validation") for c in cohorts}
    # ---------------- (A) severity ----------------
    sev_rows = []
    for c in cohorts:
        d = data[c]
        for sv in [x for x in d.columns if x.startswith("sev_") and d[x].notna().mean() >= A4.MIN_SEV_COMPLETE]:
            dd = d.dropna(subset=[sv, "ml_index", "age10", "death"]).copy()
            dd["sev_z"] = (dd[sv] - dd[sv].mean()) / dd[sv].std()
            for pred in ["ml_index", "mhc_class_ii"] + (["gc_score"] if gc_ok else []):
                r0 = assoc(dd, pred)
                r1 = assoc(dd, pred, adjust=("age10", "male", "sev_z"))
                sev_rows.append({"cohort": c, "severity": sv[4:], "predictor": pred, "n": r1["n"],
                                 "events": r1["events"], "logOR_unadj": r0["logratio"], "se_unadj": r0["se"],
                                 "logOR_adj": r1["logratio"], "se_adj": r1["se"]})
    sev = pd.DataFrame(sev_rows)
    sev_plot, sev_pool = [], {}
    if len(sev):
        for pred, g in sev.groupby("predictor"):
            for kind in ("unadj", "adj"):
                pr = pooled_row(f"Pooled {pred} ({kind})", g[f"logOR_{kind}"], g[f"se_{kind}"])
                sev_pool[(pred, kind)] = pr
        for _, r in sev[sev["predictor"] == "ml_index"].iterrows():
            for kind in ("unadj", "adj"):
                e, s = r[f"logOR_{kind}"], r[f"se_{kind}"]
                sev_plot.append({"label": f"{r['cohort']} ({'+ ' + r['severity'][:12] if kind == 'adj' else 'unadjusted'})",
                                 "value": np.exp(e), "lo": np.exp(e - 1.96 * s), "hi": np.exp(e + 1.96 * s),
                                 "role": "validation"})
        for kind in ("unadj", "adj"):
            pr = sev_pool[("ml_index", kind)]
            sev_plot.append({**pr, "label": f"Pooled ({'severity-adjusted' if kind == 'adj' else 'unadjusted'})"})
    sev.to_csv(C.TAB / "S29_severity_adjusted.tsv", sep="\t", index=False)
    if len(sev):
        pa = sev_pool[("ml_index", "adj")]
        LOG.info(f"(A) ML index adjusted for severity ({sev['cohort'].nunique()} cohorts): pooled OR "
                 f"{pa['value']:.2f} ({pa['lo']:.2f}-{pa['hi']:.2f}), P = {pa['P']:.2g}")
        gates.append({"gate": "(A) ML index adjusted for severity (pooled OR, cohorts with severity scores)",
                      "value": f"{pa['value']:.2f} ({pa['lo']:.2f}-{pa['hi']:.2f}), P = {pa['P']:.2g}, "
                               f"k = {pa['k']}", "criterion": "reported", "passed": True, "severity": "warning"})
    # ---------------- (B) incremental discrimination ----------------
    dd0 = data[disc].dropna(subset=["age10", "male", "ml_index", "death"])
    mods = {"base": ["age10", "male"], "base+ML": ["age10", "male", "ml_index"]}
    if data[disc]["sweeney"].notna().mean() > 0.9:
        mods.update({"base+Sweeney": ["age10", "male", "sweeney"], "base+Sweeney+ML": ["age10", "male", "sweeney",
                                                                                       "ml_index"]})
    fits = {k: logit_fit(dd0.dropna(subset=v)["death"].astype(int), dd0.dropna(subset=v)[v].astype(float))
            for k, v in mods.items()}
    auc_rows = []
    for c in val:
        for b_, f_ in (("base", "base+ML"), ("base+Sweeney", "base+Sweeney+ML")):
            if f_ not in fits:
                continue
            dd = data[c].dropna(subset=mods[f_] + ["death"])
            if len(dd) < 20 or dd["death"].nunique() < 2:
                continue
            y = dd["death"].astype(int).to_numpy()
            r = delong(y, design(dd, mods[b_]) @ fits[b_].params.to_numpy(),
                       design(dd, mods[f_]) @ fits[f_].params.to_numpy())
            auc_rows.append({"cohort": c, "comparison": f"{f_} vs {b_}", "n": len(dd), "events": int(y.sum()), **r})
    au = pd.DataFrame(auc_rows)
    dauc_plot = []
    if len(au):
        for comp, g in au.groupby("comparison", sort=False):
            pr = pooled_row(f"Pooled: {comp}", g["diff"], g["se"], exp_=False)
            au = pd.concat([au, pd.DataFrame([{"cohort": "POOLED", "comparison": comp, "diff": pr["est"],
                                               "se": pr["se"], "P": pr["P"], "auc1": g["auc1"].mean(),
                                               "auc2": g["auc2"].mean(), "n": g["n"].sum(),
                                               "events": g["events"].sum()}])], ignore_index=True)
            LOG.info(f"(B) {comp}: pooled delta AUC {pr['est']:+.3f} ({pr['lo']:+.3f} to {pr['hi']:+.3f}), "
                     f"P = {pr['P']:.2g}; mean AUC {g['auc1'].mean():.3f} -> {g['auc2'].mean():.3f}")
            gates.append({"gate": f"(B) external delta AUC, {comp} (pooled over validation cohorts)",
                          "value": f"{pr['est']:+.3f} ({pr['lo']:+.3f} to {pr['hi']:+.3f}), P = {pr['P']:.2g}",
                          "criterion": "reported", "passed": True, "severity": "warning"})
            if comp.startswith("base+ML"):
                for _, r in g.iterrows():
                    dauc_plot.append({"label": r["cohort"], "value": r["diff"], "lo": r["diff"] - 1.96 * r["se"],
                                      "hi": r["diff"] + 1.96 * r["se"], "role": "validation"})
            dauc_plot.append({**pr, "label": f"Pooled ({comp.split(' vs ')[1]})"})
    within = []
    for c in cohorts:
        d = data[c]
        for sv in [x for x in d.columns if x.startswith("sev_") and d[x].notna().mean() >= A4.MIN_SEV_COMPLETE]:
            dd = d.copy()
            dd["sev_z"] = (dd[sv] - dd[sv].mean()) / dd[sv].std()
            base = ["age10", "sev_z"] + (["male"] if dd["male"].nunique() > 1 else [])
            try:
                within.append({"cohort": c, "severity": sv[4:],
                               **optimism_auc(dd, base, base + ["ml_index"], A4.N_BOOT, C.SEED)})
            except Exception as exc:
                LOG.warning(f"{c}: within-cohort AUC failed ({exc})")
    au.to_csv(C.TAB / "S30_incremental_auc.tsv", sep="\t", index=False)
    pd.DataFrame(within).to_csv(C.TAB / "S30b_within_cohort_auc_severity.tsv", sep="\t", index=False)
    if within:
        LOG.info("(B) Within severity cohorts (optimism-corrected):\n" + pd.DataFrame(within).round(3).to_string(index=False))
    # ---------------- (C) Sweeney comparison ----------------
    swr = []
    for c in cohorts:
        d = data[c]
        if d["sweeney"].isna().all():
            continue
        r_sw = assoc(d, "sweeney")
        r_ml = assoc(d, "ml_index", adjust=("age10", "male", "sweeney"))
        cr = spearman_fz(d["sweeney"], d["ml_index"])
        swr.append({"cohort": c, "role": role[c], "sweeney_logOR": r_sw["logratio"], "sweeney_se": r_sw["se"],
                    "ml_adj_sweeney_logOR": r_ml["logratio"], "ml_adj_sweeney_se": r_ml["se"],
                    "rho_sweeney_ml": cr["rho"], "rho_z": cr["z"], "rho_se": cr["se"], "n": r_ml["n"]})
    sw = pd.DataFrame(swr)
    if len(sw):
        p1 = pooled_row("Sweeney score", sw["sweeney_logOR"], sw["sweeney_se"])
        p2 = pooled_row("ML index adjusted for Sweeney", sw["ml_adj_sweeney_logOR"], sw["ml_adj_sweeney_se"])
        p3 = dl_meta(sw["rho_z"].to_numpy(), sw["rho_se"].to_numpy())
        LOG.info(f"(C) Sweeney score OR {p1['value']:.2f} (P = {p1['P']:.2g}); ML index adjusted for Sweeney OR "
                 f"{p2['value']:.2f} ({p2['lo']:.2f}-{p2['hi']:.2f}), P = {p2['P']:.2g}; rho(Sweeney, ML) "
                 f"{np.tanh(p3['est']):.2f}")
        gates.append({"gate": "(C) ML index adjusted for the Sweeney mortality score (all cohorts pooled)",
                      "value": f"OR {p2['value']:.2f} ({p2['lo']:.2f}-{p2['hi']:.2f}), P = {p2['P']:.2g}; "
                               f"Sweeney alone OR {p1['value']:.2f}; rho {np.tanh(p3['est']):.2f}",
                      "criterion": "reported", "passed": True, "severity": "warning"})
    sw.to_csv(C.TAB / "S31_sweeney_comparison.tsv", sep="\t", index=False)
    # ---------------- (D) GC-response score ----------------
    gcd_plot, gca_plot, cor_plot, gres = [], [], [], []
    if gc_ok:
        for c in cohorts:
            d = data[c]
            r1 = assoc(d, "gc_score")
            r2 = assoc(d, "gc_score", adjust=("age10", "male", "ml_index", "mhc_class_ii"))
            ra = ols_age(d, "gc_score")
            cm = spearman_fz(d["gc_score"], d["mhc_class_ii"])
            cl = spearman_fz(d["gc_score"], d["ml_index"])
            gres.append({"cohort": c, "role": role[c], "n": r1["n"], "events": r1["events"],
                         "death_logOR": r1["logratio"], "death_se": r1["se"],
                         "death_adj_logOR": r2["logratio"], "death_adj_se": r2["se"], "age_b": ra["b"],
                         "age_se": ra["se"], "rho_mhc": cm["rho"], "rho_mhc_z": cm["z"], "rho_mhc_se": cm["se"],
                         "rho_ml": cl["rho"], "rho_ml_z": cl["z"], "rho_ml_se": cl["se"]})
        g = pd.DataFrame(gres)
        pdth = pooled_row("All cohorts pooled", g["death_logOR"], g["death_se"])
        pdadj = pooled_row("Pooled, adj. ML + MHC-II", g["death_adj_logOR"], g["death_adj_se"])
        page = pooled_row("All cohorts pooled", g["age_b"], g["age_se"], exp_=False)
        holm = multipletests([pdth["P"], page["P"]], method="holm")[1]
        pm_ = dl_meta(g["rho_mhc_z"].to_numpy(), g["rho_mhc_se"].to_numpy())
        pl_ = dl_meta(g["rho_ml_z"].to_numpy(), g["rho_ml_se"].to_numpy())
        summ = pd.DataFrame([
            {"test": "D1 death ~ GC score + age + sex (primary)", "estimate": pdth["value"], "lo": pdth["lo"],
             "hi": pdth["hi"], "P": pdth["P"], "P_holm": holm[0], "k": pdth["k"], "I2": pdth["I2"]},
            {"test": "D2 GC score ~ age per 10 y + sex (primary)", "estimate": page["value"], "lo": page["lo"],
             "hi": page["hi"], "P": page["P"], "P_holm": holm[1], "k": page["k"], "I2": page["I2"]},
            {"test": "death ~ GC score, adjusted for ML index + MHC-II (secondary)", "estimate": pdadj["value"],
             "lo": pdadj["lo"], "hi": pdadj["hi"], "P": pdadj["P"], "P_holm": np.nan, "k": pdadj["k"],
             "I2": pdadj["I2"]},
            {"test": "Spearman rho GC score vs MHC-II (secondary)", "estimate": np.tanh(pm_["est"]),
             "lo": np.tanh(pm_["est"] - 1.96 * pm_["se"]), "hi": np.tanh(pm_["est"] + 1.96 * pm_["se"]),
             "P": pm_["P"], "P_holm": np.nan, "k": pm_["k"], "I2": pm_["I2"]},
            {"test": "Spearman rho GC score vs ML index (secondary)", "estimate": np.tanh(pl_["est"]),
             "lo": np.tanh(pl_["est"] - 1.96 * pl_["se"]), "hi": np.tanh(pl_["est"] + 1.96 * pl_["se"]),
             "P": pl_["P"], "P_holm": np.nan, "k": pl_["k"], "I2": pl_["I2"]}])
        g.to_csv(C.TAB / "S28a_gc_score_per_cohort.tsv", sep="\t", index=False)
        summ.to_csv(C.TAB / "S28_gc_score_results.tsv", sep="\t", index=False)
        LOG.info("(D) Glucocorticoid-response score:\n" + summ.round(4).to_string(index=False))
        for t_, r_ in zip(("D1 death", "D2 age"), summ.iloc[:2].itertuples()):
            gates.append({"gate": f"(D) PRIMARY {t_}: GC-response score (all cohorts pooled)",
                          "value": f"{r_.estimate:.2f} ({r_.lo:.2f}-{r_.hi:.2f}), P = {r_.P:.2g}, Holm {r_.P_holm:.2g}, "
                                   f"I2 {r_.I2 * 100:.0f}%", "criterion": "reported (Holm < 0.05 = supported)",
                          "passed": True, "severity": "warning"})
        for _, r in g.iterrows():
            e, s = r["death_logOR"], r["death_se"]
            gcd_plot.append({"label": r["cohort"], "value": np.exp(e), "lo": np.exp(e - 1.96 * s),
                             "hi": np.exp(e + 1.96 * s), "role": r["role"]})
            gca_plot.append({"label": r["cohort"], "value": r["age_b"], "lo": r["age_b"] - 1.96 * r["age_se"],
                             "hi": r["age_b"] + 1.96 * r["age_se"], "role": r["role"]})
        gcd_plot += [{**pdth}, {**pdadj}]
        gca_plot += [{**page}]
        for key_, lab_, m_ in (("mhc", "Pooled: vs MHC class II", pm_), ("ml", "Pooled: vs ML index", pl_)):
            for _, r in g.iterrows():
                z_, s_ = r[f"rho_{key_}_z"], r[f"rho_{key_}_se"]
                cor_plot.append({"label": f"{r['cohort']} ({'MHC-II' if key_ == 'mhc' else 'ML'})",
                                 "value": r[f"rho_{key_}"], "lo": np.tanh(z_ - 1.96 * s_), "hi": np.tanh(z_ + 1.96 * s_),
                                 "role": r["role"]})
            cor_plot.append({"label": lab_, "value": np.tanh(m_["est"]), "lo": np.tanh(m_["est"] - 1.96 * m_["se"]),
                             "hi": np.tanh(m_["est"] + 1.96 * m_["se"])})
        if len(sev):
            ga = sev_pool.get(("gc_score", "adj"))
            if ga:
                LOG.info(f"(A) GC score adjusted for severity: pooled OR {ga['value']:.2f} ({ga['lo']:.2f}-"
                         f"{ga['hi']:.2f}), P = {ga['P']:.2g}")
    record_gates("analysis4", gates)
    safe_plot(plot_fig6, pd.DataFrame(sev_plot), pd.DataFrame(dauc_plot), pd.DataFrame(gcd_plot),
              pd.DataFrame(gca_plot), pd.DataFrame(cor_plot))
    # ---------------- final summary of the whole project ----------------
    rows = []

    def add(stage, finding, value):
        rows.append({"stage": stage, "finding": finding, "result": value})
    try:
        h = pd.read_csv(C.TAB / "S11_H1_ml_index_mortality.tsv", sep="\t")
        for _, r in h[h["role"].isin(["discovery", "pooled"])].iterrows():
            add("2 (H1)", f"ML index and death: {r['label']}", f"{r['ratio']:.2f} ({r['lo']:.2f}-{r['hi']:.2f}), "
                                                               f"P = {r['P']:.2g}")
    except Exception:
        pass
    try:
        p17 = pd.read_csv(C.TAB / "S17_primary_results.tsv", sep="\t")
        for _, r in p17.iterrows():
            add("2 (primary)", " | ".join(str(r[c_]) for c_ in p17.columns[:2]),
                " | ".join(str(r[c_]) for c_ in p17.columns[2:6]))
    except Exception:
        pass
    try:
        lr = pd.read_csv(C.TAB / "S21_ligand_receptor_bulk.tsv", sep="\t")
        for _, r in lr[lr["death_FDR"] < 0.05].iterrows():
            add("3 (ligand-receptor)", r["axis"], f"OR {r['death_OR']:.2f}, FDR {r['death_FDR']:.2g}; adjusted for ML "
                                                  f"{r['death_OR_adjML']:.2f}, FDR {r['adjML_FDR']:.2g}")
    except Exception:
        pass
    try:
        cd = pd.read_csv(C.TAB / "S24_op3_candidates.tsv", sep="\t")
        add("4 (drugs)", "candidates Tier 1 / Tier 2 / caution", " / ".join(str(int((cd["tier"] == t).sum()))
                                                                          for t in ("Tier 1", "Tier 2", "Caution")))
        gc_ = cd[cd["glucocorticoid"].astype(str).str.lower() == "true"]
        for _, r in gc_.iterrows():
            add("4 (drugs)", f"glucocorticoid {r['compound']} (myeloid)", f"reversal z {r['R_z_primary']:.2f}; "
                                                                           f"MHC-II z {r['mhc_z_primary']:.2f}")
    except Exception:
        pass
    for gt in _GATES:
        add("5 (closing)", gt["gate"], gt["value"])
    pd.DataFrame(rows).to_csv(C.TAB / "S32_project_summary.tsv", sep="\t", index=False)
    LOG.info(f"Project summary written: {C.TAB / 'S32_project_summary.tsv'}")
    LOG.info(f"Run time {(time.time() - t0) / 60:.1f} min")
    LOG.info(gate_summary().replace("STAGE 1 QUALITY GATES", "CLOSING ANALYSES (STAGE 05) - RESULTS AND GATES").replace(
        "Stage 2 can start after the field audit has been checked.", "Stage 05 complete - project analyses finished.")
             .replace("do not start Stage 2", "check the failed gates"))


# #####################################################################################
# PART F - STAGE 06: external single-cell validation in septic patients (GSE167363, PBMC scRNA-seq)
# #####################################################################################
import shutil
from scipy import sparse

LOG = get_logger("sepsis_singlecell")
GATES_FILE = C.TAB / "S0_singlecell_quality_gates.tsv"
_GATES.clear()


class A5:
    ACC = "GSE167363"
    TAR_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167363/suppl/GSE167363_RAW.tar"
    # sample metadata as deposited in GEO (series matrix GSE167363; Qiu et al., J Leukoc Biol 2021)
    SAMPLES = {
        "GSM5102900": ("HC1", "donor1", "healthy", "T0", "male", "35-40"),
        "GSM5102901": ("HC2", "donor2", "healthy", "T0", "female", "45-50"),
        "GSM5102902": ("NS_LS_T0", "donor3", "non-survivor", "T0", "female", "65-70"),
        "GSM5102903": ("NS_LS_T6", "donor3", "non-survivor", "T6", "female", "65-70"),
        "GSM5102904": ("S1_T0", "donor4", "survivor", "T0", "female", "65-70"),
        "GSM5102905": ("S1_T6", "donor4", "survivor", "T6", "female", "65-70"),
        "GSM5511351": ("NS_ES_T0", "donor5", "non-survivor", "T0", "male", "90-95"),
        "GSM5511352": ("NS_ES_T6", "donor5", "non-survivor", "T6", "male", "90-95"),
        "GSM5511353": ("S2_T0", "donor6", "survivor", "T0", "male", "45-50"),
        "GSM5511354": ("S2_T6", "donor6", "survivor", "T6", "male", "45-50"),
        "GSM5511355": ("S3_T0", "donor7", "survivor", "T0", "female", "70-75"),
        "GSM5511356": ("S3_T6", "donor7", "survivor", "T6", "female", "70-75"),
    }
    MIN_GENES, MAX_MITO, MIN_CELLS_GENE = 200, 20.0, 3
    RAW_BARCODES = 50000            # more barcodes than this = unfiltered matrix: keep barcodes with >= 500 UMIs
    N_HVG, N_PCS = 2000, 30
    MIN_MONO = 100                  # CD14 monocytes per sample for the within-patient correlation
    MARKERS = {
        "CD14 monocytes": ["CD14", "LYZ", "VCAN", "FCN1", "S100A12"],
        "CD16 monocytes": ["FCGR3A", "MS4A7", "CDKN1C", "RHOC"],
        "Dendritic cells": ["FCER1A", "CD1C", "CLEC10A"],
        "pDC": ["LILRA4", "CLEC4C", "IL3RA"],
        "B cells": ["MS4A1", "CD79A", "CD79B", "BANK1"],
        "Plasma cells": ["MZB1", "JCHAIN", "TNFRSF17"],
        "CD4 T cells": ["CD3E", "IL7R", "CD4", "CCR7"],
        "CD8 T cells": ["CD3E", "CD8A", "CD8B", "GZMK"],
        "NK cells": ["NKG7", "GNLY", "KLRD1", "KLRF1"],
        "Platelets": ["PPBP", "PF4", "GP9"],
        "Erythroid cells": ["HBB", "HBA1", "ALAS2", "CA1"],
        "Granulocytes": ["FCGR3B", "CSF3R", "CXCR2", "MMP9"],
    }
    MONO = "CD14 monocytes"
    CT_COLORS = ["#0072B2", "#56B4E9", "#009E73", "#8FD3B6", "#E69F00", "#F0C987", "#D55E00", "#F4A582", "#CC79A7",
                 "#999999", "#B2182B", "#6A3D9A"]


AMENDMENT_8 = """# Amendment 8 - external single-cell validation (stage 06); written before the single-cell data were downloaded

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
"""


def write_amendment_8() -> None:
    p = C.PREREG / "AMENDMENT_8.md"
    if p.exists():
        old = _verify("AMENDMENT_8", AMENDMENT_8)
        LOG.info(f"Amendment 8 present (registered {old['registered_utc']}, unchanged)")
        return
    if (C.TAB / "S35_sc_monocyte_coupling.tsv").exists():
        raise SystemExit("Single-cell results exist but Amendment 8 is missing - the amendment must precede the analysis.")
    _register("AMENDMENT_8", AMENDMENT_8)


# -------------------------------------------------------------------------------------
# data
# -------------------------------------------------------------------------------------
def sc_files() -> Path:
    import tarfile
    d = C.RAW / A5.ACC
    out = d / "mtx"
    if out.exists() and len(list(out.glob("*matrix.mtx.gz"))) >= len(A5.SAMPLES):
        return out
    tar = d / "GSE167363_RAW.tar"
    if not tar.exists():
        download_resumable(A5.TAR_URL, tar)
    out.mkdir(parents=True, exist_ok=True)
    try:
        with tarfile.open(tar) as t:
            for m in t.getmembers():
                name = Path(m.name).name
                if m.isfile() and name.startswith("GSM") and ".." not in m.name:
                    with t.extractfile(m) as src, open(out / name, "wb") as dst:
                        shutil.copyfileobj(src, dst)
    except tarfile.TarError as exc:
        tar.replace(tar.with_name(tar.name + ".broken"))
        raise SystemExit(f"{tar.name} is damaged ({exc}) - run the cell again to download it again. Manual: "
                         f"{A5.TAR_URL}")
    LOG.info(f"Extracted {len(list(out.iterdir()))} files to {out}")
    return out


def read_sample(folder: Path, gsm: str):
    import anndata as ad
    from scipy.io import mmread
    files = {f.name: f for f in folder.glob(f"{gsm}_*")}
    pick = lambda key: next((f for n, f in files.items() if key in n), None)
    mtx, bc = pick("matrix.mtx"), pick("barcodes")
    ft = pick("features") or pick("genes")
    if not (mtx and bc and ft):
        raise SystemExit(f"{gsm}: matrix / barcodes / features file missing in {folder} (found {sorted(files)})")
    X = sparse.csr_matrix(mmread(str(mtx)).T.astype(np.float32))
    barcodes = pd.read_csv(bc, header=None, sep="\t")[0].astype(str).to_numpy()
    feats = pd.read_csv(ft, header=None, sep="\t")
    sym = feats[1] if feats.shape[1] > 1 else feats[0]
    if feats.shape[1] > 2:
        keep = feats[2].astype(str).str.contains("Gene Expression", case=False).to_numpy()
        X, sym = X[:, keep], sym[keep]
    if X.shape != (len(barcodes), len(sym)):
        raise SystemExit(f"{gsm}: matrix {X.shape} does not match barcodes / features ({len(barcodes)}, {len(sym)})")
    if X.shape[0] > A5.RAW_BARCODES:
        umi = np.asarray(X.sum(axis=1)).ravel()
        X, barcodes = X[umi >= 500], barcodes[umi >= 500]
        LOG.info(f"  {gsm}: unfiltered matrix - kept {len(barcodes)} barcodes with >= 500 UMIs")
    a = ad.AnnData(X=X, obs=pd.DataFrame(index=[f"{gsm}_{b}" for b in barcodes]),
                   var=pd.DataFrame(index=pd.Index(sym.astype(str).to_numpy())))
    a.var_names_make_unique()
    name, donor, group, tp, sex, age = A5.SAMPLES[gsm]
    for k, v in (("gsm", gsm), ("sample", name), ("donor", donor), ("group", group), ("time", tp), ("sex", sex),
                 ("age_range", age)):
        a.obs[k] = v
    return a


def build_atlas(folder: Path):
    import anndata as ad
    import scanpy as sc
    cache = C.PROC / "sc_GSE167363_processed.h5ad"
    if cache.exists():
        LOG.info(f"Processed single-cell object loaded from cache ({cache.name})")
        return ad.read_h5ad(cache), {}
    parts = []
    for gsm in A5.SAMPLES:
        a = read_sample(folder, gsm)
        LOG.info(f"  {gsm} {A5.SAMPLES[gsm][0]}: {a.n_obs} cells x {a.n_vars} genes")
        parts.append(a)
    adata = ad.concat(parts, join="outer", fill_value=0)
    adata.obs_names_make_unique()
    n0 = adata.n_obs
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    adata = adata[(adata.obs["n_genes_by_counts"] >= A5.MIN_GENES) & (adata.obs["pct_counts_mt"] <= A5.MAX_MITO)].copy()
    sc.pp.filter_genes(adata, min_cells=A5.MIN_CELLS_GENE)
    LOG.info(f"QC: {adata.n_obs} of {n0} cells kept ({adata.n_vars} genes)")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=A5.N_HVG, flavor="seurat", batch_key="sample")
    ad_h = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(ad_h, max_value=10)
    sc.tl.pca(ad_h, n_comps=A5.N_PCS, random_state=C.SEED)
    adata.obsm["X_pca"] = ad_h.obsm["X_pca"]
    info = {"integration": "none"}
    rep = "X_pca"
    try:
        import harmonypy
        ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs[["sample"]], "sample", random_state=C.SEED)
        Zc = np.asarray(ho.Z_corr)
        Zc = Zc.T if Zc.shape[0] != adata.n_obs else Zc
        adata.obsm["X_pca_harmony"] = Zc
        rep, info["integration"] = "X_pca_harmony", "Harmony"
    except Exception as exc:
        LOG.warning(f"Harmony not used ({type(exc).__name__}: {exc}) - UMAP without integration (pre-specified fallback)")
    X = np.asarray(adata.obsm[rep])[:, :A5.N_PCS]
    graph = None
    try:
        sc.pp.neighbors(adata, use_rep=rep, n_neighbors=15, random_state=C.SEED)
        sc.tl.umap(adata, random_state=C.SEED)
        graph = "scanpy"
    except Exception as exc:                 # library-version mismatch (scanpy / scipy / scikit-learn)
        import umap as umap_
        from umap.umap_ import fuzzy_simplicial_set
        from sklearn.neighbors import NearestNeighbors
        LOG.warning(f"scanpy neighbours failed ({type(exc).__name__}) - kNN graph and UMAP built with "
                    f"scikit-learn / umap-learn directly (same method: 15 neighbours, UMAP connectivities)")
        dist, idx = NearestNeighbors(n_neighbors=15).fit(X).kneighbors(X)
        conn = fuzzy_simplicial_set(X, 15, random_state=np.random.RandomState(C.SEED), metric="euclidean",
                                    knn_indices=idx, knn_dists=dist)[0]
        adata.uns["_adjacency"] = sparse.csr_matrix(conn)
        adata.obsm["X_umap"] = umap_.UMAP(n_neighbors=15, random_state=C.SEED).fit_transform(X)
        graph = "manual"
    try:
        kw = dict(resolution=1.0, key_added="cluster", n_iterations=2, directed=False, random_state=C.SEED)
        if graph == "manual":
            kw["adjacency"] = adata.uns.pop("_adjacency")
        try:
            sc.tl.leiden(adata, flavor="igraph", **kw)
        except TypeError:                   # older scanpy without the 'flavor' argument
            kw.pop("n_iterations", None)
            kw.pop("directed", None)
            sc.tl.leiden(adata, **kw)
        info["clustering"] = "Leiden"
    except Exception as exc:
        adata.uns.pop("_adjacency", None)
        from sklearn.cluster import KMeans
        adata.obs["cluster"] = pd.Categorical(KMeans(25, n_init=10, random_state=C.SEED).fit_predict(X).astype(str))
        info["clustering"] = f"k-means (Leiden not available: {type(exc).__name__})"
        LOG.warning(f"Leiden not available - k-means clustering used ({exc})")
    info["graph"] = graph
    adata.uns["processing"] = info
    adata.write_h5ad(cache, compression="gzip")
    return adata, info


def annotate(adata) -> pd.DataFrame:
    import scanpy as sc
    sets = {}
    for ct, genes in A5.MARKERS.items():
        g = [x for x in genes if x in adata.var_names]
        if len(g) >= 2:
            sc.tl.score_genes(adata, g, score_name=f"mk_{ct}", random_state=C.SEED)
            sets[ct] = g
    cols = [f"mk_{ct}" for ct in sets]
    M = adata.obs.groupby("cluster", observed=True)[cols].mean()
    Z = (M - M.mean()) / M.std(ddof=0).replace(0, 1)
    lab = Z.idxmax(axis=1).str.replace("mk_", "", regex=False)
    adata.obs["cell_type"] = adata.obs["cluster"].map(lab).astype(str)
    out = Z.copy()
    out.columns = [c.replace("mk_", "") for c in out.columns]
    out.insert(0, "label", lab)
    out.insert(1, "n_cells", adata.obs["cluster"].value_counts().reindex(out.index).values)
    return out


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> float:
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))
    Z = np.column_stack([np.ones_like(rz), rz])
    ex = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]
    ey = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1])


# -------------------------------------------------------------------------------------
# figure 7
# -------------------------------------------------------------------------------------
def plot_fig7(obs: pd.DataFrame, umap: np.ndarray, delta: pd.DataFrame, coup: pd.DataFrame, donor: pd.DataFrame,
              pooled: dict) -> None:
    set_style()
    fig = plt.figure(figsize=(FIG_W, 165 * MM))
    gs = fig.add_gridspec(2, 3, hspace=0.38, wspace=0.55, left=0.06, right=0.97, top=0.95, bottom=0.09)
    rng = np.random.default_rng(C.SEED)
    idx = rng.permutation(len(obs))[:min(len(obs), 40000)]
    U, O = umap[idx], obs.iloc[idx]
    cts = list(A5.MARKERS)
    ax = fig.add_subplot(gs[0, 0])
    for i, ct in enumerate(cts):
        m = (O["cell_type"] == ct).to_numpy()
        if m.any():
            ax.scatter(U[m, 0], U[m, 1], s=0.3, color=A5.CT_COLORS[i % len(A5.CT_COLORS)], rasterized=True,
                       linewidths=0)
            cx, cy = np.median(U[m, 0]), np.median(U[m, 1])
            ax.text(cx, cy, ct.replace(" cells", "").replace(" monocytes", " mono"), fontsize=4.5, ha="center",
                    va="center", bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.7))
    ax.set_title(f"GSE167363 PBMC ({len(obs):,} cells)", fontsize=6)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    panel_label(ax, "a", x=-0.12)
    for k, (col, title, cmap, lab) in enumerate((("gc_score", "Glucocorticoid-response score", "magma", "b"),
                                                  ("mhc_score", "MHC class II score", "viridis", "c"))):
        ax = fig.add_subplot(gs[0, k + 1])
        v = O[col].to_numpy()
        lo, hi = np.nanpercentile(v, [2, 98])
        order = np.argsort(v)
        sc_ = ax.scatter(U[order, 0], U[order, 1], c=v[order], s=0.3, cmap=cmap, vmin=lo, vmax=hi, rasterized=True,
                         linewidths=0)
        cb = fig.colorbar(sc_, ax=ax, fraction=0.04, pad=0.02)
        cb.ax.tick_params(labelsize=4.5)
        ax.set_title(title, fontsize=6)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("UMAP 1")
        panel_label(ax, lab, x=-0.08)
    # d: sepsis minus healthy per cell type (donor-level)
    ax = fig.add_subplot(gs[1, 0])
    if len(delta):
        dd = delta.sort_values("delta_gc")
        y = np.arange(len(dd))
        ax.barh(y, dd["delta_gc"], color=[OKABE[1] if c == A5.MONO else "#BBBBBB" for c in dd.index], height=0.65)
        ax.set_yticks(y)
        ax.set_yticklabels(dd.index, fontsize=4.8)
        ax.axvline(0, color="#777777", lw=0.5)
    ax.set_xlabel("GC score: sepsis T0 - healthy\n(donor means)")
    ax.set_title("Where the programme rises", fontsize=6)
    panel_label(ax, "d", x=-0.55)
    # e: within-patient monocyte coupling
    ax = fig.add_subplot(gs[1, 1])
    if len(coup):
        rows = coup.copy()
        y = np.arange(len(rows) + 1)[::-1]
        for yy, (_, r) in zip(y[:-1], rows.iterrows()):
            col_ = {"non-survivor": OKABE[1], "survivor": OKABE[0], "healthy": "#999999"}[r["group"]]
            ax.plot([np.tanh(r["z"] - 1.96 * r["se"]), np.tanh(r["z"] + 1.96 * r["se"])], [yy, yy], color=col_, lw=0.8)
            ax.scatter([r["rho"]], [yy], color=col_, s=10, marker="s", zorder=3)
        if pooled:
            ax.plot([pooled["lo"], pooled["hi"]], [y[-1], y[-1]], color="black", lw=0.9)
            ax.scatter([pooled["rho"]], [y[-1]], color="black", marker="D", s=18, zorder=3)
        ax.set_yticks(y)
        tag = rows["group"].map({"healthy": "HC", "survivor": "S", "non-survivor": "NS"})
        ax.set_yticklabels([f"{d_.replace('donor', 'D')} {t_} ({g_})" for d_, t_, g_ in zip(rows["donor"], rows["time"], tag)]
                           + ["Pooled septic T0"],
                           fontsize=4.8)
        ax.axvline(0, color="#777777", lw=0.5, ls="--")
    ax.set_xlabel("Partial Spearman rho (GC vs MHC-II)")
    ax.set_title("Coupling within CD14 monocytes", fontsize=6)
    panel_label(ax, "e", x=-0.55)
    # f: donor-level monocyte scores by group
    ax = fig.add_subplot(gs[1, 2])
    if len(donor):
        for grp, col_, mk in (("healthy", "#999999", "o"), ("survivor", OKABE[0], "s"), ("non-survivor", OKABE[1], "^")):
            s = donor[(donor["group"] == grp) & (donor["time"] == "T0")]
            ax.scatter(s["mhc_score"], s["gc_score"], color=col_, marker=mk, s=22, label=grp, zorder=3)
            for _, r in s.iterrows():
                ax.annotate(r["donor"].replace("donor", "D"), (r["mhc_score"], r["gc_score"]), fontsize=4.5,
                            xytext=(3, 2), textcoords="offset points")
        ax.legend(frameon=False, fontsize=4.8)
    ax.set_xlabel("Monocyte MHC class II score (donor mean)")
    ax.set_ylabel("Monocyte GC score (donor mean)")
    ax.set_title("Donors at T0 (descriptive)", fontsize=6)
    panel_label(ax, "f", x=-0.3)
    save_figure(fig, "Fig7_single_cell_validation")


# -------------------------------------------------------------------------------------
# run
# -------------------------------------------------------------------------------------
def run_singlecell() -> None:
    import scanpy as sc
    t0 = time.time()
    warnings.filterwarnings("ignore")
    LOG.info(f"Project folder: {C.PROJECT_DIR}; scanpy {sc.__version__}")
    for f, w in ((C.TAB / "S27_gc_response_signature.tsv", "05_closing_analyses.py"),
                 (C.TAB / "S20_reversal_signature_prespecified.tsv", "03_analysis_pathways_LR.py")):
        if not f.exists():
            raise SystemExit(f"{f.name} missing - run {w} first.")
    write_preregistration()
    write_amendment()
    write_amendment_3()
    write_amendment_4()
    write_amendment_5()
    write_amendment_7()
    write_amendment_8()
    gates = []
    adata, info = build_atlas(sc_files())
    info = dict(adata.uns.get("processing", info))
    gates.append({"gate": "single-cell: cells after QC", "value": f"{adata.n_obs} ({info.get('integration')}, "
                  f"{info.get('clustering')})", "criterion": ">= 20000", "passed": bool(adata.n_obs >= 20000),
                  "severity": "critical"})
    if adata.n_obs < 20000:
        record_gates("singlecell", gates)
        raise SystemExit("Too few cells after QC - check the downloaded files.")
    ann = annotate(adata)
    ann.round(2).to_csv(C.TAB / "S33_sc_cluster_annotation.tsv", sep="\t")
    counts = pd.crosstab(adata.obs["sample"], adata.obs["cell_type"])
    counts.to_csv(C.TAB / "S33b_sc_cell_counts.tsv", sep="\t")
    LOG.info("Cell types per sample:\n" + counts.to_string())
    # scores
    sig = pd.read_csv(C.TAB / "S27_gc_response_signature.tsv", sep="\t")
    s20 = pd.read_csv(C.TAB / "S20_reversal_signature_prespecified.tsv", sep="\t")
    cov = {}

    def score(name, up, dn=None):
        u = [g for g in up if g in adata.var_names]
        d = [g for g in (dn or []) if g in adata.var_names]
        cov[name] = (len(u), len(up), len(d), len(dn or []))
        sc.tl.score_genes(adata, u, score_name="_u", random_state=C.SEED)
        v = adata.obs["_u"].to_numpy()
        if d:
            sc.tl.score_genes(adata, d, score_name="_d", random_state=C.SEED)
            v = v - adata.obs["_d"].to_numpy()
        adata.obs[name] = v
    score("gc_score", sig.loc[sig["direction"] == "UP", "gene"].tolist(), sig.loc[sig["direction"] == "DOWN", "gene"].tolist())
    score("mhc_score", C.GENE_SETS["mhc_class_ii"])
    score("s20_score", s20.loc[s20["direction"] == "UP", "gene"].tolist(), s20.loc[s20["direction"] == "DOWN", "gene"].tolist())
    LOG.info("Score genes found (UP found/total, DOWN found/total): " + "; ".join(f"{k} {v}" for k, v in cov.items()))
    gc_cov = (cov["gc_score"][0] + cov["gc_score"][2]) / max(cov["gc_score"][1] + cov["gc_score"][3], 1)
    gates.append({"gate": "GC-response genes detected in the single-cell data", "value": f"{gc_cov:.0%}",
                  "criterion": ">= 80%", "passed": bool(gc_cov >= 0.8), "severity": "warning"})
    adata.obs["log_counts"] = np.log1p(adata.obs["total_counts"].astype(float))
    obs = adata.obs
    # (sc1) within-patient coupling in CD14 monocytes
    rows = []
    for smp, g in obs[obs["cell_type"] == A5.MONO].groupby("sample", observed=True):
        if len(g) < A5.MIN_MONO:
            LOG.warning(f"{smp}: {len(g)} CD14 monocytes (< {A5.MIN_MONO}) - not used for the coupling test")
            continue
        rho = partial_spearman(g["gc_score"].to_numpy(), g["mhc_score"].to_numpy(), g["log_counts"].to_numpy())
        rows.append({"sample": smp, "donor": g["donor"].iloc[0], "group": g["group"].iloc[0], "time": g["time"].iloc[0],
                     "n_cells": len(g), "rho": rho, "z": float(np.arctanh(np.clip(rho, -0.999, 0.999))),
                     "se": 1 / np.sqrt(len(g) - 3)})
    coup = pd.DataFrame(rows)
    coup.to_csv(C.TAB / "S35_sc_monocyte_coupling.tsv", sep="\t", index=False)
    pooled = {}
    for tp in ("T0", "T6"):
        s = coup[(coup["group"] != "healthy") & (coup["time"] == tp)] if len(coup) else coup
        if len(s) >= 2:
            m = dl_meta(s["z"].to_numpy(), s["se"].to_numpy())
            res = {"rho": float(np.tanh(m["est"])), "lo": float(np.tanh(m["est"] - 1.96 * m["se"])),
                   "hi": float(np.tanh(m["est"] + 1.96 * m["se"])), "P": m["P"], "k": m["k"],
                   "n_negative": int((s["rho"] < 0).sum())}
            if tp == "T0":
                pooled = res
            supported = res["rho"] < 0 and res["P"] < 0.05 and res["n_negative"] >= max(len(s) - 1, 1)
            res["P_sign"] = float(stats.binomtest(res["n_negative"], len(s), 0.5).pvalue)
            LOG.info(f"(sc1{' sensitivity' if tp == 'T6' else ''}) {tp}: pooled partial rho {res['rho']:.2f} "
                     f"({res['lo']:.2f} to {res['hi']:.2f}), P = {max(res['P'], 1e-300):.2g}; {res['n_negative']}/{len(s)} "
                     f"patients negative (patient-level sign test P = {res['P_sign']:.2g}) -> "
                     f"{'supported' if supported else 'not supported'}")
            gates.append({"gate": f"{'PRIMARY (sc1)' if tp == 'T0' else 'Sensitivity (sc1, T6)'}: GC vs MHC-II "
                                  f"coupling in patients' CD14 monocytes",
                          "value": f"rho {res['rho']:.2f} ({res['lo']:.2f} to {res['hi']:.2f}), P = "
                                   f"{max(res['P'], 1e-300):.2g}, {res['n_negative']}/{len(s)} negative (sign test P = "
                                   f"{res['P_sign']:.2g})",
                          "criterion": "rho < 0, P < 0.05, >= 4 of 5 patients negative", "passed": bool(supported),
                          "severity": "warning"})
    hc = coup[coup["group"] == "healthy"] if len(coup) else coup
    if len(hc):
        LOG.info("Healthy donors (reference): " + ", ".join(f"{r.sample} rho {r.rho:.2f}" for r in hc.itertuples()))
    # donor-level means (sc2, sc3)
    donor = obs.groupby(["donor", "group", "time", "cell_type"], observed=True)[["gc_score", "mhc_score", "s20_score"]] \
        .mean().reset_index()
    donor.to_csv(C.TAB / "S34_sc_donor_celltype_scores.tsv", sep="\t", index=False)
    t0d = donor[donor["time"] == "T0"]
    sep_ = t0d[t0d["group"] != "healthy"].groupby("cell_type")[["gc_score", "mhc_score", "s20_score"]].mean()
    hea_ = t0d[t0d["group"] == "healthy"].groupby("cell_type")[["gc_score", "mhc_score", "s20_score"]].mean()
    n_ct = obs.groupby("cell_type", observed=True).size()
    delta = (sep_ - hea_).dropna().rename(columns=lambda c: "delta_" + c.replace("_score", ""))
    delta = delta[delta.index.map(lambda c: n_ct.get(c, 0) >= 200)]
    delta.to_csv(C.TAB / "S34b_sc_sepsis_minus_healthy.tsv", sep="\t")
    if len(delta):
        top = delta["delta_gc"].idxmax()
        rank = int(delta["delta_gc"].rank(ascending=False)[A5.MONO]) if A5.MONO in delta.index else None
        LOG.info("(sc2) GC score, septic T0 minus healthy (donor means), by cell type:\n" +
                 delta.sort_values("delta_gc", ascending=False).round(3).to_string())
        gates.append({"gate": "(sc2, descriptive) cell type with the largest GC-score rise in sepsis",
                      "value": f"{top}; CD14 monocytes rank {rank} of {len(delta)}", "criterion": "reported",
                      "passed": True, "severity": "warning"})
    s20_ct = obs[obs["group"] != "healthy"].groupby("cell_type", observed=True)["s20_score"].mean().sort_values(
        ascending=False)
    LOG.info("(sc3) death-age signature score by cell type (septic patients):\n" + s20_ct.round(3).to_string())
    mono_d = donor[donor["cell_type"] == A5.MONO]
    LOG.info("(sc3) CD14 monocytes, donor means at T0:\n" + mono_d[mono_d["time"] == "T0"][
        ["donor", "group", "gc_score", "mhc_score"]].round(3).to_string(index=False))
    record_gates("singlecell", gates)
    obs[["sample", "donor", "group", "time", "cell_type", "cluster", "gc_score", "mhc_score", "s20_score",
         "total_counts", "n_genes_by_counts", "pct_counts_mt"]].to_csv(C.TAB / "S36_sc_cell_scores.tsv.gz", sep="\t")
    safe_plot(plot_fig7, obs.reset_index(drop=True), np.asarray(adata.obsm["X_umap"]), delta,
              coup.sort_values(["time", "group", "sample"]).reset_index(drop=True) if len(coup) else coup,
              mono_d, pooled)
    LOG.info(f"Run time {(time.time() - t0) / 60:.1f} min")
    LOG.info(gate_summary().replace("STAGE 1 QUALITY GATES", "STAGE 06 - SINGLE-CELL VALIDATION").replace(
        "Stage 2 can start after the field audit has been checked.", "Stage 06 complete.").replace(
        "do not start Stage 2", "check the failed gates"))


# #####################################################################################
# PART G - FIGURE POLISH (Fig 4a, Fig 5b, Fig 7a/f): plots ONLY from the saved result tables; no statistic is
#          recomputed except the pooled single-cell estimate of Fig 7e, which is recomputed from S35 exactly as in 06
# #####################################################################################
LOG = get_logger("sepsis_figures")
GATES_FILE = C.TAB / "S0_figures_polish.tsv"
_GATES.clear()


def repel_labels(fig, ax, xs, ys, texts, fontsize=4.6, color="black", line_color="#777777", max_r=10,
                 obstacles=None, avoid=(), box=False, start=1):
    """Greedy non-overlapping placement: each label tries growing rings of offsets around its point and is placed
    where its box overlaps no earlier label, no other point and stays inside the axes; a thin leader line links
    label and point."""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    axbox = ax.get_window_extent(rend)
    pts = ax.transData.transform(np.column_stack([xs, ys]))
    obs_pts = ax.transData.transform(np.asarray(obstacles, float)) if obstacles is not None else pts
    boxes = [a_.get_window_extent(rend).expanded(1.05, 1.1) for a_ in avoid]
    dirs = [(1, 1), (1, -1), (-1, 1), (-1, -1), (1, 0), (-1, 0), (0, 1), (0, -1)]
    centre = obs_pts.mean(axis=0)
    own = ax.transData.transform(np.column_stack([xs, ys]))
    for k, (x, y, t) in enumerate(zip(xs, ys, texts)):
        placed, best = None, None
        v = own[k] - centre                                  # prefer directions pointing away from the data cloud
        order = sorted(dirs, key=lambda d: -(d[0] * v[0] + d[1] * v[1]) / (np.hypot(*d) * (np.hypot(*v) + 1e-9)))
        for r in range(start, max_r + 1):
            for dx, dy in order:
                off = (dx * (4 + 7 * r), dy * (3 + 6 * r))
                kw = dict(fontsize=fontsize, color=color, textcoords="offset points", zorder=5,
                          bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85) if box else None,
                          ha="left" if dx > 0 else ("right" if dx < 0 else "center"),
                          va="bottom" if dy > 0 else ("top" if dy < 0 else "center"))
                probe = ax.annotate(t, (x, y), xytext=off, **kw)          # text only (no leader line) to measure
                bb = probe.get_window_extent(rend).expanded(1.06, 1.15)
                probe.remove()
                inside = (bb.x0 >= axbox.x0 and bb.x1 <= axbox.x1 and bb.y0 >= axbox.y0 and bb.y1 <= axbox.y1)
                near = (obs_pts[:, 0] >= bb.x0 - 3) & (obs_pts[:, 0] <= bb.x1 + 3) &                        (obs_pts[:, 1] >= bb.y0 - 3) & (obs_pts[:, 1] <= bb.y1 + 3)
                clear = not any(bb.overlaps(b) for b in boxes)
                if inside and clear and not near.any():
                    placed = (off, kw, bb)
                    break
                if inside and clear and best is None:
                    best = (off, kw, bb)                                # fallback: no label overlap, may cover a dot
            if placed is not None:
                break
        placed = placed or best
        if placed is None:
            off, kw = (4, 3), dict(fontsize=fontsize, color=color, textcoords="offset points", zorder=5)
            bb = None
        else:
            off, kw, bb = placed
        ann = ax.annotate(t, (x, y), xytext=off, arrowprops=dict(arrowstyle="-", color=line_color, lw=0.35,
                                                                   shrinkA=0, shrinkB=1.5), **kw)
        boxes.append(bb if bb is not None else ann.get_window_extent(rend))


def plot_fig4_v2(hall_d: pd.DataFrame, hall_a: pd.DataFrame, react_d: pd.DataFrame, lr: pd.DataFrame) -> None:
    set_style()
    fig = plt.figure(figsize=(FIG_W, 125 * MM))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1], height_ratios=[1, 0.85], wspace=1.45, hspace=0.6,
                          left=0.2, right=0.95, top=0.93, bottom=0.08)
    # a: Hallmark pathways - death and age side by side (colour = NES, size = -log10 FDR, ring = FDR < 0.05)
    ax = fig.add_subplot(gs[0, 0])
    h = hall_d.merge(hall_a[["pathway", "NES", "FDR"]], on="pathway", suffixes=("", "_age"))
    h = h[(h["FDR"] < A2.FDR) | (h["FDR_age"] < A2.FDR)]
    h = pd.concat([h.nsmallest(10, "NES"), h.nlargest(14, "NES")]).drop_duplicates("pathway").sort_values("NES")
    y = np.arange(len(h))
    lim = float(np.ceil(np.nanmax(np.abs(h[["NES", "NES_age"]].to_numpy()))))
    size = lambda f: 5 + 20 * np.clip(-np.log10(np.clip(f.to_numpy(float), 1e-4, 1)), 0, 4) / 4
    for xi, (nes, fdr) in enumerate((("NES", "FDR"), ("NES_age", "FDR_age"))):
        sig = (h[fdr] < A2.FDR).to_numpy()
        sc_ = ax.scatter(np.full(len(h), xi), y, s=size(h[fdr]), c=h[nes], cmap="RdBu_r", vmin=-lim, vmax=lim,
                         edgecolors=np.where(sig, "black", "none"), linewidths=0.6, zorder=3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Death", "Age"], fontsize=5.5)
    ax.set_xlim(-0.7, 1.7)
    ax.set_ylim(-0.8, len(h) - 0.2)
    for yy in y:
        ax.axhline(yy, color="#EEEEEE", lw=0.4, zorder=0)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(axis="x", length=0)
    ax.set_yticks(y)
    ax.set_yticklabels([(p[0].upper() + p[1:]).replace("_", " ")[:40] for p in h["pathway"]], fontsize=5)
    pos = ax.get_position()
    ax.set_position([pos.x0, pos.y0, pos.width * 0.34, pos.height])
    pos = ax.get_position()
    cax = fig.add_axes([pos.x1 + 0.012, pos.y0 + pos.height * 0.58, 0.008, pos.height * 0.38])
    cb = fig.colorbar(sc_, cax=cax)
    cb.set_label("NES", fontsize=5)
    cb.ax.tick_params(labelsize=4.5)
    from matplotlib.lines import Line2D
    hs = [Line2D([], [], marker="o", ls="", markerfacecolor="#BBBBBB", markeredgecolor="none",
                 markersize=np.sqrt(size(pd.Series([f]))[0]), label=l_) for f, l_ in ((0.05, "0.05"), (0.01, "0.01"),
                                                                                    (1e-4, "<= 1e-4"))]
    hs.append(Line2D([], [], marker="o", ls="", markerfacecolor="white", markeredgecolor="black", markersize=4,
                     label="FDR < 0.05"))
    ax.legend(handles=hs, title="FDR", title_fontsize=4.8, fontsize=4.6, frameon=False, loc="upper left",
              bbox_to_anchor=(1.05, 0.5), handletextpad=0.3, labelspacing=0.5)
    ax.set_title("Hallmark pathways:\ndeath and age", fontsize=6)
    panel_label(ax, "a", x=-0.9 / 0.34, y=1.03)
    # b: top Reactome death pathways
    ax = fig.add_subplot(gs[0, 1])
    r = react_d[react_d["FDR"] < A2.FDR].copy()
    r = pd.concat([r.nlargest(10, "NES"), r.nsmallest(10, "NES")]).drop_duplicates("pathway").sort_values("NES")
    y = np.arange(len(r))
    sc_ = ax.scatter(r["NES"], y, s=np.clip(r["size"] / 4, 6, 40), c=-np.log10(r["FDR"].clip(lower=1e-10)),
                     cmap="viridis", linewidths=0)
    ax.axvline(0, color="#777777", lw=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([re.sub(r"\s*R-HSA-\d+", "", p)[:46] for p in r["pathway"]], fontsize=4.8)
    ax.set_xlabel("NES (death)")
    ax.set_title("Reactome pathways and death\n(FDR < 0.05; top 10 per direction)", fontsize=6)
    cb = fig.colorbar(sc_, ax=ax, fraction=0.05, pad=0.03)
    cb.set_label("-log10 FDR", fontsize=5)
    cb.ax.tick_params(labelsize=5)
    panel_label(ax, "b", x=-1.35, y=1.03)
    # c, d: ligand-receptor axes
    l_ = lr.sort_values("death_OR")
    y = np.arange(len(l_))
    ax = fig.add_subplot(gs[1, 0])
    for off, col_, est, lo, hi, lab in ((0.12, OKABE[0], "death_OR", "death_lo", "death_hi", "age, sex adjusted"),
                                        (-0.12, OKABE[1], "death_OR_adjML", "adjML_lo", "adjML_hi", "+ ML index")):
        ax.errorbar(l_[est], y + off, xerr=[l_[est] - l_[lo], l_[hi] - l_[est]], fmt="o", ms=2.5, color=col_,
                    elinewidth=0.6, capsize=0, label=lab)
    ax.axvline(1, color="#777777", lw=0.5, ls="--")
    log_axis(ax)
    ax.set_yticks(y)
    ax.set_yticklabels(l_["axis"], fontsize=5)
    ax.set_xlabel("Pooled OR per SD (all cohorts)")
    ax.set_title("Ligand-receptor axes and death", fontsize=6)
    ax.legend(frameon=False, fontsize=5, loc="upper center", bbox_to_anchor=(0.45, -0.22), ncol=2)
    panel_label(ax, "c", x=-0.9, y=1.03)
    ax = fig.add_subplot(gs[1, 1])
    ax.errorbar(l_["age_b"], y, xerr=1.96 * l_["age_se"], fmt="o", ms=2.5, color=OKABE[2], elinewidth=0.6, capsize=0)
    ax.axvline(0, color="#777777", lw=0.5, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_xlabel("SD change per 10 years of age")
    ax.set_title("Ligand-receptor axes and age", fontsize=6)
    panel_label(ax, "d", x=-0.12, y=1.03)
    save_figure(fig, "Fig4_pathways_ligand_receptor")


def plot_fig5_v2(genes_tab: pd.DataFrame, res: pd.DataFrame, cand: pd.DataFrame, pc: pd.DataFrame, dg: pd.DataFrame,
              primary: str) -> None:
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    set_style()
    fig = plt.figure(figsize=(FIG_W, 140 * MM))
    outer = fig.add_gridspec(2, 1, height_ratios=[1, 1.1], hspace=0.55, left=0.08, right=0.98, top=0.95, bottom=0.08)
    top_ = outer[0].subgridspec(1, 2, width_ratios=[1, 1.9], wspace=0.35)
    bot_ = outer[1].subgridspec(1, 3, width_ratios=[1.05, 1, 1], wspace=1.25)
    # a: which signature genes are testable in PBMCs
    ax = fig.add_subplot(top_[0, 0])
    cats = [("in OP3 (testable)", OKABE[2]), ("erythroid (absent)", "#E8A0A0"), ("neutrophil (absent)", "#F2C879"),
            ("other (absent)", "#CCCCCC")]
    for i, dirn in enumerate(("UP", "DOWN")):
        s = genes_tab[genes_tab["direction"] == dirn]
        n = len(s)
        vals = [s["in_OP3"].sum(), ((~s["in_OP3"]) & (s["lineage"] == "erythroid")).sum(),
                ((~s["in_OP3"]) & (s["lineage"] == "neutrophil")).sum(), ((~s["in_OP3"]) & (s["lineage"] == "")).sum()]
        left = 0
        for v, (_, col) in zip(vals, cats):
            ax.barh(1 - i, v / n, left=left, color=col, height=0.55)
            left += v / n
        ax.text(1.02, 1 - i, f"{vals[0]}/{n}", va="center", fontsize=5)
    ax.axvline(A3.FULL_COVERAGE, color="#555555", lw=0.5, ls="--")
    ax.set_yticks([1, 0])
    ax.set_yticklabels(["UP in non-survivors", "DOWN in non-survivors"])
    ax.set_xlim(0, 1.18)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xlabel("Fraction of signature genes")
    ax.legend(handles=[Patch(color=c_, label=l_) for l_, c_ in cats], frameon=False, fontsize=4.5, ncol=2,
              loc="upper center", bbox_to_anchor=(0.45, -0.25))
    ax.set_title("Testable part of the death signature\n(OP3 PBMC gene universe)", fontsize=6)
    panel_label(ax, "a", x=-0.45)
    # b: primary cell type, reversal vs MHC-II (non-overlapping labels)
    ax = fig.add_subplot(top_[0, 1])
    p = res[res["cell_type"] == primary].copy()
    if len(p) and "mhc_z" in p:
        has_r = p["R_z"].notna().any()
        xv = p["R_z"] if has_r else pd.Series(0.0, index=p.index)
        cs = cand.set_index("compound") if len(cand) else pd.DataFrame(columns=["tier"])
        tier = p["sm_name"].map(cs["tier"]) if len(cand) else pd.Series(np.nan, index=p.index)
        colr = np.array(["#C8C8C8"] * len(p), dtype=object)
        colr[(tier == "Tier 2").to_numpy()] = OKABE[4]
        colr[(tier == "Tier 1").to_numpy()] = OKABE[2]
        colr[(tier == "Caution").to_numpy()] = OKABE[1]
        is_gc = p["sm_name"].str.contains(A3.GLUCOCORTICOIDS, case=False, regex=True).to_numpy()
        xl = np.nanmax(np.abs(xv)) * 1.12
        yl = np.nanmax(np.abs(p["mhc_z"])) * 1.12
        ax.add_patch(plt.Rectangle((-xl, 0), xl, yl, color=OKABE[2], alpha=0.06, lw=0, zorder=0))
        ax.add_patch(plt.Rectangle((0, -yl), xl, yl, color=OKABE[1], alpha=0.06, lw=0, zorder=0))
        q1 = ax.text(-xl * 0.97, yl * 0.95, "reverses signature\n+ restores MHC-II", fontsize=4.6, color=OKABE[2],
                     va="top", ha="left", style="italic")
        q2 = ax.text(xl * 0.97, -yl * 0.95, "mimics signature\n+ lowers MHC-II", fontsize=4.6, color=OKABE[1],
                     va="bottom", ha="right", style="italic")
        ax.scatter(xv, p["mhc_z"], s=np.where(colr == "#C8C8C8", 7, 14), c=list(colr), edgecolor="white",
                   linewidth=0.3, zorder=3)
        ax.scatter(xv[is_gc], p.loc[is_gc, "mhc_z"], s=26, facecolors="none", edgecolors="black", linewidths=0.7,
                   zorder=4)
        ax.axhline(0, color="#777777", lw=0.5)
        ax.axvline(0, color="#777777", lw=0.5)
        ax.set_xlim(-xl, xl)
        ax.set_ylim(-yl, yl)
        # labels: candidates relevant to the primary cell type, glucocorticoids, and the 3 strongest MHC-II lowerers
        rel = p["sm_name"].map(lambda c: bool(c in cs.index and (str(cs.loc[c].get("reversal_primary")) == "True" or
                                                                  str(cs.loc[c].get("mhc_up_primary")) == "True")))
        caut = p[(tier == "Caution").to_numpy()].reindex(p.loc[(tier == "Caution").to_numpy(), "mhc_z"].sort_values()
                                                          .index).head(3).index
        show = p.index[rel.to_numpy() | is_gc] .union(caut)
        show = (p.loc[show, "mhc_z"].abs() + xv.loc[show].abs()).sort_values(ascending=False).index
        names = [str(p.loc[i_, "sm_name"]).split(" (")[0].split(";")[0][:20] for i_ in show]
        allpts = np.column_stack([xv.to_numpy(float), p["mhc_z"].to_numpy(float)])
        allpts = allpts[np.isfinite(allpts).all(axis=1)]
        repel_labels(fig, ax, xv.loc[show].to_numpy(float), p.loc[show, "mhc_z"].to_numpy(float), names,
                     obstacles=allpts, avoid=(q1, q2))
        ax.set_xlabel("Signature reversal z (< 0: reverses)" if has_r else "Signature reversal not testable (gate)")
        ax.set_ylabel("MHC class II shift z (> 0: restores)")
        hs = [Line2D([], [], marker="o", ls="", color=c_, markersize=3, label=l_) for c_, l_ in
              ((OKABE[2], "Tier 1"), (OKABE[4], "Tier 2"), (OKABE[1], "Caution"), ("#C8C8C8", "Other"))]
        hs.append(Line2D([], [], marker="o", ls="", markerfacecolor="none", markeredgecolor="black", markersize=4,
                         label="Glucocorticoid"))
        ax.legend(handles=hs, frameon=False, fontsize=4.5, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    ax.set_title(f"{primary}: {len(p)} compounds (OP3, 24 h)", fontsize=6)
    panel_label(ax, "b", x=-0.12)
    # c: candidates x cell types
    ax = fig.add_subplot(bot_[0, 0])
    cts = sorted(res["cell_type"].unique())
    if res["R_z"].notna().any():
        top = (cand.loc[cand["tier"] != "Caution", "compound"].head(15).tolist() if len(cand) else []) or \
            res[res["cell_type"] == primary].nsmallest(15, "R_z")["sm_name"].tolist()
        Mz = res.pivot_table(index="sm_name", columns="cell_type", values="R_z").reindex(index=top, columns=cts)
        Fz = res.pivot_table(index="sm_name", columns="cell_type", values="R_FDR").reindex(index=top, columns=cts)
        lim = np.nanmax(np.abs(Mz.to_numpy())) if np.isfinite(Mz.to_numpy()).any() else 1
        im = ax.imshow(Mz.to_numpy(), cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto", interpolation="none")
        for (i, j), v in np.ndenumerate(Fz.to_numpy()):
            if np.isfinite(v) and v < A3.FDR:
                ax.text(j, i, "*", ha="center", va="center", fontsize=6)
        ax.set_xticks(range(len(cts)))
        ax.set_xticklabels([short_ct(c) for c in cts], rotation=45, ha="right")
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels([t[:24] for t in top], fontsize=4.5)
        cb = fig.colorbar(im, ax=ax, fraction=0.06, pad=0.04, shrink=0.8)
        cb.set_label("Reversal z", fontsize=5)
        cb.ax.tick_params(labelsize=4.5)
    else:
        ax.text(0.5, 0.5, "Reversal not testable\n(expression gate)", ha="center", va="center", fontsize=6)
        ax.axis("off")
    ax.set_title("Reversal across PBMC cell types\n(* FDR < 0.05)", fontsize=6)
    panel_label(ax, "c", x=-0.75)
    # d: positive controls
    ax = fig.add_subplot(bot_[0, 1])
    if len(pc):
        yy = np.arange(len(pc))[::-1]
        ax.barh(yy, pc["value"], color=[OKABE[2] if ok else "#BBBBBB" for ok in pc["passed"]], height=0.6)
        ax.set_yticks(yy)
        ax.set_yticklabels([s[:34] for s in pc["label"]], fontsize=4.5)
        ax.axvline(0, color="#777777", lw=0.5)
    else:
        ax.text(0.5, 0.5, "No control testable", ha="center", va="center", fontsize=6)
    ax.set_xlabel("Effect (n genes, z or -log10 adj. P; S25)")
    ax.set_title("Positive controls\n(green: as expected)", fontsize=6)
    panel_label(ax, "d", x=-1.05)
    # e: DGIdb
    ax = fig.add_subplot(bot_[0, 2])
    if len(dg):
        t = dg[dg["approved"]].groupby("gene").agg(n=("drug", "nunique"), match=("matches_desired_direction", "sum"))
        order = [g for g in list(A3.LR_TARGETS) + list(A3.CONTROL_RECEPTORS) if g in t.index]
        t = t.reindex(order)
        yy = np.arange(len(t))[::-1]
        ax.barh(yy, t["n"], color="#CCCCCC", height=0.6, label="approved drugs")
        ax.barh(yy, t["match"], color=OKABE[0], height=0.6, label="desired direction")
        ax.set_yticks(yy)
        ax.set_yticklabels([f"{g} ({(A3.LR_TARGETS | A3.CONTROL_RECEPTORS)[g][:3]}.)" for g in t.index], fontsize=4.5)
        ax.legend(frameon=False, fontsize=4.5, loc="upper center", bbox_to_anchor=(0.4, -0.16), ncol=2)
    else:
        ax.text(0.5, 0.5, "DGIdb not reachable", ha="center", va="center", fontsize=6)
    ax.set_xlabel("Approved drugs (DGIdb)")
    ax.set_title("Druggable axis genes", fontsize=6)
    panel_label(ax, "e", x=-0.75)
    save_figure(fig, "Fig5_drug_repurposing")


def plot_fig7_v2(obs: pd.DataFrame, umap: np.ndarray, delta: pd.DataFrame, coup: pd.DataFrame, donor: pd.DataFrame,
              pooled: dict) -> None:
    set_style()
    fig = plt.figure(figsize=(FIG_W, 165 * MM))
    gs = fig.add_gridspec(2, 3, hspace=0.38, wspace=0.55, left=0.06, right=0.97, top=0.95, bottom=0.09)
    rng = np.random.default_rng(C.SEED)
    idx = rng.permutation(len(obs))[:min(len(obs), 40000)]
    U, O = umap[idx], obs.iloc[idx]
    cts = list(A5.MARKERS)
    lab_x, lab_y, lab_t = [], [], []
    ax = fig.add_subplot(gs[0, 0])
    for i, ct in enumerate(cts):
        m = (O["cell_type"] == ct).to_numpy()
        if m.any():
            ax.scatter(U[m, 0], U[m, 1], s=0.3, color=A5.CT_COLORS[i % len(A5.CT_COLORS)], rasterized=True,
                       linewidths=0)
            lab_x.append(np.median(U[m, 0]))
            lab_y.append(np.median(U[m, 1]))
            lab_t.append(ct.replace(" cells", "").replace(" monocytes", " mono"))
    repel_labels(fig, ax, np.array(lab_x), np.array(lab_y), lab_t, fontsize=4.8, line_color="#444444", box=True,
                 start=0, max_r=6)
    ax.set_title(f"GSE167363 PBMC ({len(obs):,} cells)", fontsize=6)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    panel_label(ax, "a", x=-0.12)
    for k, (col, title, cmap, lab) in enumerate((("gc_score", "Glucocorticoid-response score", "magma", "b"),
                                                  ("mhc_score", "MHC class II score", "viridis", "c"))):
        ax = fig.add_subplot(gs[0, k + 1])
        v = O[col].to_numpy()
        lo, hi = np.nanpercentile(v, [2, 98])
        order = np.argsort(v)
        sc_ = ax.scatter(U[order, 0], U[order, 1], c=v[order], s=0.3, cmap=cmap, vmin=lo, vmax=hi, rasterized=True,
                         linewidths=0)
        cb = fig.colorbar(sc_, ax=ax, fraction=0.04, pad=0.02)
        cb.ax.tick_params(labelsize=4.5)
        ax.set_title(title, fontsize=6)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("UMAP 1")
        panel_label(ax, lab, x=-0.08)
    # d: sepsis minus healthy per cell type (donor-level)
    ax = fig.add_subplot(gs[1, 0])
    if len(delta):
        dd = delta.sort_values("delta_gc")
        y = np.arange(len(dd))
        ax.barh(y, dd["delta_gc"], color=[OKABE[1] if c == A5.MONO else "#BBBBBB" for c in dd.index], height=0.65)
        ax.set_yticks(y)
        ax.set_yticklabels(dd.index, fontsize=4.8)
        ax.axvline(0, color="#777777", lw=0.5)
    ax.set_xlabel("GC score: sepsis T0 - healthy\n(donor means)")
    ax.set_title("Where the programme rises", fontsize=6)
    panel_label(ax, "d", x=-0.55)
    # e: within-patient monocyte coupling
    ax = fig.add_subplot(gs[1, 1])
    if len(coup):
        rows = coup.copy()
        y = np.arange(len(rows) + 1)[::-1]
        for yy, (_, r) in zip(y[:-1], rows.iterrows()):
            col_ = {"non-survivor": OKABE[1], "survivor": OKABE[0], "healthy": "#999999"}[r["group"]]
            ax.plot([np.tanh(r["z"] - 1.96 * r["se"]), np.tanh(r["z"] + 1.96 * r["se"])], [yy, yy], color=col_, lw=0.8)
            ax.scatter([r["rho"]], [yy], color=col_, s=10, marker="s", zorder=3)
        if pooled:
            ax.plot([pooled["lo"], pooled["hi"]], [y[-1], y[-1]], color="black", lw=0.9)
            ax.scatter([pooled["rho"]], [y[-1]], color="black", marker="D", s=18, zorder=3)
        ax.set_yticks(y)
        tag = rows["group"].map({"healthy": "HC", "survivor": "S", "non-survivor": "NS"})
        ax.set_yticklabels([f"{d_.replace('donor', 'D')} {t_} ({g_})" for d_, t_, g_ in zip(rows["donor"], rows["time"], tag)]
                           + ["Pooled septic T0"],
                           fontsize=4.8)
        ax.axvline(0, color="#777777", lw=0.5, ls="--")
    ax.set_xlabel("Partial Spearman rho (GC vs MHC-II)")
    ax.set_title("Coupling within CD14 monocytes", fontsize=6)
    panel_label(ax, "e", x=-0.55)
    # f: donor-level CD14-monocyte scores by outcome group, T0 -> T6 (descriptive)
    sub = gs[1, 2].subgridspec(1, 2, wspace=0.75)
    groups = [("healthy", "Healthy", "#8C8C8C", "o"), ("survivor", "Survivor", OKABE[0], "s"),
              ("non-survivor", "Non-survivor", OKABE[1], "^")]
    for j, (col, ttl) in enumerate((("mhc_score", "MHC class II"), ("gc_score", "GC response"))):
        ax = fig.add_subplot(sub[0, j])
        if len(donor):
            for gx, (grp, glab, col_, mk) in enumerate(groups):
                s = donor[donor["group"] == grp]
                donors_ = sorted(s["donor"].unique())
                for k, d_ in enumerate(donors_):
                    sd = s[s["donor"] == d_].set_index("time")
                    x0 = gx + (k - (len(donors_) - 1) / 2) * 0.22 - 0.05
                    if "T0" in sd.index and "T6" in sd.index:
                        ax.plot([x0, x0 + 0.1], [sd.loc["T0", col], sd.loc["T6", col]], color=col_, lw=0.6, alpha=0.7)
                    if "T0" in sd.index:
                        ax.scatter([x0], [sd.loc["T0", col]], color=col_, marker=mk, s=16, zorder=3,
                                   edgecolors="white", linewidths=0.3)
                    if "T6" in sd.index:
                        ax.scatter([x0 + 0.1], [sd.loc["T6", col]], facecolors="white", edgecolors=col_, marker=mk,
                                   s=12, linewidths=0.7, zorder=3)
                t0v = s.loc[s["time"] == "T0", col]
                if len(t0v):
                    ax.plot([gx - 0.3, gx + 0.3], [t0v.mean()] * 2, color=col_, lw=1.0, zorder=2)
            ax.set_xticks(range(len(groups)))
            ax.set_xticklabels([g[1] for g in groups], rotation=40, ha="right", fontsize=5)
            ax.set_xlim(-0.6, len(groups) - 0.4)
        ax.set_title(ttl, fontsize=5.8)
        ax.set_ylabel("CD14 monocytes, donor mean" if j == 0 else "", fontsize=5)
        if j == 0:
            panel_label(ax, "f", x=-0.75)
    from matplotlib.lines import Line2D
    hs = [Line2D([], [], marker="o", ls="", color="#555555", markersize=3, label="T0 (filled)"),
          Line2D([], [], marker="o", ls="", markerfacecolor="white", markeredgecolor="#555555", markersize=3,
                 label="T6 (open)"), Line2D([], [], color="#555555", lw=1.0, label="group mean, T0")]
    ax.legend(handles=hs, frameon=False, fontsize=4.4, loc="upper center", bbox_to_anchor=(-0.45, -0.24), ncol=3,
              handletextpad=0.3, columnspacing=0.8)
    save_figure(fig, "Fig7_single_cell_validation")


def run_figures() -> None:
    import anndata as ad
    t0 = time.time()
    warnings.filterwarnings("ignore")
    LOG.info(f"Project folder: {C.PROJECT_DIR}")
    bk = C.FIG / "before_polish"
    bk.mkdir(exist_ok=True)
    for stem in ("Fig4_pathways_ligand_receptor", "Fig5_drug_repurposing", "Fig7_single_cell_validation"):
        for ext in ("pdf", "png", "tiff"):
            f = C.FIG / f"{stem}.{ext}"
            if f.exists() and not (bk / f.name).exists():
                shutil.copy2(f, bk / f.name)
    LOG.info(f"Previous versions kept in {bk}")
    done = []
    # ---------------- Fig 4 ----------------
    if (C.TAB / "S19_gsea.tsv").exists() and (C.TAB / "S21_ligand_receptor_bulk.tsv").exists():
        s = pd.read_csv(C.TAB / "S19_gsea.tsv", sep="\t")
        lr = pd.read_csv(C.TAB / "S21_ligand_receptor_bulk.tsv", sep="\t")
        pick = lambda r, c: s[(s["ranking"] == r) & (s["collection"] == c)]
        plot_fig4_v2(pick("death", "Hallmark"), pick("age", "Hallmark"), pick("death", "Reactome"), lr)
        done.append("Fig4")
    else:
        LOG.warning("S19 / S21 missing - Fig 4 skipped")
    # ---------------- Fig 5 ----------------
    need5 = ["S22b_op3_signature_genes.tsv", "S23_op3_compound_results.tsv", "S24_op3_candidates.tsv",
             "S25_positive_controls.tsv"]
    if all((C.TAB / f).exists() for f in need5):
        gtab = pd.read_csv(C.TAB / need5[0], sep="\t", keep_default_na=False, na_values=[""])
        gtab["lineage"] = gtab["lineage"].fillna("")
        gtab["in_OP3"] = gtab["in_OP3"].astype(str).str.lower() == "true"
        res = pd.read_csv(C.TAB / need5[1], sep="\t")
        cand = pd.read_csv(C.TAB / need5[2], sep="\t")
        pc = pd.read_csv(C.TAB / need5[3], sep="\t")
        pc["passed"] = pc["passed"].astype(str).str.lower() == "true"
        dg = pd.read_csv(C.TAB / "S26_dgidb_targets.tsv", sep="\t") if (C.TAB / "S26_dgidb_targets.tsv").exists() \
            else pd.DataFrame(columns=["gene", "role", "drug", "approved", "matches_desired_direction"])
        if len(dg):
            for c_ in ("approved", "matches_desired_direction"):
                dg[c_] = dg[c_].astype(str).str.lower() == "true"
        primary = next(c for c in sorted(res["cell_type"].unique()) if re.search(A3.PRIMARY_CT, c, re.I))
        pc_plot = pc[pc["control"] != "(a) GR targets induced"]
        pc_a = pc[pc["control"] == "(a) GR targets induced"].sort_values("value", ascending=False) \
            .drop_duplicates("compound")
        plot_fig5_v2(gtab, res, cand, pd.concat([pc_a, pc_plot]).head(14), dg, primary)
        done.append("Fig5")
    else:
        LOG.warning("stage-4 tables missing - Fig 5 skipped")
    # ---------------- Fig 7 ----------------
    cache = C.PROC / "sc_GSE167363_processed.h5ad"
    need7 = ["S36_sc_cell_scores.tsv.gz", "S35_sc_monocyte_coupling.tsv", "S34_sc_donor_celltype_scores.tsv",
             "S34b_sc_sepsis_minus_healthy.tsv"]
    if cache.exists() and all((C.TAB / f).exists() for f in need7):
        a = ad.read_h5ad(cache, backed="r")
        um = pd.DataFrame(np.asarray(a.obsm["X_umap"]), index=a.obs_names.astype(str))
        obs = pd.read_csv(C.TAB / need7[0], sep="\t", index_col=0)
        obs.index = obs.index.astype(str)
        um = um.reindex(obs.index)
        if um.isna().any().any():
            raise SystemExit("Cells of S36 and the cached single-cell object do not match - re-run 06 first.")
        coup = pd.read_csv(C.TAB / need7[1], sep="\t")
        donor = pd.read_csv(C.TAB / need7[2], sep="\t")
        delta = pd.read_csv(C.TAB / need7[3], sep="\t", index_col=0)
        s_ = coup[(coup["group"] != "healthy") & (coup["time"] == "T0")]
        pooled = {}
        if len(s_) >= 2:
            m = dl_meta(s_["z"].to_numpy(), s_["se"].to_numpy())
            pooled = {"rho": float(np.tanh(m["est"])), "lo": float(np.tanh(m["est"] - 1.96 * m["se"])),
                      "hi": float(np.tanh(m["est"] + 1.96 * m["se"]))}
        plot_fig7_v2(obs.reset_index(drop=True), um.to_numpy(), delta,
                     coup.sort_values(["time", "group", "sample"]).reset_index(drop=True),
                     donor[donor["cell_type"] == A5.MONO], pooled)
        done.append("Fig7")
    else:
        LOG.warning("single-cell outputs missing - Fig 7 skipped (run 06 first)")
    LOG.info(f"Redrawn from the saved tables: {', '.join(done) or 'none'} ({(time.time() - t0):.0f} s). "
             f"No statistic was recomputed; previous versions are in {bk}.")


run_figures()
