# =====================================================================================
# # SEPSIS PROJECT - 08 ADDITIONAL POST HOC ANALYSES (after an internal review of the draft)
#
# HOW TO RUN: paste the whole file into ONE Jupyter cell and run it after stages 01-06.
# It reads only result tables written by the earlier stages (no new downloads):
#   S10_immune_scores.tsv (patient-level scores + clinical data), S28b_patient_scores_stage05.tsv,
#   S11_H1_ml_index_mortality.tsv, S13_scores_mortality_meta.tsv, S28a_gc_score_per_cohort.tsv,
#   S31_sweeney_comparison.tsv, S18_gene_meta.tsv.gz, S20_reversal_signature_prespecified.tsv,
#   S22b_op3_signature_genes.tsv, S36_sc_cell_scores.tsv.gz
# Optional (if present in data/raw/GSE54514/): GSE54514_series_matrix.txt.gz -> measured neutrophil proportion.
#
# # All analyses here were added after an internal review of the draft and were specified AFTER all earlier results were known.
# They are labelled post hoc in the manuscript. Output: tables S37-S44 (S45 if the GSE54514 series matrix is present) and a log.
# =====================================================================================
from __future__ import annotations

import gzip
import json
import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from scipy.special import expit

warnings.filterwarnings("ignore")
PROJECT_DIR = Path(os.environ.get("SEPSIS_PROJECT_DIR", Path.cwd() / "SEPSIS_immuneage")).resolve()
TAB = Path(os.environ.get("SEPSIS_TAB_DIR", PROJECT_DIR / "results" / "tables"))
RAW = PROJECT_DIR / "data" / "raw"
SEED = 0
N_BOOT = 1000
MIN_EVENTS_INTERACTION = 40
DISC = "GSE65682"
VAL = ["E-MTAB-4421", "E-MTAB-4451", "E-MTAB-5273", "GSE95233", "GSE54514", "GSE185263"]
COHORTS = [DISC] + VAL
NO_RISK_SCALE = ("E-MTAB-4451",)
WHOLE_BLOOD = ["GSE65682", "GSE95233", "GSE54514", "GSE185263"]
LEUKOCYTES = ["E-MTAB-4421", "E-MTAB-4451", "E-MTAB-5273"]
DAY28 = ["GSE65682", "E-MTAB-4421", "E-MTAB-4451", "E-MTAB-5273", "GSE54514"]
ERYTHROID = set("RHAG GYPB GYPE GYPA SPTA1 HBQ1 HBD HBB HBA1 HBA2 HBM AHSP CA1 KLF1 GATA1 ALAS2 EPB42 ANK1 "
                "SLC4A1 TMCC2 XK TRIM10 DMTN BPGM TSPO2 GLRX5 TFRC HMBS RBM38 MPP1 TFDP2 E2F2".split())
LOG_LINES: list[str] = []


def log(msg: str) -> None:
    print(msg)
    LOG_LINES.append(msg)


def find(name: str) -> Path:
    p = TAB / name
    if not p.exists():
        raise SystemExit(f"Missing input table {p} - run stages 01-06 first.")
    return p


# -------------------------------------------------------------------------------------
# statistics (identical to stages 02-05)
# -------------------------------------------------------------------------------------
def dl_meta(est, se) -> dict:
    est, se = np.asarray(est, float), np.asarray(se, float)
    ok = np.isfinite(est) & np.isfinite(se) & (se > 0)
    est, se = est[ok], se[ok]
    k = len(est)
    if k == 0:
        return {"k": 0, "est": np.nan, "se": np.nan, "P": np.nan, "I2": np.nan, "tau2": np.nan, "pi_lo": np.nan,
                "pi_hi": np.nan}
    w = 1 / se ** 2
    fe = np.sum(w * est) / np.sum(w)
    q = np.sum(w * (est - fe) ** 2)
    tau2 = max(0.0, (q - (k - 1)) / (np.sum(w) - np.sum(w ** 2) / np.sum(w))) if k > 1 else 0.0
    wr = 1 / (se ** 2 + tau2)
    m = np.sum(wr * est) / np.sum(wr)
    s = np.sqrt(1 / np.sum(wr))
    i2 = max(0.0, (q - (k - 1)) / q) if k > 1 and q > 0 else 0.0
    # 95% prediction interval (Riley et al., BMJ 2011): t distribution with k-2 df, requires k >= 3
    if k >= 3:
        h = stats.t.ppf(0.975, k - 2) * np.sqrt(tau2 + s ** 2)
        pi_lo, pi_hi = m - h, m + h
    else:
        pi_lo = pi_hi = np.nan
    return {"k": k, "est": m, "se": s, "P": 2 * stats.norm.sf(abs(m / s)), "I2": i2, "tau2": tau2,
            "pi_lo": pi_lo, "pi_hi": pi_hi}


def logit(y, X):
    r = sm.Logit(np.asarray(y, int), sm.add_constant(np.asarray(X, float), has_constant="add")).fit(
        disp=0, method="newton", maxiter=200)
    if not r.mle_retvals.get("converged", True):
        raise RuntimeError("no convergence")
    return r


def assoc(d: pd.DataFrame, preds: list, adjust=("age10", "male")) -> dict:
    cols = list(preds) + [a for a in adjust if a not in preds and not (a == "male" and d["male"].nunique() < 2)]
    dd = d.dropna(subset=cols + ["death"])
    try:
        r = logit(dd["death"], dd[cols])
        b, se = np.asarray(r.params)[1:len(preds) + 1], np.asarray(r.bse)[1:len(preds) + 1]
    except Exception:
        b = se = np.full(len(preds), np.nan)
    return {"n": len(dd), "events": int(dd["death"].sum()), "b": b, "se": se}


def g_formula(d: pd.DataFrame, meds: list, covs: list) -> dict:
    """Natural effects (risk difference) of +10 years of age; mediators jointly; as in stage 02, with optional
    extra covariates in the mediator and outcome models (covs always includes sex when it varies)."""
    y = d["death"].astype(int).to_numpy()
    a = d["age10"].to_numpy(float)
    Z = d[covs].to_numpy(float) if covs else np.zeros((len(d), 0))
    M = d[meds].to_numpy(float)
    Xm = np.column_stack([np.ones(len(a)), a, Z])
    shift = np.linalg.lstsq(Xm, M, rcond=None)[0][1]
    inter = y.sum() >= MIN_EVENTS_INTERACTION
    design = (lambda age_, Mx: np.column_stack([np.ones(len(age_)), age_, Mx, age_[:, None] * Mx, Z])) if inter \
        else (lambda age_, Mx: np.column_stack([np.ones(len(age_)), age_, Mx, Z]))
    r = sm.Logit(y, design(a, M)).fit(disp=0, method="newton", maxiter=200)
    if not r.mle_retvals.get("converged", True) or not np.all(np.isfinite(r.params)):
        raise RuntimeError("no convergence")
    g = r.params

    def risk(age_, Mx):
        return float(np.mean(expit(design(age_, Mx) @ g)))
    y11, y10, y00 = risk(a + 1, M + shift), risk(a + 1, M), risk(a, M)
    te, nde, nie = y11 - y00, y10 - y00, y11 - y10
    return {"TE": te, "NDE": nde, "NIE": nie}


def mediation(d: pd.DataFrame, meds: list, seed: int, extra: tuple = ()) -> dict:
    covs = (["male"] if d["male"].nunique() > 1 else []) + list(extra)
    d = d.dropna(subset=meds + covs + ["age10", "death"])
    est = g_formula(d, meds, covs)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(N_BOOT):
        idx = rng.integers(0, len(d), len(d))
        try:
            boots.append(g_formula(d.iloc[idx], meds, covs))
        except Exception:
            continue
    B = pd.DataFrame(boots)
    out = {"n": len(d), "events": int(d["death"].sum()), "boots_ok": len(B)}
    for k in ("TE", "NDE", "NIE"):
        out[k], out[f"{k}_se"] = est[k], float(B[k].std(ddof=1))
    return out


def pooled_ratio(nie: dict, te: dict, n: int = 200000, seed: int = 0):
    rng = np.random.default_rng(seed)
    a, b = rng.normal(nie["est"], nie["se"], n), rng.normal(te["est"], te["se"], n)
    r = a / b
    return nie["est"] / te["est"], *np.nanpercentile(r[np.abs(b) > 1e-6], [2.5, 97.5])


def evalue_or(or_: float, common: bool = True) -> float:
    """E-value (VanderWeele & Ding 2017) for an odds ratio; common outcome (>15%): RR ~ sqrt(OR)."""
    rr = np.sqrt(or_) if common else or_
    rr = rr if rr >= 1 else 1 / rr
    return float(rr + np.sqrt(rr * (rr - 1)))


def auc(y, p) -> float:
    y, p = np.asarray(y, int), np.asarray(p, float)
    pos, neg = p[y == 1], p[y == 0]
    return float((np.sum(pos[:, None] > neg[None, :]) + 0.5 * np.sum(pos[:, None] == neg[None, :]))
                 / (len(pos) * len(neg)))


def fisher_pool(r, n) -> dict:
    r, n = np.asarray(r, float), np.asarray(n, float)
    z, se = np.arctanh(r), 1 / np.sqrt(n - 3)
    m = dl_meta(z, se)
    return {"rho": np.tanh(m["est"]), "lo": np.tanh(m["est"] - 1.96 * m["se"]),
            "hi": np.tanh(m["est"] + 1.96 * m["se"]), "P": m["P"], "I2": m["I2"], "k": m["k"]}


def fmt_or(m: dict) -> str:
    return (f"{np.exp(m['est']):.2f} ({np.exp(m['est'] - 1.96 * m['se']):.2f}-{np.exp(m['est'] + 1.96 * m['se']):.2f}); "
            f"95% PI {np.exp(m['pi_lo']):.2f}-{np.exp(m['pi_hi']):.2f}; P = {m['P']:.2g}; I2 {m['I2'] * 100:.0f}%; k = {m['k']}")


# -------------------------------------------------------------------------------------
# data
# -------------------------------------------------------------------------------------
def load() -> dict:
    s = pd.read_csv(find("S10_immune_scores.tsv"), sep="\t", index_col=0, low_memory=False)
    s = s[s["eligible"].astype(str).str.lower().isin(["true", "1"])]
    extra = pd.read_csv(find("S28b_patient_scores_stage05.tsv"), sep="\t", index_col=0)
    s = s.join(extra[["sweeney", "gc_score"]], how="left")
    s = s.rename(columns={"sev_sofa_24h_post_admisssion": "sofa", "sev_severity_apacheii_": "apache"})
    return {c: s[s["cohort"] == c].copy() for c in COHORTS}


# -------------------------------------------------------------------------------------
# R1 - H1 robustness, components, subgroups, prediction intervals, E-value
# -------------------------------------------------------------------------------------
def r1_h1(data: dict) -> pd.DataFrame:
    log("\n== R1: H1 robustness (logistic, age + sex adjusted; DerSimonian-Laird; 95% prediction intervals)")
    per = {c: assoc(data[c], ["ml_index"]) for c in COHORTS}
    s11 = pd.read_csv(find("S11_H1_ml_index_mortality.tsv"), sep="\t").set_index("label")
    for c in VAL:                                       # reproduction check against stage 02
        assert abs(per[c]["b"][0] - s11.loc[c, "logratio"]) < 1e-6, f"H1 reproduction failed for {c}"
    log("  per-cohort logistic estimates reproduce S11 exactly")
    rows = []

    def add(label, cs, note=""):
        m = dl_meta([per[c]["b"][0] for c in cs], [per[c]["se"][0] for c in cs])
        rows.append({"analysis": label, "cohorts": ", ".join(cs), "note": note, **m,
                     "OR": np.exp(m["est"]), "lo": np.exp(m["est"] - 1.96 * m["se"]),
                     "hi": np.exp(m["est"] + 1.96 * m["se"]), "PI_lo": np.exp(m["pi_lo"]), "PI_hi": np.exp(m["pi_hi"])})
        log(f"  {label}: {fmt_or(m)}")
    add("Validation (primary H1)", VAL)
    add("All cohorts", COHORTS, "discovery logistic")
    add("Whole-blood cohorts", WHOLE_BLOOD, "discovery logistic")
    add("Isolated-leukocyte cohorts", LEUKOCYTES)
    add("28-day mortality cohorts only", DAY28, "discovery logistic")
    add("Hospital / in-hospital mortality cohorts", ["GSE95233", "GSE185263"])
    # components: each alone and mutually adjusted in one model
    for lab, preds in (("immature neutrophil (alone)", ["immature_neutrophil"]), ("T cell (alone)", ["t_cell"]),
                       ("immature neutrophil (adjusted for T cell)", ["immature_neutrophil", "t_cell"]),
                       ("T cell (adjusted for immature neutrophil)", ["t_cell", "immature_neutrophil"])):
        r_ = [assoc(data[c], preds) for c in COHORTS]
        m = dl_meta([x["b"][0] for x in r_], [x["se"][0] for x in r_])
        rows.append({"analysis": f"Component: {lab}", "cohorts": "all 7", **m, "OR": np.exp(m["est"]),
                     "lo": np.exp(m["est"] - 1.96 * m["se"]), "hi": np.exp(m["est"] + 1.96 * m["se"]),
                     "PI_lo": np.exp(m["pi_lo"]), "PI_hi": np.exp(m["pi_hi"])})
        log(f"  Component {lab}: {fmt_or(m)}")
    # correlation of the two components within cohorts
    rr = [stats.spearmanr(data[c]["immature_neutrophil"], data[c]["t_cell"], nan_policy="omit")[0] for c in COHORTS]
    fp = fisher_pool(rr, [len(data[c]) for c in COHORTS])
    log(f"  Spearman rho immature neutrophil vs T cell (pooled): {fp['rho']:.2f} ({fp['lo']:.2f} to {fp['hi']:.2f})")
    rows.append({"analysis": "Spearman rho immature neutrophil vs T cell", "cohorts": "all 7", "est": fp["rho"],
                 "lo": fp["lo"], "hi": fp["hi"], "P": fp["P"], "I2": fp["I2"], "k": fp["k"]})
    # E-value for the primary validation estimate
    v = rows[0]
    ev, ev_ci = evalue_or(v["OR"]), evalue_or(v["lo"])
    log(f"  E-value (validation OR {v['OR']:.2f}, common-outcome approximation): {ev:.2f}; for the lower CI limit {ev_ci:.2f}")
    rows.append({"analysis": "E-value for validation OR (point / lower CI)", "est": ev, "lo": ev_ci})
    # prediction intervals for the other pooled estimates already reported
    s13 = pd.read_csv(find("S13_scores_mortality_meta.tsv"), sep="\t")
    for _, r in s13.iterrows():
        k, tau2, se = int(r["k"]), float(r["tau2"]), float(r["se"])
        h = stats.t.ppf(0.975, k - 2) * np.sqrt(tau2 + se ** 2)
        rows.append({"analysis": f"Score vs death (S13): {r['score']}", "cohorts": "all 7", "k": k, "est": r["est"],
                     "OR": np.exp(r["est"]), "PI_lo": np.exp(r["est"] - h), "PI_hi": np.exp(r["est"] + h),
                     "I2": r["I2"], "P": r["P"]})
    gc = pd.read_csv(find("S28a_gc_score_per_cohort.tsv"), sep="\t")
    m = dl_meta(gc["death_logOR"], gc["death_se"])
    rows.append({"analysis": "GC-response score vs death (D1)", "cohorts": "all 7", **m, "OR": np.exp(m["est"]),
                 "PI_lo": np.exp(m["pi_lo"]), "PI_hi": np.exp(m["pi_hi"])})
    log(f"  GC score D1: {fmt_or(m)}")
    sw = pd.read_csv(find("S31_sweeney_comparison.tsv"), sep="\t")
    for col, lab in (("sweeney_logOR", "Sweeney score vs death"), ("ml_adj_sweeney_logOR", "ML index adjusted for Sweeney")):
        m = dl_meta(sw[col], sw[col.replace("logOR", "se")])
        rows.append({"analysis": lab, "cohorts": "all 7", **m, "OR": np.exp(m["est"]), "PI_lo": np.exp(m["pi_lo"]),
                     "PI_hi": np.exp(m["pi_hi"])})
        log(f"  {lab}: {fmt_or(m)}")
    out = pd.DataFrame(rows)
    out.to_csv(TAB / "S37_H1_robustness_prediction_intervals.tsv", sep="\t", index=False)
    return out


# -------------------------------------------------------------------------------------
# R2 - mediation sensitivity
# -------------------------------------------------------------------------------------
def r2_mediation(data: dict) -> pd.DataFrame:
    log("\n== R2: mediation sensitivity (g-computation as in stage 02; 1000 bootstraps per cohort)")
    specs = [("Primary: ML index, MHC-II, inflammation", ["ml_index", "mhc_class_ii", "inflammation_hallmark"], ()),
             ("Extended: + B cell + T cell", ["ml_index", "mhc_class_ii", "inflammation_hallmark", "b_cell", "t_cell"], ()),
             ("Lymphoid only: B cell + T cell", ["b_cell", "t_cell"], ())]
    rows = []
    for name, meds, extra in specs:
        per = []
        for i, c in enumerate(COHORTS):
            if c in NO_RISK_SCALE:
                continue
            r = mediation(data[c], meds, SEED + i, extra)
            per.append({"model": name, "cohort": c, **r})
        P = pd.DataFrame(per)
        v = P[P["cohort"] != DISC]
        nie, te, nde = (dl_meta(v[k], v[f"{k}_se"]) for k in ("NIE", "TE", "NDE"))
        pm = pooled_ratio(nie, te)
        d = P[P["cohort"] == DISC].iloc[0]
        rows += per
        rows.append({"model": name, "cohort": "Validation pooled", "TE": te["est"], "TE_se": te["se"],
                     "NDE": nde["est"], "NDE_se": nde["se"], "NIE": nie["est"], "NIE_se": nie["se"], "NIE_P": nie["P"],
                     "PM": pm[0], "PM_lo": pm[1], "PM_hi": pm[2], "boots_min": int(v["boots_ok"].min())})
        log(f"  {name}: validation NIE {nie['est'] * 100:.2f} pp ({(nie['est'] - 1.96 * nie['se']) * 100:.2f} to "
            f"{(nie['est'] + 1.96 * nie['se']) * 100:.2f}), P = {nie['P']:.2g}; PM {pm[0] * 100:.0f}% "
            f"({pm[1] * 100:.0f} to {pm[2] * 100:.0f}); TE {te['est'] * 100:.2f} pp; discovery NIE {d['NIE'] * 100:.2f} pp; "
            f"min bootstraps {int(v['boots_ok'].min())}")
    # severity-adjusted mediation in the one cohort with SOFA (adds SOFA to mediator and outcome models)
    d185 = data["GSE185263"].copy()
    d185["sofa_z"] = (d185["sofa"] - d185["sofa"].mean()) / d185["sofa"].std()
    for name, extra in (("GSE185263 complete cases, no SOFA", ()), ("GSE185263 complete cases, + SOFA", ("sofa_z",))):
        dd = d185.dropna(subset=["sofa_z"])
        r = mediation(dd, ["ml_index", "mhc_class_ii", "inflammation_hallmark"], SEED + 6, extra)
        rows.append({"model": name, "cohort": "GSE185263", **r})
        log(f"  {name}: n {r['n']}, TE {r['TE'] * 100:.2f} pp, NIE {r['NIE'] * 100:.2f} pp "
            f"({(r['NIE'] - 1.96 * r['NIE_se']) * 100:.2f} to {(r['NIE'] + 1.96 * r['NIE_se']) * 100:.2f})")
    out = pd.DataFrame(rows)
    out.to_csv(TAB / "S38_mediation_sensitivity.tsv", sep="\t", index=False)
    # correlations among mediators (within cohort, pooled)
    meds = ["ml_index", "mhc_class_ii", "inflammation_hallmark", "b_cell", "t_cell"]
    cor = []
    for i, a in enumerate(meds):
        for b in meds[i + 1:]:
            rr = [stats.spearmanr(data[c][a], data[c][b], nan_policy="omit")[0] for c in COHORTS]
            fp = fisher_pool(rr, [len(data[c]) for c in COHORTS])
            cor.append({"a": a, "b": b, **fp})
    cor = pd.DataFrame(cor)
    cor.to_csv(TAB / "S39_mediator_correlations.tsv", sep="\t", index=False)
    log("  mediator correlations (pooled Spearman):\n" + cor[["a", "b", "rho", "lo", "hi", "I2"]].round(2).to_string(index=False))
    return out


# -------------------------------------------------------------------------------------
# R3 - external calibration and decision curves
# -------------------------------------------------------------------------------------
def r3_calibration(data: dict) -> pd.DataFrame:
    log("\n== R3: external calibration and net benefit (models fitted in GSE65682, applied unchanged)")
    disc = data[DISC].dropna(subset=["age10", "male", "ml_index", "sweeney", "death"])
    models = {"base": ["age10", "male"], "base+ML": ["age10", "male", "ml_index"],
              "base+Sweeney": ["age10", "male", "sweeney"], "base+Sweeney+ML": ["age10", "male", "sweeney", "ml_index"]}
    fits = {k: logit(disc["death"], disc[v]) for k, v in models.items()}
    rows, pooled = [], []
    for c in VAL:
        d = data[c].dropna(subset=["age10", "male", "ml_index", "sweeney", "death"])
        y = d["death"].astype(int).to_numpy()
        for k, v in models.items():
            lp = np.asarray(sm.add_constant(d[v].astype(float), has_constant="add") @ fits[k].params)
            p = expit(lp)
            slope = sm.Logit(y, sm.add_constant(lp)).fit(disp=0)
            citl = sm.GLM(y, np.ones((len(y), 1)), family=sm.families.Binomial(), offset=lp).fit()
            rows.append({"cohort": c, "model": k, "n": len(y), "events": int(y.sum()), "AUC": auc(y, p),
                         "cal_intercept": float(citl.params[0]), "cal_intercept_se": float(citl.bse[0]),
                         "cal_slope": float(slope.params[1]), "cal_slope_se": float(slope.bse[1]),
                         "brier": float(np.mean((p - y) ** 2)), "observed_risk": float(y.mean()),
                         "mean_predicted_risk": float(p.mean())})
            if c not in NO_RISK_SCALE:
                pooled.append(pd.DataFrame({"model": k, "y": y, "p": p, "cohort": c}))
    out = pd.DataFrame(rows)
    out.to_csv(TAB / "S40_external_calibration.tsv", sep="\t", index=False)
    s = out.groupby("model")[["AUC", "cal_slope", "brier"]].mean()
    log("  mean over the six validation cohorts:\n" + s.round(3).to_string())
    for k in models:
        mi = dl_meta(out.loc[(out.model == k) & (~out.cohort.isin(NO_RISK_SCALE)), "cal_intercept"],
                     out.loc[(out.model == k) & (~out.cohort.isin(NO_RISK_SCALE)), "cal_intercept_se"])
        ms = dl_meta(out.loc[out.model == k, "cal_slope"], out.loc[out.model == k, "cal_slope_se"])
        log(f"  {k}: pooled calibration intercept {mi['est']:.2f} ({mi['est'] - 1.96 * mi['se']:.2f} to "
            f"{mi['est'] + 1.96 * mi['se']:.2f}; E-MTAB-4451 excluded); pooled slope {ms['est']:.2f} "
            f"({ms['est'] - 1.96 * ms['se']:.2f} to {ms['est'] + 1.96 * ms['se']:.2f})")
    # decision curve (patients of the five risk-scale validation cohorts pooled)
    P = pd.concat(pooled)
    dca = []
    for t in (0.10, 0.15, 0.20, 0.25, 0.30):
        prev = P.loc[P.model == "base", "y"].mean()
        dca.append({"threshold": t, "model": "treat all", "net_benefit": prev - (1 - prev) * t / (1 - t)})
        for k in models:
            q = P[P.model == k]
            tp = np.sum((q.p >= t) & (q.y == 1)) / len(q)
            fp = np.sum((q.p >= t) & (q.y == 0)) / len(q)
            dca.append({"threshold": t, "model": k, "net_benefit": tp - fp * t / (1 - t)})
    dca = pd.DataFrame(dca)
    dca.to_csv(TAB / "S41_decision_curve.tsv", sep="\t", index=False)
    log("  net benefit (five validation cohorts pooled, n = %d):\n" % (len(P) // len(models))
        + dca.pivot(index="threshold", columns="model", values="net_benefit").round(3).to_string())
    return out


# -------------------------------------------------------------------------------------
# R4 - severity: correlation with the scores, joint adjustment, timing information
# -------------------------------------------------------------------------------------
def r4_severity(data: dict) -> pd.DataFrame:
    log("\n== R4: severity")
    rows = []
    for c, sev in (("GSE185263", "sofa"), ("GSE54514", "apache")):
        d = data[c].dropna(subset=[sev]).copy()
        d["sev_z"] = (d[sev] - d[sev].mean()) / d[sev].std()
        for s in ("ml_index", "mhc_class_ii", "gc_score", "immature_neutrophil", "t_cell"):
            rho, p = stats.spearmanr(d[s], d[sev])
            rows.append({"cohort": c, "severity": sev, "score": s, "n": len(d), "spearman_rho": rho, "P": p})
            log(f"  {c}: Spearman rho({s}, {sev}) = {rho:.2f} (P = {p:.2g}, n = {len(d)})")
        for lab, preds, adj in (("ML index, + severity", ["ml_index"], ("age10", "male", "sev_z")),
                                ("ML index, + severity + MHC-II", ["ml_index"], ("age10", "male", "sev_z", "mhc_class_ii")),
                                ("MHC-II, + severity", ["mhc_class_ii"], ("age10", "male", "sev_z")),
                                ("severity, + age + sex", ["sev_z"], ("age10", "male"))):
            r = assoc(d, preds, adj)
            rows.append({"cohort": c, "severity": sev, "score": lab, "n": r["n"], "OR": np.exp(r["b"][0]),
                         "lo": np.exp(r["b"][0] - 1.96 * r["se"][0]), "hi": np.exp(r["b"][0] + 1.96 * r["se"][0])})
            log(f"  {c}: {lab}: OR {np.exp(r['b'][0]):.2f} ({np.exp(r['b'][0] - 1.96 * r['se'][0]):.2f}-"
                f"{np.exp(r['b'][0] + 1.96 * r['se'][0]):.2f})")
    out = pd.DataFrame(rows)
    out.to_csv(TAB / "S42_severity_sensitivity.tsv", sep="\t", index=False)
    return out


# -------------------------------------------------------------------------------------
# R5 - erythroid content of the gene-level results and of the death-age signature
# -------------------------------------------------------------------------------------
def r5_erythroid() -> pd.DataFrame:
    log("\n== R5: erythroid genes")
    g = pd.read_csv(find("S18_gene_meta.tsv.gz"), sep="\t")
    sig = pd.read_csv(find("S20_reversal_signature_prespecified.tsv"), sep="\t")
    op3 = pd.read_csv(find("S22b_op3_signature_genes.tsv"), sep="\t")
    ery = g["gene"].isin(ERYTHROID)
    fdr = g["death_FDR"] < 0.05
    rows = [{"item": "erythroid marker genes measured in all cohorts", "value": int(ery.sum())},
            {"item": "death FDR < 0.05, all genes (higher / lower)",
             "value": f"{int(fdr.sum())} ({int((fdr & (g.death_z > 0)).sum())} / {int((fdr & (g.death_z < 0)).sum())})"},
            {"item": "death FDR < 0.05, erythroid genes", "value": int((fdr & ery).sum())},
            {"item": "death FDR < 0.05 after excluding erythroid genes (higher / lower)",
             "value": f"{int((fdr & ~ery).sum())} ({int((fdr & ~ery & (g.death_z > 0)).sum())} / "
                      f"{int((fdr & ~ery & (g.death_z < 0)).sum())})"},
            {"item": "signature UP genes / erythroid among them",
             "value": f"{int((sig.direction == 'UP').sum())} / {int(((sig.direction == 'UP') & sig.gene.isin(ERYTHROID)).sum())}"},
            {"item": "signature after excluding erythroid genes (UP / DOWN)",
             "value": f"{int(((sig.direction == 'UP') & ~sig.gene.isin(ERYTHROID)).sum())} / "
                      f"{int(((sig.direction == 'DOWN') & ~sig.gene.isin(ERYTHROID)).sum())}"},
            {"item": "erythroid UP genes measured in OP3 / absent",
             "value": f"{int((op3.gene.isin(ERYTHROID) & op3.in_OP3).sum())} / {int((op3.gene.isin(ERYTHROID) & ~op3.in_OP3).sum())}"}]
    for name, gs in (("ML index genes", ["CEACAM8", "MMP8", "LCN2", "LTF", "CAMP", "OLFM4", "DEFA4", "BPI", "CRISP3", "CD24",
                                         "CD3E", "CD3D", "CD3G", "CD247", "CD2", "LCK", "CD28", "IL7R", "CCR7", "TRAT1"]),
                     ("MHC class II genes", ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DMA", "HLA-DMB", "CD74"])):
        rows.append({"item": f"erythroid genes in {name}", "value": len(set(gs) & ERYTHROID)})
    out = pd.DataFrame(rows)
    for _, r in out.iterrows():
        log(f"  {r['item']}: {r['value']}")
    out.to_csv(TAB / "S43_erythroid_sensitivity.tsv", sep="\t", index=False)
    return out


# -------------------------------------------------------------------------------------
# R6 - single cell: patient-level inference
# -------------------------------------------------------------------------------------
def r6_singlecell() -> pd.DataFrame:
    log("\n== R6: single-cell coupling, patient-level inference")
    c = pd.read_csv(find("S36_sc_cell_scores.tsv.gz"), sep="\t")
    m = c[c["cell_type"] == "CD14 monocytes"].copy()
    m["logc"] = np.log(m["total_counts"])
    rows = []
    for s, d in m.groupby("sample"):
        if len(d) < 100:
            rows.append({"sample": s, "n": len(d), "rho": np.nan, "note": "< 100 monocytes"})
            continue
        r = d[["gc_score", "mhc_score", "logc"]].rank()
        res = lambda a: r[a] - sm.OLS(r[a], sm.add_constant(r["logc"])).fit().fittedvalues
        rho = float(np.corrcoef(res("gc_score"), res("mhc_score"))[0, 1])
        # cell bootstrap (assumes independent cells within a sample; reported for transparency)
        rng = np.random.default_rng(SEED)
        bs = []
        for _ in range(500):
            dd = d.iloc[rng.integers(0, len(d), len(d))]
            rr = dd[["gc_score", "mhc_score", "logc"]].rank()
            e1 = rr["gc_score"] - sm.OLS(rr["gc_score"], sm.add_constant(rr["logc"])).fit().fittedvalues
            e2 = rr["mhc_score"] - sm.OLS(rr["mhc_score"], sm.add_constant(rr["logc"])).fit().fittedvalues
            bs.append(np.corrcoef(e1, e2)[0, 1])
        rows.append({"sample": s, "group": d["group"].iloc[0], "time": d["time"].iloc[0], "n": len(d), "rho": rho,
                     "boot_lo": np.percentile(bs, 2.5), "boot_hi": np.percentile(bs, 97.5)})
    out = pd.DataFrame(rows)
    t0 = out[(out["group"] != "healthy") & (out["time"] == "T0") & out["rho"].notna()]
    z = np.arctanh(t0["rho"])
    tt = stats.ttest_1samp(z, 0.0)
    log(out.round(3).to_string(index=False))
    log(f"  septic T0 patients (k = {len(t0)}): mean Fisher z {z.mean():.3f}; one-sample t test on patient-level "
        f"z values (patients as units) P = {tt.pvalue:.2g}")
    out.to_csv(TAB / "S44_sc_patient_level_coupling.tsv", sep="\t", index=False)
    return out


# -------------------------------------------------------------------------------------
# R7 - optional: ML index vs measured neutrophil proportion (GSE54514)
# -------------------------------------------------------------------------------------
def r7_cbc(data: dict) -> None:
    p = RAW / "GSE54514" / "GSE54514_series_matrix.txt.gz"
    if not p.exists():
        log("\n== R7: GSE54514 series matrix not found - measured neutrophil proportion not analysed")
        return
    with gzip.open(p, "rt", errors="ignore") as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh if l.startswith("!Sample_")]
    meta = {}
    for l in lines:
        meta.setdefault(l[0], []).append([x.strip('"') for x in l[1:]])
    gsm = meta["!Sample_geo_accession"][0]
    neut = None
    for row in meta.get("!Sample_characteristics_ch1", []):
        if any("neutrophil proportion" in x for x in row):
            neut = pd.Series([float(x.split(":")[1]) if ":" in x and x.split(":")[1].strip() not in ("", "NA") else np.nan
                              for x in row], index=gsm)
    if neut is None:
        log("\n== R7: neutrophil proportion field not found")
        return
    d = data["GSE54514"].join(neut.rename("neut_prop"), how="left")
    log(f"\n== R7: GSE54514 measured neutrophil proportion available for {int(d['neut_prop'].notna().sum())} patients")
    rows = []
    for s in ("ml_index", "immature_neutrophil", "mature_neutrophil", "t_cell", "mhc_class_ii"):
        dd = d.dropna(subset=["neut_prop", s])
        rho, pv = stats.spearmanr(dd[s], dd["neut_prop"])
        rows.append({"score": s, "n": len(dd), "rho": rho, "P": pv})
        log(f"  Spearman rho({s}, measured neutrophil proportion) = {rho:.2f} (P = {pv:.2g})")
    pd.DataFrame(rows).to_csv(TAB / "S45_cbc_validation.tsv", sep="\t", index=False)


def run() -> None:
    t0 = time.time()
    data = load()
    for c in COHORTS:
        log(f"{c}: {len(data[c])} patients, {int(data[c]['death'].sum())} deaths")
    r1_h1(data)
    r2_mediation(data)
    r3_calibration(data)
    r4_severity(data)
    r5_erythroid()
    r6_singlecell()
    r7_cbc(data)
    log(f"\nRun time {(time.time() - t0) / 60:.1f} min; python {sys.version.split()[0]}, pandas {pd.__version__}, "
        f"numpy {np.__version__}, statsmodels {sm.__version__}")
    (PROJECT_DIR / "logs").mkdir(parents=True, exist_ok=True)
    (PROJECT_DIR / "logs" / "sepsis_revision.log").write_text("\n".join(LOG_LINES))


run()
