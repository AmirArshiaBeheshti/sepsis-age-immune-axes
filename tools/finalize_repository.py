# =====================================================================================
# FINALIZE THE REPOSITORY - run once, after scripts 01-08 have finished
#
# HOW TO RUN: paste the whole file into ONE Jupyter cell (or: python tools/finalize_repository.py) in the folder
#   that contains BOTH the project folder (SEPSIS_immuneage) and the unzipped repository (sepsis-age-immune-axes).
#   Change the two paths below if your folders are elsewhere. Only the Python standard library is needed.
#
# What it does:
#   1  copies the registered protocol and amendments (.md + .json) from SEPSIS_immuneage/preregistration/
#      into sepsis-age-immune-axes/protocol/
#   2  recomputes the SHA-256 hash of every .md file exactly as the analysis scripts do and compares it with the
#      registered hash; writes protocol/registration_index.tsv (document, UTC time stamp, hash, status)
#   3  records the software and hardware versions -> environment_versions.txt (in the repository)
#   4  reads the dates of the downloaded public data (file dates in SEPSIS_immuneage/data/raw/)
#   5  prints a ready-to-paste text for Methods section 2.14 and writes it to methods_2.14_software.txt
#   Nothing in SEPSIS_immuneage is changed.
# =====================================================================================
from __future__ import annotations

import datetime
import hashlib
import json
import os
import platform
import shutil
import sys
from importlib import metadata
from pathlib import Path

PROJECT_DIR = Path(os.environ.get("SEPSIS_PROJECT_DIR", Path.cwd() / "SEPSIS_immuneage")).resolve()
REPO_DIR = Path(os.environ.get("SEPSIS_REPO_DIR", Path.cwd() / "sepsis-age-immune-axes")).resolve()

PACKAGES = ["numpy", "pandas", "scipy", "statsmodels", "matplotlib", "requests", "anndata", "h5py", "scanpy",
            "harmonypy", "umap-learn", "scikit-learn", "python-igraph", "igraph", "leidenalg", "openpyxl"]
# names used in the Methods text (package -> label)
METHODS_PACKAGES = [("numpy", "NumPy"), ("pandas", "pandas"), ("scipy", "SciPy"), ("statsmodels", "statsmodels"),
                    ("scikit-learn", "scikit-learn"), ("matplotlib", "matplotlib"), ("scanpy", "Scanpy"),
                    ("anndata", "anndata"), ("harmonypy", "harmonypy"), ("umap-learn", "umap-learn"),
                    ("leidenalg", "leidenalg"), ("requests", "requests")]
EXPECTED_DOCS = ["PREREG_v1", "AMENDMENT_1", "AMENDMENT_3", "AMENDMENT_4", "AMENDMENT_5", "AMENDMENT_6",
                 "AMENDMENT_7", "AMENDMENT_8"]


def version(pkg: str) -> str | None:
    try:
        return metadata.version(pkg)
    except metadata.PackageNotFoundError:
        return None


# -------------------------------------------------------------------------------------
# 1-2  protocol: copy and verify
# -------------------------------------------------------------------------------------
def copy_and_verify_protocol() -> bool:
    src = PROJECT_DIR / "preregistration"
    dst = REPO_DIR / "protocol"
    if not src.is_dir():
        raise SystemExit(f"Folder not found: {src}\nSet PROJECT_DIR at the top of this cell.")
    if not REPO_DIR.is_dir():
        raise SystemExit(f"Repository folder not found: {REPO_DIR}\nUnzip sepsis-age-immune-axes.zip here or set REPO_DIR.")
    dst.mkdir(parents=True, exist_ok=True)
    files = sorted([p for p in src.iterdir() if p.suffix in (".md", ".json")])
    for p in files:
        shutil.copy2(p, dst / p.name)                     # copy2 keeps the original file dates
    print(f"1  copied {len(files)} files from {src}\n   to {dst}")

    rows, ok_all = [], True
    print("\n2  verification (SHA-256 of each .md file vs its registered .json record)")
    for rec in sorted(dst.glob("*.json")):
        meta = json.loads(rec.read_text())
        doc = dst / meta.get("file", rec.with_suffix(".md").name)
        if not doc.exists():
            status, digest = "MISSING DOCUMENT", ""
        else:
            digest = hashlib.sha256(doc.read_text().encode()).hexdigest()   # identical to the analysis scripts
            status = "OK" if digest == meta.get("sha256") else "HASH MISMATCH"
        ok_all &= status == "OK"
        rows.append((doc.name, str(meta.get("registered_utc", "")), str(meta.get("sha256", "")), status))
        print(f"   {status:16s} {doc.name:16s} registered {meta.get('registered_utc', '?')}  "
              f"sha256 {str(meta.get('sha256', ''))[:16]}...")
    found = {Path(r[0]).stem for r in rows}
    missing = [d for d in EXPECTED_DOCS if d not in found]
    if missing:
        print(f"   NOTE: no record for {', '.join(missing)} (expected for documents that were never written)")
    (dst / "registration_index.tsv").write_text(
        "document\tregistered_utc\tsha256\tstatus\n" + "\n".join("\t".join(r) for r in rows) + "\n")
    print(f"   {sum(r[3] == 'OK' for r in rows)} of {len(rows)} documents verified -> protocol/registration_index.tsv")
    if not ok_all:
        print("   WARNING: at least one document does not match its record - do not publish before checking it.")
    return ok_all


# -------------------------------------------------------------------------------------
# 3  environment
# -------------------------------------------------------------------------------------
def record_environment() -> dict:
    env = {"python": sys.version.split()[0], "platform": platform.platform(), "machine": platform.machine(),
           "processor": platform.processor() or "n/a", "cpu_count": str(os.cpu_count())}
    try:
        env["memory_GB"] = f"{os.sysconf('SC_PHYS_PAGES') * os.sysconf('SC_PAGE_SIZE') / 1e9:.1f}"
    except (ValueError, OSError, AttributeError):
        env["memory_GB"] = "n/a"
    for p in PACKAGES:
        env[p] = version(p) or "not installed"
    out = REPO_DIR / "environment_versions.txt"
    out.write_text("\n".join(f"{k}\t{v}" for k, v in env.items()) + "\n")
    print(f"\n3  environment recorded -> {out.name}")
    for k, v in env.items():
        print(f"   {k:14s} {v}")
    return env


# -------------------------------------------------------------------------------------
# 4  data retrieval dates
# -------------------------------------------------------------------------------------
def retrieval_dates() -> tuple[str, str] | None:
    raw = PROJECT_DIR / "data" / "raw"
    if not raw.is_dir():
        print(f"\n4  {raw} not found - retrieval dates not available")
        return None
    per_source, all_dates = {}, []
    for top in sorted(p for p in raw.iterdir() if p.is_dir()):
        ds = [datetime.date.fromtimestamp(f.stat().st_mtime) for f in top.rglob("*")
              if f.is_file() and not f.name.endswith(".part")]
        if ds:
            per_source[top.name] = (min(ds), max(ds))
            all_dates += ds
    if not all_dates:
        print("\n4  no downloaded files found")
        return None
    print("\n4  dates of the downloaded files (data/raw), first - last:")
    for k, (a, b) in per_source.items():
        print(f"   {k:14s} {a} - {b}")
    first, last = min(all_dates), max(all_dates)
    print(f"   overall: {first} - {last}")
    return first.strftime("%-d %B %Y") if os.name != "nt" else first.strftime("%#d %B %Y"), \
        last.strftime("%-d %B %Y") if os.name != "nt" else last.strftime("%#d %B %Y")


# -------------------------------------------------------------------------------------
# 5  Methods text
# -------------------------------------------------------------------------------------
def methods_text(env: dict, dates) -> str:
    pk = ", ".join(f"{label} {env[p]}" for p, label in METHODS_PACKAGES if env.get(p) not in (None, "not installed"))
    ig = env.get("python-igraph") if env.get("python-igraph") != "not installed" else env.get("igraph")
    if ig and ig != "not installed":
        pk += f", python-igraph {ig}"
    text = (f"All analyses were run in Python {env['python']} ({pk}) on {env['platform']} "
            f"({env['cpu_count']} logical CPUs, {env['memory_GB']} GB memory).")
    if dates:
        text += (f" Public data were retrieved on {dates[0]}." if dates[0] == dates[1]
                 else f" Public data were retrieved between {dates[0]} and {dates[1]}.")
    text += (" The random seed was fixed at 0. Eight scripts (01–08), each executable as a single notebook cell, "
             "reproduce every result from the public sources. Each script writes quality gates, tables and figures "
             "and refuses to run if a registered protocol text has changed.")
    (REPO_DIR / "methods_2.14_software.txt").write_text(text + "\n")
    print("\n5  text for Methods 2.14 (also saved to methods_2.14_software.txt):\n")
    print(text)
    return text


def main() -> None:
    print(f"project folder:    {PROJECT_DIR}\nrepository folder: {REPO_DIR}\n")
    ok = copy_and_verify_protocol()
    env = record_environment()
    dates = retrieval_dates()
    methods_text(env, dates)
    print("\nDONE - " + ("all protocol records verified." if ok else "check the protocol warnings above."))
    print("Note: run this cell in the same Python environment that was used for scripts 01-08, "
          "otherwise the recorded versions describe a different environment.")


main()
