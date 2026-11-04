"""Record the computing environment used to run the analysis (for the Methods section and the repository).

Writes environment_versions.txt with the operating system, processor, Python version and the versions of all
packages used by scripts 01-08.

Usage:  python tools/record_environment.py [output file]
"""
from __future__ import annotations

import os
import platform
import sys
from importlib import metadata

PACKAGES = ["numpy", "pandas", "scipy", "statsmodels", "matplotlib", "requests", "anndata", "h5py", "scanpy",
            "harmonypy", "umap-learn", "scikit-learn", "python-igraph", "igraph", "leidenalg", "openpyxl"]


def main() -> None:
    out = sys.argv[1] if len(sys.argv) > 1 else "environment_versions.txt"
    lines = [f"python\t{sys.version.split()[0]}",
             f"platform\t{platform.platform()}",
             f"machine\t{platform.machine()}",
             f"processor\t{platform.processor() or 'n/a'}",
             f"cpu_count\t{os.cpu_count()}"]
    try:
        pages, size = os.sysconf("SC_PHYS_PAGES"), os.sysconf("SC_PAGE_SIZE")
        lines.append(f"memory_GB\t{pages * size / 1e9:.1f}")
    except (ValueError, OSError, AttributeError):
        lines.append("memory_GB\tn/a")
    for p in PACKAGES:
        try:
            lines.append(f"{p}\t{metadata.version(p)}")
        except metadata.PackageNotFoundError:
            lines.append(f"{p}\tnot installed")
    open(out, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nwritten to {out}")


if __name__ == "__main__":
    main()
