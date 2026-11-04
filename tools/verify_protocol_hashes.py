"""Verify the registered protocol and amendments.

For every <name>.json record written by the analysis scripts, recompute the SHA-256 hash of <name>.md exactly as the
scripts do (UTF-8 text) and compare it with the registered hash. Prints one line per document and writes
registration_index.tsv next to the files.

Usage:  python tools/verify_protocol_hashes.py [folder]
        (default folder: protocol/ ; the scripts write the files to SEPSIS_immuneage/preregistration/)
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


def main() -> int:
    folder = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[1] / "protocol"
    records = sorted(folder.glob("*.json"))
    if not records:
        print(f"No .json registration records found in {folder}")
        return 1
    rows, failed = [], 0
    for rec in records:
        meta = json.loads(rec.read_text())
        doc = folder / meta.get("file", rec.with_suffix(".md").name)
        if not doc.exists():
            status, digest = "MISSING DOCUMENT", ""
        else:
            digest = hashlib.sha256(doc.read_text().encode()).hexdigest()
            status = "OK" if digest == meta.get("sha256") else "HASH MISMATCH"
        failed += status != "OK"
        rows.append((doc.name, meta.get("registered_utc", ""), meta.get("sha256", ""), status))
        print(f"{status:16s} {doc.name:22s} registered {meta.get('registered_utc', '?')}  sha256 {meta.get('sha256', '')[:16]}...")
    out = folder / "registration_index.tsv"
    out.write_text("document\tregistered_utc\tsha256\tstatus\n" + "\n".join("\t".join(r) for r in rows) + "\n")
    print(f"\n{len(rows) - failed} of {len(rows)} documents verified; index written to {out}")
    return 0 if failed == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
