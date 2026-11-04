# Registered protocol and amendments

This folder contains the registered texts and their records:

- `PREREG_v1.md`, `AMENDMENT_1.md`, `AMENDMENT_3.md` to `AMENDMENT_8.md`: the registered texts
- `PREREG_v1.json`, `AMENDMENT_*.json`: the records written by the analysis scripts, holding the UTC time stamp
  (`registered_utc`) and the SHA-256 hash (`sha256`) of each text; `AMENDMENT_6.json` also holds the audit
  showing that no drug result existed when Amendment 6 was written
- `registration_index.tsv`: one line per document (time stamp, full hash, verification status)
- `REGISTRATION_NOTE.md`: original registration times of protocol v1 and Amendment 1

Verify the texts against their records with:

    python tools/verify_protocol_hashes.py protocol

Do not edit the `.md` files: the analysis scripts stop if a registered text differs from its record.
The protocol was time-stamped and hash-locked locally; it was not deposited in a public registry before the analyses.
