# Note on the registration time stamps

All texts were written and hashed before the analyses they govern. The registration times on 2 November 2026 (UTC) were:

| Document | Registered (UTC) | Stored record (UTC) | Source of the original time |
| --- | --- | --- | --- |
| Protocol v1 (PREREG_v1) | 05:01:05 | 09:59:49 | log of stage 02 (`sepsis_analysis1.log`, first run) |
| Amendment 1 | 05:45:15 | 09:59:49 | log of stage 02 (`sepsis_analysis1.log`, first run) |
| Amendment 3 | 09:55:59 | 09:55:59 | record |
| Amendment 4 | 10:31:27 | 10:31:27 | record |
| Amendment 5 (superseded) | 11:03:00 | 11:03:00 | record |
| Amendment 6 | 11:25:36 | 11:25:36 | record |
| Amendment 7 | 16:51:45 | 16:51:45 | record |
| Amendment 8 | 17:29:21 | 17:29:21 | record |

The data-preparation stage (script 01) was re-run at 09:59 UTC, after Amendment 3 had been registered. During that run, the records of protocol v1 and Amendment 1 were regenerated with identical texts. The SHA-256 hashes were unchanged: 9a65487b... for protocol v1, the same value that the first analysis run had verified at 09:55 UTC. For this reason the stored time stamps of these two records are later than their original registration. The first log of stage 02 shows the original times:

    Pre-registration present: registered 2026-11-02T05:01:05+00:00, SHA-256 9a65487b85bbf8a9... (unchanged)
    Amendment 1 present (registered 2026-11-02T05:45:15+00:00, unchanged)
    Amendment 3 written (SHA-256 4c7f7d12d0ca67ae...) - before any outcome analysis

Log times are local time (UTC+03:30); record times are UTC.

The protocol was time-stamped and hash-locked locally; it was not deposited in a public registry before the analyses.
