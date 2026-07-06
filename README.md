# AoU COPE-R85H / IRAK3 — analysis code

Population-genetics analysis of **COPE p.R85H** (rs34510432) and a candidate second hit
(**IRAK3** loss-of-function, and the wider innate-immune pathway) in the NIH *All of Us* Researcher
Workbench (v9, CDR `C2025Q4R6`), for the COPA-syndrome-spectrum study. Plan → `PREREGISTRATION.md`.

## Pipeline (phases)
| Phase | Script | Environment | Status |
|-------|--------|-------------|--------|
| **0–1 · Recon + descriptive** | `aou_driver_v9.py` | standard JupyterLab (BigQuery + plink2) | ✅ runs |
| **1 · Hail genomics** | `hail_genomics_v9.py` | **Hail Genomic Analysis (Dataproc/Spark)** | scaffold — validate on first run |
| **1 · STING phasing** | `sting_haplotypes_v9.py` | bcftools/tabix on long-read jointcall | to add |
| **2 · Confirmatory (SAIGE)** | *added at lock* | Dataproc | **firewall-locked** — built + git-timestamped at SHA-lock, run once |

## Run order in the AoU workspace
1. Clone / pull this repo, then update before each run:
   ```
   cd aou-cope-r85h && git pull
   ```
2. **Standard JupyterLab** (BigQuery + plink2):
   ```
   %run aou-cope-r85h/aou_driver_v9.py
   ```
   → environment check, carrier counts, analysis-ready cohort (QC + relatedness + ancestry),
   per-ancestry R85H, RNA-seq gating, R85H het/hom zygosity. Prints one compact SUMMARY.
3. **Hail Genomic Analysis app** (for VDS/MatrixTable work):
   - Environments → delete current → create → **App = "JupyterLab Spark cluster"** (this is the Hail/Dataproc env;
     the standard "JupyterLab" app has no pyspark).
   - Preflight in a terminal: `python -c "import hail, pyspark; print(hail.__version__, pyspark.__version__)"`
   - ```
     %run aou-cope-r85h/hail_genomics_v9.py
     ```
     → R85H zygosity/HWE per ancestry, STING diplotypes, innate-pathway burden mask.
     *Hail reads `gs://` paths (Dataproc workers can't see the FUSE mount).*

## Key resolved facts (recon, 2026-07-06)
- CDR (BigQuery): `wb-silky-artichoke-2408.C2025Q4R6` (Verily sets no `WORKSPACE_CDR`).
- Genomic files: `gs://vwb-aou-datasets-controlled/v9/...` (FUSE-mounted at `$HOME/workspace/vwb-aou-datasets-controlled-v9/`).
- `cb_variant_to_person.person_ids` is `STRUCT<list ARRAY<STRUCT<element INT64>>>` → `UNNEST(person_ids.list).element`.
- R85H = 3,135 carriers (82.5% AFR); IRAK3 L210* = 6 (ultra-rare); AFR R85H∩RNA-seq = 39.

## Discipline
- **Recon / descriptive** code runs freely (burns no confirmatory standing).
- The **confirmatory** pipeline is added and git-timestamped at SHA-lock, then run **once** on the confirmatory partition — no peeking.
- **No participant data leaves the workspace.** This repository is code only.
