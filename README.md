# AoU COPE-R85H / IRAK3 — analysis code

Population-genetics analysis of **COPE p.R85H** (rs34510432) and a candidate second hit
(**IRAK3** loss-of-function) in the NIH *All of Us* Researcher Workbench (v9, CDR `C2025Q4R6`),
for the COPA-syndrome-spectrum study.

## Use in the AoU workspace
1. Clone this repo into the workspace (Verily Workbench → **Clone a repo**), or it auto-clones on app start.
2. In a JupyterLab cell: `%run aou_driver_v9.py`
3. Paste the printed **SUMMARY** block back for interpretation.

## Contents
- `aou_driver_v9.py` — consolidated **recon / descriptive** driver: environment check, carrier counts
  (COPE R85H, STING, IRAK3), analysis-ready cohort (QC + relatedness + ancestry), per-ancestry R85H,
  and the RNA-seq gating number. Bounded, fail-safe; prints one compact SUMMARY.

## Discipline (pre-registration)
- This holds **recon / descriptive** code only until the pre-registration is SHA-locked.
- The **confirmatory** pipeline is added and git-timestamped at lock, then run **once** on the
  confirmatory partition — no peeking.
- **No participant data leaves the workspace.** This repository is code only.
