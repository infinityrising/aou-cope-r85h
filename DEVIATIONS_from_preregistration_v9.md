# Deviations from the preregistration — v9 analysis (transparency ledger)

*2026-07-06. Companion to `PREREGISTRATION_AMENDMENT_1.md`. Every departure from `preregistration_v5` (and the v4
design) is listed here — documented deviations are legitimate; hidden ones are not. Governing plan = v5 + Amendment #1.*

| # | Pre-reg specified | What we did | Material? | Handling |
|---|---|---|---|---|
| D1 | **Engine = SAIGE mixed-model** (GRM for relatedness + SPA) | OLS / logistic + ancestry PCs | Yes (confirmatory) | Disclosed in Amendment #1 §4. SAIGE = planned sensitivity, not run. |
| D2 | **Data partition**: v8-overlap=discovery, v9-net-new held-out=confirmatory | Used the whole RNA∩PAV set for discovery; no held-out | Yes | Amendment #1: RNA arm (46 carriers) too small to partition → declared discovery in full. |
| D3 | **H_primary = AQ-strict PROTECTIVE** | Discovery found AQ-strict **RISK** (opposite) | Yes | Amendment #1: reclassified discovery; reversed direction; locked a fresh confirmatory (phenome). |
| D4 | Sentinel (1348881) EXCLUDED from association tests | Not excluded | **Re-evaluated → NOT a required exclusion** | The pre-reg exclusion assumed 1348881 was the convergence-target proband AND a sentinel-DERIVED phenotype. Both changed: **OA37 (external) is the real proband (1348881 ≠ OA37)**, and the v9 phenotype is **proband-BLIND ICD**. So 1348881 is a **valid population participant → KEEP in primary.** Residual = she was pre-identified as R85H+AQ+COPA, so *if* she is one of the 4 R85H+AQ+ cases, **disclose it + report a with/without sensitivity** (`sentinel_check_v9.py`). Do NOT reflexively exclude. |
| D5 | **HLA/MHC as covariate/differential** (RA→DRB1, ANCA→HLA-DP/DQ; net-zero confirmatory) | **Not examined at all** | Low (cross-chromosome interaction) | **OPEN** — optional add as a phenome covariate; won't change the underpowered verdict. Document. |
| D6 | **Layered L0→L5** with per-layer GO/NO-GO; L0 = allelic landscape + per-ancestry HWE/F_IS + founder-haplotype cis-enumeration + off-founder EUR/AMR/MID carrier harvest | Ran an abbreviated **molecular-first** arm (STING×ISG) + phenome; skipped the full L0 landscape, the de-confounding GATE, EUR/AMR harvest | Partial | Discovery scope was narrowed to the STING/IFN spine. Document; the skipped layers remain available. |
| D7 | **Specified covariates/controls**: APOL1 (renal), IFIH1 rs1990760 interaction, autoimmune-PRS *negative control*, MUC5B/telomere (ILD), DADA2 differential | None of these included | Low–moderate | Not run. The founder-tilt was instead controlled by 16 PCs + the R85H− penetrance control (n=399) + cell-composition. Document. |
| D8 | **SAIGE-GENE+ COPI/innate burden** (R85H-excluded) | Ran a simple BQ/OLS damaging-COPI burden × HAQ (null) + IRAK panel | Partial | Exploratory burden only; not the registered gene-based test. Document. |
| D9 | **Self-validating ladder SAVI↑ / AQ↓** as score validation | AQ end null; **SAVI end untestable (0 carriers in RNA)** | Yes | Disclosed: STING-specific validation of the ISG score is **incomplete** (score is an externally-validated IFN signature, but STING-directionality not shown here). |

## Net assessment

The v9 work is **honest DISCOVERY that diverged from the detailed v5 confirmatory SAP.** That is acceptable *because*
we reclassified it as discovery (Amendment #1) and locked the one confirmatory (phenome) before its result — but the
divergence is broader than Amendment #1 alone stated, hence this ledger. **Two items warrant action before any
write-up: D4 (exclude the sentinel — circularity) and, optionally, D5 (HLA covariate).** D1–D2, D6–D8 are scope/engine
choices appropriate to a discovery phase and are simply disclosed. Nothing here rescues the AQ signal — it remains
underpowered/exploratory regardless.
