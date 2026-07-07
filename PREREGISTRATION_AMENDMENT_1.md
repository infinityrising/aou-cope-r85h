# Preregistration Amendment #1 — v9 STING × COPE-R85H conditional penetrance

**Date:** 2026-07-06 · **Amends:** `preregistration_v5` (Shum Lab, internal lock 2026-05-27) ·
**Locked by:** git commit + GitHub push (server-side push timestamp) **BEFORE the confirmatory phenome analysis was run
or observed.** This is an INTERNAL statistical-analysis-plan lock, NOT a formal public preregistration — mirror to
OSF/AsPredicted before manuscript for a citable time-of-record (do not claim a formal PS4/time-stamped registration).

---

## 1. Reason for amendment (full transparency)

The v9 long-read + RNA-seq analysis (2026-07-06) **inverted the registered primary direction.** v5 registered
**H_primary: AQ-strict STING1 is PROTECTIVE** (dampens interferon / COPA-spectrum among COPI/R85H carriers, β_INT<0).
The discovery analysis found the **opposite**: AQ-strict is **risk-permissive** — it *raises* the 6-gene ISG in R85H
carriers (exploratory rank p=0.014; joint pooled-adjusted β +0.80, p=0.008; interferon-activation penetrance 36% [4/11]
vs ~14% baseline). A registered directional hypothesis **cannot be "confirmed" by the opposite result**, and these
discovery p-values are exploratory (the model was iterated after seeing data — added AQ after HAQ, added covariates,
corrected a reference-group contamination, added penetrance). Therefore this amendment (a) **reclassifies all v9
analyses to date as DISCOVERY**, and (b) **pre-specifies a NEW confirmatory test of the REVISED direction, on an
independent outcome, before it is run.**

## 2. Reclassification of prior work

All v9 analyses run 2026-07-06 — `copi_sting_analyze_v9`, `irak_rna_v9`, `flagship_robustness_v9` (analysis-log runs
018–027; raw outputs in `RAW_OUTPUTS_v9.md`) — are **DISCOVERY / hypothesis-generating (Layer L0).** No confirmatory
claim attaches to any of them. The RNA∩PAV cohort (n=8,327; 46 R85H carriers) is too small to partition and is
designated the **discovery set in full.**

## 3. Revised hypotheses

- **H1 (CONFIRMATORY, directional — REVERSED from v5):** In the long-read cohort, cis-phased **AQ-strict STING1
  increases** odds of COPA-spectrum disease among R85H carriers — interaction **OR > 1**.
- **H2 (secondary):** IRAK3 loss-of-function increases COPA-spectrum odds (main effect; the parallel innate-brake route).
- **Positive control (GATING):** SAVI (STING1 gain-of-function) carriers show elevated ISG (validates that the 6-gene
  score reads STING *activation*). **If this positive control fails, the STING-ISG mechanism claim is withheld**
  regardless of the phenome result.

## 4. Confirmatory analysis plan (pre-specified; locked at this commit)

- **Cohort:** full long-read PAV cohort (~13,252; R85H carriers ~138). Outcome is EHR phenotype — a *different data
  modality* from the discovery outcome (ISG).
- **Exposure:** R85H carrier (`cb_variant_to_person`, vid `19-18911007-C-T`).
- **Modifier:** cis-phased STING1 **AQ-strict** (`AQ_d>0`) from PAV; HAQ modeled jointly (common WT reference).
- **PRIMARY outcome:** COPA-spectrum composite = any of **inflammatory arthritis** (ICD-10 M05/M06/M07/M08; ICD-9 714),
  **ILD** (J84/515/516/J98.4/D86), **vasculitis** (M31.3/M31.7/M30/M31/446/I77.6). Proband-blind, standard ICD;
  **osteoarthritis (M15–M19 / 715) excluded by construction.**
- **Secondary outcomes:** the three domains separately.
- **Model:** logistic `COPA ~ R85H * AQ + 16 ancestry PCs + age + sex` (+ ancestry category, pooled). Primary estimand
  = **R85H:AQ interaction OR.** Also the within-R85H-carrier `COPA ~ AQ` test.
- **Direction / threshold:** one-sided **AQ-risk (OR>1)** given the directional discovery prior; **α = 0.05**
  (two-sided also reported).
- **Independence sensitivity (pre-specified):** repeat the primary on long-read participants **NOT in the RNA discovery
  set** (~net-new subset) for cleaner discovery/confirmatory separation.
- **Small-cell / power rule:** AoU suppression (<20 → suppressed). **If R85H+AQ+ cases < 5, the interaction is reported
  as UNDERPOWERED / descriptive only** (point estimate + CI, no confirmatory claim).
- **Engine note:** OLS/logistic + PC adjustment (not the v5-registered SAIGE); relatedness via AoU unrelated set + PCs.
  A SAIGE re-run is a planned **sensitivity**, not the primary.
- **Execution:** `phenome_r85h_sting_v9.py` at this commit SHA. **Run ONCE.** Technical fixes to the phenotype query
  (OMOP/ICD schema) are permitted; **model, endpoints, covariates, and direction are frozen** — no changes after the
  result is seen.

## 5. Interpretation rules (pre-committed, before the result)

1. **AQ-risk OR>1, p<0.05, AND SAVI validates** → AQ-risk conditional-penetrance model **CONFIRMED on an independent
   clinical endpoint** (direction now established; explicitly reverses v5).
2. **Underpowered (cases <5)** → directionally-consistent-but-underpowered; **not confirmed**; flagged for a
   larger/replication cohort.
3. **Null or opposite** → the discovery AQ-risk signal **does not replicate** on the clinical endpoint; reported as such.
4. **SAVI fails** → STING-ISG mechanism validation withheld regardless of the phenome result.

## 6. Honesty ledger (limitations disclosed)

- The lock is a git+GitHub server-side push timestamp prior to running the confirmatory — **internal**, not a formal
  public registration. Mirror to OSF before manuscript.
- **Partial independence only:** the confirmatory outcome (EHR) is independent of the discovery outcome (ISG), but the
  long-read cohort *contains* the RNA discovery participants (same STING genotypes). True independence needs a separate
  cohort; the §4 net-new sensitivity is the best available mitigation within AoU, and this overlap is disclosed.
- The v5-registered SAIGE engine and v9-net-new held-out partition were **not** used for the RNA discovery arm; this
  amendment does not retroactively claim them.
