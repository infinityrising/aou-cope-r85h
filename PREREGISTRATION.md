# Pre-registration summary — COPE-R85H / IRAK3 oligogenic hypothesis (AoU v9)

Concise, public summary of the governing analysis plan (full draft held privately as `preregistration_v5_DRAFT.md`).
**Confirmatory results are not produced until this is SHA-locked and run once on the confirmatory partition — no peeking.**

## Questions
- **Q1** — Is COPE p.R85H pathogenic *alone* across genetic ancestries? (Prior AoU v8: marginal null.)
- **Q2** — Does its pathogenicity *require a second hit* (oligogenic), e.g. an innate-immune (IRAK3/TLR/NF-κB) loss-of-function?

## Engine
SAIGE mixed model (genetic-relatedness matrix + genetic-ancestry PCs + SPA), never self-reported race.
Interaction estimand = **additive-scale RERI/AP** (not the multiplicative product term, which two additive main
effects can inflate). Primary de-confounded effect = **within-AFR-local-ancestry stratified** (not nested adjustment).

## Endpoints (tiered)
- **Primary:** the 6-gene type-I IFN score (RNA-seq) as the molecular pillar — *AFR R85H∩RNA-seq = 39*, moderate power —
  paired with the proband-blind, **externally-templated** COPA-spectrum PheRS (Watkin/Vece, not case-derived).
- **Secondary (FDR):** COPA composite + components; STING-haplotype decomposition; R85H×{STING-R220H, AQ, COPI-burden,
  IFIH1, **innate-pathway burden**} interactions; trans-ancestry replication (thin off-founder: EUR=23 → caps at haplotype-association).

## De-confounding (the founder tilt is real: R85H 82.5% AFR; PC1-10 does NOT collapse it)
The **decisive test** = the R85H × second-hit interaction, within-AFR, benchmarked against a pre-registered panel of
**≥20 frequency/ancestry-matched negative-control-variant PAIRS**, with SAVI-STING as a positive control and local ancestry
at both loci. A signal that collapses in the AFR-only refit, or sits inside the negative-control-pair band, is declared **null**.

## Oligogenic arm — falsifiability lock
ONE primary confirmatory second-hit mask (the innate-pathway rare-damaging burden). Directional decision, pre-committed:
declare oligogenic effect **only if** RERI>0 AND survives within-AFR AND beats the pair band AND SAVI fires — **else R85H = benign
African founder marker**. IRAK3 L210* is PM2-only (LoF-tolerant gene, LOEUF 1.26) → VUS/modifier, no "driver" language.

## Firewall / partition
Discovery = v8-overlap participants (exploratory). Confirmatory = v9-net-new. Recon + endpoint selection use the
**discovery partition only**; the primary is SHA-locked before any v9-net-new query; the confirmatory pipeline runs once.

## No data egress
Only code lives here. No participant-level data leaves the AoU workspace.
