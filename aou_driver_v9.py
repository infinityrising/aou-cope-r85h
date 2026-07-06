#!/usr/bin/env python
"""AoU v9 COPE-R85H / IRAK3 — consolidated RECON driver (Verily Workbench RW 2.0).
Upload to the AoU Jupyter file browser, then in a cell run:   %run aou_driver_v9.py
Bounded + fail-safe. Reproduces recon runs 001-009 and prints ONE compact SUMMARY to paste back.
RECON / DESCRIPTIVE ONLY — it stops before any confirmatory association test (the pre-registered
confirmatory analysis is SHA-locked first, then run once). Ends with 'run complete' / 'run failed'.
"""
import json, pandas as pd
from google.cloud import bigquery

# ---- constants resolved during recon (edit only if the workspace CDR/mount changes) ----
CDR   = "wb-silky-artichoke-2408.C2025Q4R6"
PROJ, DS = CDR.split(".", 1)
MNT   = "/home/jupyter/workspace/vwb-aou-datasets-controlled-v9/v9"
AUX   = f"{MNT}/wgs/short_read/snpindel/aux"
ANC_F = f"{AUX}/ancestry/ancestry_preds.tsv"
QC_F  = f"{AUX}/qc/flagged_samples.tsv"
REL_F = f"{AUX}/relatedness/relatedness_flagged_samples.tsv"
RNA_F = f"{MNT}/multiomics/rnaseq/manifest.tsv"
VIDS  = {'COPE_R85H':'19-18911007-C-T','STING_R71H':'5-139481493-C-T','STING_G230A':'5-139478340-C-G',
         'STING_R293Q':'5-139477397-C-T','STING_R220H':'5-139478370-C-T','IRAK3_L210*':'12-66217211-T-A'}

bq = bigquery.Client()
def q1(sql): return list(bq.query(sql))[0]
def carrier_ids(vid):
    return set(r.element for r in bq.query(
        f"SELECT DISTINCT e.element FROM `{PROJ}.{DS}.cb_variant_to_person`, UNNEST(person_ids.list) e WHERE vid='{vid}'"))
S = {}

try:
    print("== 0. environment ==")
    assert q1(f"SELECT COUNT(*) c FROM `{PROJ}.{DS}`.INFORMATION_SCHEMA.TABLES WHERE table_name='cb_variant_to_person'").c, \
        "CDR/cb_variant_to_person unreachable — check the workspace CDR attach"
    print(f"  CDR OK: {CDR}")

    print("== 1. carrier counts ==")
    for lab, vid in VIDS.items():
        S[lab] = q1(f"SELECT COUNT(DISTINCT e.element) n FROM `{PROJ}.{DS}.cb_variant_to_person`, "
                    f"UNNEST(person_ids.list) e WHERE vid='{vid}'").n
        print(f"  {lab:12s} {vid:20s} = {S[lab]:,}")

    print("== 2. ancestry + analysis-ready cohort ==")
    anc = pd.read_csv(ANC_F, sep="\t", usecols=['research_id', 'ancestry_pred'])
    anc['research_id'] = anc.research_id.astype('int64'); research = set(anc.research_id)
    def firstcol_ids(path):
        df = pd.read_csv(path, sep="\t"); s = pd.to_numeric(df[df.columns[0]], errors='coerce').dropna()
        return set(s.astype('int64')) if len(s) else set(df[df.columns[0]].astype(str))
    qc_ids, rel_ids = firstcol_ids(QC_F), firstcol_ids(REL_F)
    id_ok = len(qc_ids & research) > 0.5 * max(len(qc_ids), 1)
    ready = anc[~anc.research_id.isin(qc_ids)]; unrel = ready[~ready.research_id.isin(rel_ids)]
    S.update(cohort_called=len(anc), cohort_ready=len(ready), cohort_unrelated=len(unrel),
             id_scheme_ok=bool(id_ok), ready_by_anc=ready.ancestry_pred.value_counts().to_dict())
    print(f"  ancestry-called={len(anc):,} QC-passed={len(ready):,} unrelated={len(unrel):,} id_scheme_ok={id_ok}")
    print(f"  ready by ancestry: {S['ready_by_anc']}")

    print("== 3. R85H ancestry + RNA-seq gating ==")
    r85h = carrier_ids('19-18911007-C-T'); a = anc[anc.research_id.isin(r85h)]
    rna = set(pd.read_csv(RNA_F, sep="\t", usecols=['research_id']).research_id.astype('int64')); ar = a[a.research_id.isin(rna)]
    S.update(R85H_by_anc=a.ancestry_pred.value_counts().to_dict(), R85H_AFR_frac=round(a.ancestry_pred.eq('afr').mean(), 3),
             RNAseq_N=len(rna), R85H_RNAseq=len(ar), R85H_RNAseq_by_anc=ar.ancestry_pred.value_counts().to_dict(),
             R85H_ready_AFR=len(set(ready[ready.ancestry_pred.eq('afr')].research_id) & r85h))
    print(f"  R85H by ancestry: {S['R85H_by_anc']} (AFR {S['R85H_AFR_frac']*100:.1f}%)")
    print(f"  RNA-seq N={len(rna):,}  R85H∩RNA-seq={S['R85H_RNAseq']} by anc {S['R85H_RNAseq_by_anc']}")

    print("== 4. R85H het/hom zygosity (exome pgen, plink2) ==")
    import glob, subprocess, os
    EXP = f"{MNT}/wgs/short_read/snpindel/exome/pgen"
    pgens = glob.glob(f"{EXP}/*chr19*.pgen") or sorted(glob.glob(f"{EXP}/*.pgen"))
    if pgens:
        pref = sorted(pgens)[0][:-5]
        subprocess.run(['bash', '-lc',
            f'plink2 --pfile "{pref}" --chr 19 --from-bp 18911007 --to-bp 18911007 --export A --out /tmp/r85h 2>/tmp/plink.log'],
            capture_output=True, text=True)
        if os.path.exists("/tmp/r85h.raw"):
            rw = pd.read_csv("/tmp/r85h.raw", sep="\t")
            dose_col = rw.columns[-1]                       # last column = the variant dosage
            counted = dose_col.rsplit('_', 1)[-1]           # allele plink2 counted (suffix after last '_')
            v = rw[['IID', dose_col]].copy(); v.columns = ['IID', 'd']
            v['IID'] = pd.to_numeric(v['IID'], errors='coerce')
            v['d'] = pd.to_numeric(v['d'], errors='coerce')
            v = v.dropna()
            # R85H = the T (ALT) allele; if plink2 counted the ref/other allele, invert to count T
            v['alt'] = (v['d'] if counted == 'T' else 2 - v['d']).round().astype(int)
            m = v.merge(anc, left_on='IID', right_on='research_id', how='inner')
            overall = {int(k): int(n) for k, n in m['alt'].value_counts().items()}
            hom = m.loc[m['alt'] == 2].groupby('ancestry_pred').size().astype(int).to_dict()
            het = m.loc[m['alt'] == 1].groupby('ancestry_pred').size().astype(int).to_dict()
            S.update(R85H_counted_allele=counted, R85H_geno_overall=overall, R85H_hom_by_anc=hom,
                     R85H_het_by_anc=het, R85H_carriers_from_geno=int(overall.get(1, 0) + overall.get(2, 0)))
            print(f"  dose col {dose_col} | counted allele={counted} | overall {{R85H_alleles:N}}: {overall}")
            print(f"  carriers from geno (het+hom) = {S['R85H_carriers_from_geno']:,}  (should ≈ 3,135)")
            print(f"  ★ HOM R85H by ancestry: {hom}")
            print(f"  HET R85H by ancestry: {het}")
        else:
            print("  no .raw — plink2 log:",
                  subprocess.run(['bash', '-lc', 'tail -3 /tmp/plink.log'], capture_output=True, text=True).stdout[-200:])
    else:
        print("  no exome pgen found in", EXP)

    print("\n===== SUMMARY (paste this block back) =====")
    print(json.dumps(S, indent=1))
    print("run complete")
except Exception as e:
    print(f"run failed: {type(e).__name__}: {e}")
