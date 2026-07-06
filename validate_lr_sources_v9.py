#!/usr/bin/env python
"""AoU v9 — VALIDATE per-sample long-read sources for STING phasing (task1) + immune-gene assembly seq (task3).
The manifest (15,424 samples) is the path index; PHASED per-sample products exist for the FULL cohort:
  grch38_deepvariant_phased_vcf  (read-backed phased DV VCF, GRCh38)  -> STING cis-haplotypes (task1)
  grch38_pav_vcf                 (Phased Assembly Variants; ZS's assembly choice) -> immune-gene 2nd-hit (task3)
  assembly_hap1_fa / hap2_fa     (raw haplotype contigs)               -> literal full-length gene sequence (task3)
Confirms: FUSE-mount path mapping, that a source is PHASED at the 4 STING positions (GT/phased/PS = same block),
per-source availability, and POWER (RNA + R85H-carrier overlap with each LR source). Small validation, no full-scans.
Ends 'run complete' / 'run failed'.
"""
import os, sys, subprocess, json
import pandas as pd
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); GS="gs://vwb-aou-datasets-controlled/"
LR=f"{MNT}/v9/wgs/long_read"; RNADIR=f"{MNT}/v9/multiomics/rnaseq"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1); R85H_VID='19-18911007-C-T'
STING=[('R71H',139481493,'C','T'),('G230A',139478340,'C','G'),('R293Q',139477397,'C','T'),('R220H',139478370,'C','T')]
def gs2mnt(p): return (MNT+"/"+p[len(GS):]) if isinstance(p,str) and p.startswith(GS) else None
S={}
try:
    man=pd.read_csv(f"{LR}/manifest.tsv",sep="\t"); man['research_id']=man.research_id.astype(str)
    KEY=['grch38_deepvariant_phased_vcf','grch38_pav_vcf','assembly_hap1_fa','assembly_hap2_fa','assembly_hap1_aln2_hg38_bam','grch38_haplotagged_bam']
    print("== A. per-source availability (non-null / 15424) ==")
    for c in KEY: S[f'avail_{c}']=int(man[c].notna().sum()); print(f"   {c:34s}: {S[f'avail_{c}']}")
    print("== B. gs->mount existence (first non-null each) ==")
    for c in ['grch38_deepvariant_phased_vcf','grch38_pav_vcf','assembly_hap1_fa','assembly_hap1_aln2_hg38_bam']:
        v=man[c].dropna(); v=v[v.astype(str).str.startswith('gs://')]
        if len(v):
            mp=gs2mnt(v.iloc[0]); ok=bool(mp and os.path.exists(mp)); S[f'exists_{c}']=ok
            print(f"   {c:34s}: exists={ok}  ...{v.iloc[0][len(GS):][:70]}")
    print("== C. PHASED at STING? (DV-phased vs PAV; want phased=True, shared PS) ==")
    subprocess.run([sys.executable,'-m','pip','install','-q','pysam'],capture_output=True)
    import pysam
    def fetch1(v,pos):
        for c in ('chr5','5'):
            try: return list(v.fetch(c,pos-1,pos))
            except (ValueError,KeyError): continue
        return []
    def check(mp):
        v=pysam.VariantFile(mp); sm=list(v.header.samples); out={}
        for (nm,pos,ref,alt) in STING:
            hit='(no record)'
            for rec in fetch1(v,pos):
                g=rec.samples[sm[0]]; hit=f"GT={g['GT']} phased={g.phased} PS={g.get('PS')} alt={rec.alts}"
            out[nm]=hit
        return sm[0],out
    for c in ['grch38_deepvariant_phased_vcf','grch38_pav_vcf']:
        sub=man[man[c].notna()]; done=False
        for _,r in sub.head(6).iterrows():
            mp=gs2mnt(r[c])
            if not (mp and os.path.exists(mp)): continue
            try: s,out=check(mp)
            except Exception as ex: print(f"   {c}: <{str(ex)[:50]}>"); done=True; break
            print(f"   [{c}] sample {r['research_id']} ({r['center']}/{r['platform']}):")
            for nm in [x[0] for x in STING]: print(f"       {nm}: {out[nm]}")
            S[f'phase_{c}']={'sample':r['research_id'],**out}; done=True; break
        if not done: print(f"   {c}: (none accessible in first 6)")
    print("== D. POWER: RNA + R85H-carrier overlap per source ==")
    rna=set(pd.read_csv(f"{RNADIR}/manifest.tsv",sep="\t",usecols=['research_id']).research_id.astype(str))
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    r85h=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid='{R85H_VID}'"))
    S['r85h_total']=len(r85h)
    for c in ['grch38_deepvariant_phased_vcf','grch38_pav_vcf','assembly_hap1_fa']:
        ids=set(man[man[c].notna()].research_id)
        S[f'rna_cap_{c}']=len(ids&rna); S[f'r85h_cap_{c}']=len(ids&r85h)
        print(f"   {c:34s}: source {len(ids):5d} | RNA∩ {len(ids&rna):4d} | R85H∩ {len(ids&r85h):3d}")
    print("\n===== LR SOURCE VALIDATION (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
