#!/usr/bin/env python
"""AoU v9 — PHASED STING haplotypes from long-read ASSEMBLY alignments (gold-standard cis-phase, per plan).
Reads each person's hap1/hap2 assembly-aligned-to-hg38 BAM at the 4 STING positions (pysam) -> true cis haplotypes
-> HAQ/AQ/H/R232/R220H diplotypes. VALIDATION MODE (subset, N_SAMPLE) to confirm method + measure timing before scaling.
BAMs read via the FUSE mount (whole controlled bucket mounted). Ends 'run complete' / 'run failed'.
"""
import os, sys, subprocess, json, time
from collections import Counter
import pandas as pd
subprocess.run([sys.executable,'-m','pip','install','-q','pysam'],capture_output=True)
import pysam
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9")
LR=f"{MNT}/v9/wgs/long_read"; RNADIR=f"{MNT}/v9/multiomics/rnaseq"
GS="gs://vwb-aou-datasets-controlled/"
STING=[('R71H',139481493,'C','T'),('G230A',139478340,'C','G'),('R293Q',139477397,'C','T'),('R220H',139478370,'C','T')]
N_SAMPLE=300
S={}
def gs2mnt(p): return (MNT+"/"+p[len(GS):]) if isinstance(p,str) and p.startswith(GS) else None
def allele_at(bam, pos1):
    try:
        for col in bam.pileup('chr5', pos1-1, pos1, truncate=True, min_base_quality=0, max_depth=60):
            if col.reference_pos==pos1-1:
                for p in col.pileups:
                    if p.query_position is not None:
                        return p.alignment.query_sequence[p.query_position].upper()
    except Exception:
        return None
    return None
def classify(h):
    a=lambda n: h.get(n)=='alt'
    if any(h.get(n) is None for n in ('R71H','G230A','R293Q')): return 'incomplete'
    if a('R220H'): return 'R220H'
    if a('R71H') and a('G230A') and a('R293Q'): return 'HAQ'
    if a('G230A') and a('R293Q') and not a('R71H'): return 'AQ'
    if a('R71H') and not a('G230A') and not a('R293Q'): return 'H'
    if not a('R71H') and not a('G230A') and not a('R293Q'): return 'R232'
    return 'other'
try:
    man=pd.read_csv(f"{LR}/manifest.tsv",sep="\t").dropna(subset=['assembly_hap1_aln2_hg38_bam','assembly_hap2_aln2_hg38_bam'])
    man['research_id']=man.research_id.astype('int64')
    rna=set(pd.read_csv(f"{RNADIR}/manifest.tsv",sep="\t",usecols=['research_id']).research_id.astype('int64'))
    pool=man[man.research_id.isin(rna)]
    print(f"LR-with-assembly {len(man)} | RNA∩(LR-with-asm) {len(pool)} | validating on {min(N_SAMPLE,len(pool))}")
    sub=pool.head(N_SAMPLE)
    hapc=Counter(); diplo=Counter(); ok=0; t0=time.time()
    for _,r in sub.iterrows():
        p1,p2=gs2mnt(r['assembly_hap1_aln2_hg38_bam']),gs2mnt(r['assembly_hap2_aln2_hg38_bam'])
        if not (p1 and os.path.exists(p1) and p2 and os.path.exists(p2)): continue
        try:
            b1=pysam.AlignmentFile(p1); b2=pysam.AlignmentFile(p2)
        except Exception: continue
        h1={}; h2={}
        for (name,pos,ref,alt) in STING:
            x1,x2=allele_at(b1,pos),allele_at(b2,pos)
            h1[name]='alt' if x1==alt else ('ref' if x1==ref else None)
            h2[name]='alt' if x2==alt else ('ref' if x2==ref else None)
        b1.close(); b2.close()
        c1,c2=classify(h1),classify(h2); hapc[c1]+=1; hapc[c2]+=1; diplo[tuple(sorted([c1,c2]))]+=1; ok+=1
    dt=time.time()-t0
    S['n_processed']=ok; S['sec_per_person']=round(dt/max(ok,1),3)
    S['haplotype_freq']=dict(hapc.most_common())
    S['diplotypes']={f"{a}/{b}":n for (a,b),n in diplo.most_common(12)}
    S['proj_full_min']=round(len(pool)*S['sec_per_person']/60,1)
    print(f"processed {ok} in {dt:.0f}s = {S['sec_per_person']}s/person -> full {len(pool)} ≈ {S['proj_full_min']} min")
    print("haplotype freq:", S['haplotype_freq'])
    print("diplotypes:", S['diplotypes'])
    print("\n===== PHASED-STING VALIDATION (paste back) ====="); print(json.dumps(S,indent=1)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
