#!/usr/bin/env python
"""AoU v9 — EXPLORE the long-read phasing SOURCE for STING haplotypes (per plan; ZS chose assembly gold-standard).
The jointcall VCF is UNPHASED (run 017) and the manifest assembly-BAM paths are stale. This maps what per-sample
cis-phased source actually exists + how many samples it covers, so the real phasing module targets the right file.
Checks: (A) manifest columns [enumerates every per-sample file type], (B) LR dir tree, (C) single_sample_vcf +
assembly dirs + counts, (D) is a per-sample VCF PHASED at the 4 STING positions, (E) RNA-intersection power.
No full-scans (bounded listdir/glob only). Ends 'run complete' / 'run failed'.
"""
import os, sys, glob, subprocess, json
import pandas as pd
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9")
LR=f"{MNT}/v9/wgs/long_read"; RNADIR=f"{MNT}/v9/multiomics/rnaseq"
STING=[('R71H',139481493,'C','T'),('G230A',139478340,'C','G'),('R293Q',139477397,'C','T'),('R220H',139478370,'C','T')]
S={}
def ls(d):
    try: return sorted(os.listdir(d))
    except Exception as e: return [f"<{type(e).__name__}:{str(e)[:40]}>"]
try:
    print("== A. LR manifest columns (enumerates per-sample file types) ==")
    man=pd.read_csv(f"{LR}/manifest.tsv",sep="\t"); print("   n rows:", len(man)); print("   cols:", list(man.columns))
    S['manifest_cols']=list(man.columns); S['manifest_n']=len(man)
    pcols=[c for c in man.columns if man[c].astype(str).str.contains('gs://').any()]
    print("   gs://-path cols:", pcols)
    for c in pcols[:12]:
        ex=man[c].dropna().astype(str); ex=ex[ex.str.startswith('gs://')]
        if len(ex): print(f"     {c}: {ex.iloc[0]}")
    print("== B. LR dir tree (depth 2) ==")
    for d in ls(LR):
        p=f"{LR}/{d}"
        if os.path.isdir(p): print(f"   {d}/ -> {ls(p)[:8]}")
    print("== C. locate single_sample_vcf + assembly ==")
    cands=glob.glob(f"{LR}/single_sample_vcf")+glob.glob(f"{LR}/*/single_sample_vcf")+glob.glob(f"{LR}/*/*/single_sample_vcf")
    acands=glob.glob(f"{LR}/assembly")+glob.glob(f"{LR}/*/assembly")+glob.glob(f"{LR}/*/*/assembly")
    print("   single_sample_vcf:", cands); print("   assembly:", acands)
    S['ssv_dirs']=cands; S['asm_dirs']=acands
    ssv_ids=set(); asm_ids=set()
    if cands:
        ssv=cands[0]; samps=ls(ssv); ssv_ids=set(x for x in samps if not x.startswith('<'))
        print(f"   {ssv}: {len(ssv_ids)} samples, e.g. {samps[:3]}"); S['ssv_n']=len(ssv_ids)
        s0=samps[0]; print(f"   {s0}/GRCh38: {ls(f'{ssv}/{s0}/GRCh38')}")
    if acands:
        asm_ids=set(x for x in ls(acands[0]) if not x.startswith('<')); S['asm_n']=len(asm_ids)
        print(f"   {acands[0]}: {len(asm_ids)} samples, e.g. {sorted(asm_ids)[:3]}")
    print("== D. is a per-sample VCF PHASED at STING? ==")
    vcf=[]
    if cands: vcf=glob.glob(f"{cands[0]}/{ls(cands[0])[0]}/GRCh38/*.vcf.gz")
    if vcf:
        subprocess.run([sys.executable,'-m','pip','install','-q','pysam'],capture_output=True)
        import pysam
        v=pysam.VariantFile(vcf[0]); sm=list(v.header.samples); print(f"   VCF={os.path.basename(vcf[0])} sample={sm}")
        S['ssv_phased']={}
        for (nm,pos,ref,alt) in STING:
            hit='(no record)'
            try:
                for rec in v.fetch('chr5',pos-1,pos):
                    g=rec.samples[sm[0]]; hit=f"GT={g['GT']} phased={g.phased} ref={rec.ref} alt={rec.alts}"; S['ssv_phased'][nm]=bool(g.phased)
            except Exception as ex: hit=f"<{str(ex)[:40]}>"
            print(f"   {nm} @chr5:{pos}: {hit}")
    else: print("   (no per-sample VCF found at expected path)")
    print("== E. RNA-intersection power (which phased source is well-powered?) ==")
    rna=set(pd.read_csv(f"{RNADIR}/manifest.tsv",sep="\t",usecols=['research_id']).research_id.astype(str))
    S['rna_n']=len(rna)
    if ssv_ids: S['ssv_cap_rna']=len({i.split('.')[0] for i in ssv_ids} & rna); print(f"   RNA {len(rna)} | single_sample_vcf {len(ssv_ids)} | ∩ {S['ssv_cap_rna']}")
    if asm_ids: S['asm_cap_rna']=len({i.split('.')[0] for i in asm_ids} & rna); print(f"   RNA {len(rna)} | assembly {len(asm_ids)} | ∩ {S['asm_cap_rna']}")
    print("\n===== LR-PHASE EXPLORE (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
