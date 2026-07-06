#!/usr/bin/env python
"""AoU v9 COPE-R85H — PHASED STING haplotypes from long-read PAV (Phased Assembly Variants) -> IFN. STANDARD app.
PAV gives assembly-level cis-phase: GT hap1|hap2 consistent across positions, phased=True, PS=None (global phase).
=> TRUE HAQ/AQ cis-haplotypes, NOT the unphased carrier proxies that failed (run 016). Source: manifest
grch38_pav_vcf (13,252 samples; RNA∩ 8,327; R85H∩ 138). Reads the 4 STING SNPs per haplotype -> classify
HAQ/AQ/H/R232 (+R220H dose) -> per-person dosage. Then the SELF-VALIDATING LADDER (STING haplotype -> IFN;
AQ/HAQ hypomorph should LOWER IFN) + MODULATION (R85H × STING -> IFN). Within-AFR primary. Threaded PAV reads.
Ends 'run complete' / 'run failed'.
"""
import os, sys, subprocess, io, json, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import Counter
import numpy as np, pandas as pd
subprocess.run([sys.executable,'-m','pip','install','-q','pysam'],capture_output=True)
import pysam
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); GS="gs://vwb-aou-datasets-controlled/"
LR=f"{MNT}/v9/wgs/long_read"; RNADIR=f"{MNT}/v9/multiomics/rnaseq"
TPM=f"{RNADIR}/rsem/aou_rnaseq_20260415.rsem_genes_tpm.txt.gz"
ANC=f"{MNT}/v9/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1); R85H_VID='19-18911007-C-T'
STING=[('R71H',139481493,'C','T'),('G230A',139478340,'C','G'),('R293Q',139477397,'C','T'),('R220H',139478370,'C','T')]
IFN={'IFI27':'ENSG00000165949','IFI44L':'ENSG00000137959','IFIT1':'ENSG00000185745','ISG15':'ENSG00000187608','RSAD2':'ENSG00000134321','SIGLEC1':'ENSG00000088827'}
NW=12; S={}
def sh(c): return subprocess.run(['bash','-lc',c],capture_output=True,text=True).stdout
def gs2mnt(p): return (MNT+"/"+p[len(GS):]) if isinstance(p,str) and p.startswith(GS) else None
def fetch1(v,pos):
    for c in ('chr5','5'):
        try: return list(v.fetch(c,pos-1,pos))
        except (ValueError,KeyError): continue
    return []
def classify(h):
    R,G,Q=h['R71H'],h['G230A'],h['R293Q']
    if R and G and Q: return 'HAQ'
    if G and Q and not R: return 'AQ'
    if R and not G and not Q: return 'H'
    if not R and not G and not Q: return 'R232'
    return 'partial'
def hap_call(row):
    try:
        mp=gs2mnt(row['grch38_pav_vcf'])
        if not (mp and os.path.exists(mp)): return None
        v=pysam.VariantFile(mp); sm=list(v.header.samples); h1={}; h2={}
        for (nm,pos,ref,alt) in STING:
            a1=a2=0
            for rec in fetch1(v,pos):
                if rec.pos!=pos: continue
                alts=[str(x) for x in (rec.alts or [])]
                if alt not in alts: continue
                ai=alts.index(alt)+1; gt=rec.samples[sm[0]]['GT']
                if gt and len(gt)==2:
                    a1=1 if gt[0]==ai else 0; a2=1 if gt[1]==ai else 0
            h1[nm]=a1; h2[nm]=a2
        v.close(); c1,c2=classify(h1),classify(h2)
        return {'research_id':row['research_id'],'hap1':c1,'hap2':c2,
                'R71H_d':h1['R71H']+h2['R71H'],'G230A_d':h1['G230A']+h2['G230A'],
                'R293Q_d':h1['R293Q']+h2['R293Q'],'R220H_d':h1['R220H']+h2['R220H'],
                'HAQ_d':int(c1=='HAQ')+int(c2=='HAQ'),'AQ_d':int(c1=='AQ')+int(c2=='AQ')}
    except Exception: return None
try:
    man=pd.read_csv(f"{LR}/manifest.tsv",sep="\t"); man['research_id']=man.research_id.astype(str)
    man=man[man.grch38_pav_vcf.notna()]
    rna_ids=set(pd.read_csv(f"{RNADIR}/manifest.tsv",sep="\t",usecols=['research_id']).research_id.astype(str))
    work=man[man.research_id.isin(rna_ids)].drop_duplicates('research_id')
    print(f"PAV∩RNA to phase: {len(work)} (threads={NW})")
    rows=[]; t0=time.time(); done=0
    with ThreadPoolExecutor(max_workers=NW) as ex:
        futs=[ex.submit(hap_call,r) for _,r in work.iterrows()]
        for f in as_completed(futs):
            r=f.result(); done+=1
            if r: rows.append(r)
            if done==200 or done%1500==0: print(f"   {done}/{len(work)} ok={len(rows)} {(time.time()-t0)/done:.3f}s/samp -> ~{(time.time()-t0)/done*len(work)/60:.1f}min total")
    hap=pd.DataFrame(rows); print(f"phased {len(hap)} in {time.time()-t0:.0f}s")
    hap.to_csv(os.path.expanduser("~/sting_haplotypes_v9.csv"),index=False)
    S['n_phased']=len(hap); S['haplotype_freq']=dict(Counter(list(hap.hap1)+list(hap.hap2)))
    S['diplo_top']={f"{min(a,b)}/{max(a,b)}":n for (a,b),n in Counter(zip(hap.hap1,hap.hap2)).most_common(8)}
    print("hap freq:",S['haplotype_freq']); print("diplo:",S['diplo_top'])
    raw=sh(f'zcat "{TPM}" 2>/dev/null | awk -F"\\t" \'NR==1 || $1 ~ /{"|".join(IFN.values())}/\'')
    m=pd.read_csv(io.StringIO(raw),sep="\t"); m['ensg']=m.gene_id.str.split('.').str[0]
    m=m[m.ensg.isin(IFN.values())].set_index('ensg').drop(columns=['gene_id','transcript_id(s)'])
    ex2=m.T.astype(float); ex2.index=ex2.index.astype(str); lg=np.log2(ex2+1)
    ifn=((lg-lg.mean())/lg.std()).mean(axis=1); ifndf=pd.DataFrame({'research_id':ifn.index,'ifn':ifn.values})
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    r85h=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid='{R85H_VID}'"))
    anc=pd.read_csv(ANC,sep="\t",usecols=['research_id','ancestry_pred']); anc['research_id']=anc.research_id.astype(str)
    d=hap.merge(ifndf,on='research_id',how='inner').merge(anc,on='research_id',how='left'); d['R85H']=d.research_id.isin(r85h).astype(int)
    a=d[d.ancestry_pred=='afr']
    S['n_ifn_join']=len(d); S['n_afr']=len(a); S['R85H_afr_n']=int(a.R85H.sum()); S['afr_hap_freq']=dict(Counter(list(a.hap1)+list(a.hap2)))
    print(f"IFN-joined {len(d)} | AFR {len(a)} | R85H(AFR) {int(a.R85H.sum())} | AFR hap freq {S['afr_hap_freq']}")
    import statsmodels.formula.api as smf
    def ols(f,dd,term):
        try: r=smf.ols(f,data=dd.dropna(subset=['ifn'])).fit(); return {'beta':round(float(r.params[term]),4),'p':round(float(r.pvalues[term]),4),'n':int(r.nobs)}
        except Exception as ex3: return {'error':str(ex3)[:60]}
    print("== LADDER: STING -> IFN (AQ/HAQ hypomorph => NEGATIVE beta = validates the score) ==")
    for term in ['HAQ_d','AQ_d','R220H_d','R71H_d','G230A_d','R293Q_d']:
        S[f'ladder_afr_{term}']=ols(f'ifn ~ {term}',a,term); S[f'ladder_all_{term}']=ols(f'ifn ~ {term} + C(ancestry_pred)',d,term)
        print(f"   AFR ifn~{term}: {S[f'ladder_afr_{term}']} | ALL(+anc): {S[f'ladder_all_{term}']}")
    print("== MODULATION: R85H × STING -> IFN (AFR) ==")
    for hh in ['AQ_d','HAQ_d','R220H_d']:
        a2=a.copy(); a2['inter']=a2.R85H*a2[hh]; S[f'mod_{hh}']=ols(f'ifn ~ R85H + {hh} + inter',a2,'inter'); print(f"   R85H×{hh}: {S[f'mod_{hh}']}")
    S['R85H_main_afr']=ols('ifn ~ R85H',a,'R85H'); print("   R85H main (AFR):",S['R85H_main_afr'])
    print("\n===== PHASED STING×IFN (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1400:])
