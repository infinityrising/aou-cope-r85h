#!/usr/bin/env python
"""AoU v9 COPE-R85H — STING-axis x IFN score (pre-registered MOLECULAR co-pillar + self-validating ladder).
STANDARD app (plink2 + BigQuery + pandas).
 LADDER (validation): STING alleles -> IFN. AQ hypomorph should LOWER IFN (validates the 6-gene score reads STING biology).
 MODULATION: R85H x STING -> IFN. Does the STING axis modulate R85H's IFN effect (AQ-protection / R220H-risk)?
Within-AFR. UNPHASED first pass (true HAQ/AQ cis-haplotypes need long-read). Ends 'run complete' / 'run failed'.
"""
import os, subprocess, glob, io, json
import numpy as np, pandas as pd
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
BASE=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel")
RNADIR=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9/v9/multiomics/rnaseq")
TPM=f"{RNADIR}/rsem/aou_rnaseq_20260415.rsem_genes_tpm.txt.gz"; ANC=f"{BASE}/aux/ancestry/ancestry_preds.tsv"
STING={'R71H':(139481493,'C','T'),'G230A':(139478340,'C','G'),'R293Q':(139477397,'C','T'),'R220H':(139478370,'C','T')}
IFN={'IFI27':'ENSG00000165949','IFI44L':'ENSG00000137959','IFIT1':'ENSG00000185745','ISG15':'ENSG00000187608','RSAD2':'ENSG00000134321','SIGLEC1':'ENSG00000088827'}
R85H_VID='19-18911007-C-T'; S={}
def sh(c): return subprocess.run(['bash','-lc',c],capture_output=True,text=True).stdout
try:
    print("== A. IFN score ==")
    raw=sh(f'zcat "{TPM}" 2>/dev/null | awk -F"\\t" \'NR==1 || $1 ~ /{"|".join(IFN.values())}/\'')
    m=pd.read_csv(io.StringIO(raw),sep="\t"); m['ensg']=m.gene_id.str.split('.').str[0]
    m=m[m.ensg.isin(IFN.values())].set_index('ensg').drop(columns=['gene_id','transcript_id(s)'])
    expr=m.T.astype(float); expr.index=expr.index.astype('int64')
    lg=np.log2(expr+1); ifn=((lg-lg.mean())/lg.std()).mean(axis=1)
    d=pd.DataFrame({'ifn':ifn,'research_id':ifn.index}); print(f"   {len(d)} samples")

    print("== B. STING dosages (plink2 exome chr5) ==")
    pref=sorted(glob.glob(f"{BASE}/exome/pgen/*chr5*.pgen"))[0][:-5]
    sh(f'plink2 --pfile "{pref}" --chr 5 --from-bp 139477000 --to-bp 139482000 --export A --out /tmp/sting 2>/tmp/plog')
    r5=pd.read_csv("/tmp/sting.raw",sep="\t"); st=pd.DataFrame({'research_id':r5['IID'].astype('int64')})
    for name,(pos,ref,alt) in STING.items():
        col=[c for c in r5.columns if str(pos) in c]
        if col:
            c=col[0]; counted=c.rsplit('_',1)[-1]; dose=pd.to_numeric(r5[c],errors='coerce')
            st[name]= dose if counted==alt else (2-dose)      # orient to ALT (variant) allele
        else: st[name]=np.nan; print(f"   {name}: NOT FOUND")
    print("   ALT dose means:", {k:round(float(st[k].mean()),3) for k in STING})

    print("== C. join R85H + ancestry ==")
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    r85h=set(int(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid='{R85H_VID}'"))
    anc=pd.read_csv(ANC,sep="\t",usecols=['research_id','ancestry_pred']); anc['research_id']=anc.research_id.astype('int64')
    d=d.merge(st,on='research_id',how='inner').merge(anc,on='research_id',how='left')
    d['R85H']=d.research_id.isin(r85h).astype(int)
    d['AQ']=((d.G230A>0)&(d.R293Q>0)&(d.R71H==0)).astype(int)   # AQ-strict proxy (unphased)
    print(f"   RNA∩STING {len(d)} | R85H {int(d.R85H.sum())} | AQ-proxy {int(d.AQ.sum())}")

    import statsmodels.formula.api as smf
    a=d[d.ancestry_pred=='afr'].dropna(subset=['ifn'])
    def ols(f,dd,term):
        try: r=smf.ols(f,data=dd).fit(); return {'beta':round(float(r.params[term]),4),'p':round(float(r.pvalues[term]),4),'n':int(len(dd))}
        except Exception as ex: return {'error':str(ex)[:60]}
    print("== D. LADDER: STING -> IFN (AFR). AQ/hypomorph NEGATIVE = validates the score ==")
    for v in ['R71H','G230A','R293Q','R220H','AQ']:
        S[f'ladder_{v}']=ols(f'ifn ~ {v}',a,v); print(f"   IFN ~ {v}: {S[f'ladder_{v}']}")
    print("== E. MODULATION: R85H x AQ -> IFN (AFR) ==")
    a2=a.copy(); a2['RxAQ']=a2.R85H*a2.AQ
    S['R85H_x_AQ']=ols('ifn ~ R85H + AQ + RxAQ',a2,'RxAQ')
    S['R85H_main']=ols('ifn ~ R85H',a,'R85H')
    print(f"   R85H×AQ: {S['R85H_x_AQ']} | R85H main (AFR): {S['R85H_main']}")
    print("\n===== STING-IFN SUMMARY (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1000:])
