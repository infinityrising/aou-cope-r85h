#!/usr/bin/env python
"""AoU v9 — PHASE 2: COPI+STING phased variants (from copi_sting_extract) -> ISG stratified by COPI-mutation × STING.
Loads ~/copi_sting_pav_v9.csv (phased per-patient from PHASE 1), joins ISG (6-gene RSEM) + ancestry + R85H(BQ), then:
 LADDER (validation): HAQ/AQ -> ISG (hypomorph => NEGATIVE). FLAGSHIP: R85H × HAQ 2x2 + interaction (HAQ protective
 => R85H+HAQ+ lowest, interaction negative). BROADER: any damaging-COPI × HAQ. Fast + re-runnable (no PAV reads).
Can run on a partial CSV mid-extraction (reports n so far). STANDARD app. Ends 'run complete' / 'run failed'.
"""
import os, subprocess, io, json
from collections import Counter
import numpy as np, pandas as pd
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9")
RNADIR=f"{MNT}/v9/multiomics/rnaseq"; SNP=f"{MNT}/v9/wgs/short_read/snpindel"
TPM=f"{RNADIR}/rsem/aou_rnaseq_20260415.rsem_genes_tpm.txt.gz"; ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1); R85H_VID='19-18911007-C-T'
CSV=os.path.expanduser("~/copi_sting_pav_v9.csv")
IFN={'IFI27':'ENSG00000165949','IFI44L':'ENSG00000137959','IFIT1':'ENSG00000185745','ISG15':'ENSG00000187608','RSAD2':'ENSG00000134321','SIGLEC1':'ENSG00000088827'}
S={}
def sh(c): return subprocess.run(['bash','-lc',c],capture_output=True,text=True).stdout
try:
    g=pd.read_csv(CSV); g['research_id']=g.research_id.astype(str)
    for c in ['HAQ_d','AQ_d','R220H_d','R85H_d','copi_dmg']: g[c]=pd.to_numeric(g[c],errors='coerce').fillna(0).astype(int)
    S['n_loaded']=len(g); S['hap_freq']=dict(Counter(list(g.hap1)+list(g.hap2)))
    S['R85H_from_PAV']=int((g.R85H_d>0).sum()); S['copi_dmg_carriers']=int((g.copi_dmg>0).sum()); S['novel_carriers']=int((g.pav_novel>0).sum())
    print(f"loaded {len(g)} | hap freq {S['hap_freq']} | R85H(PAV) {S['R85H_from_PAV']} | COPI-dmg {S['copi_dmg_carriers']}")
    print("== ISG (6-gene) + ancestry + R85H(BQ) ==")
    raw=sh(f'zcat "{TPM}" 2>/dev/null | awk -F"\\t" \'NR==1 || $1 ~ /{"|".join(IFN.values())}/\'')
    m=pd.read_csv(io.StringIO(raw),sep="\t"); m['ensg']=m.gene_id.str.split('.').str[0]
    m=m[m.ensg.isin(IFN.values())].set_index('ensg').drop(columns=['gene_id','transcript_id(s)'])
    ex=m.T.astype(float); ex.index=ex.index.astype(str); lg=np.log2(ex+1)
    isg=((lg-lg.mean())/lg.std()).mean(axis=1); isgdf=pd.DataFrame({'research_id':isg.index,'ifn':isg.values})
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    r85h_bq=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid='{R85H_VID}'"))
    anc=pd.read_csv(ANC,sep="\t",usecols=['research_id','ancestry_pred']); anc['research_id']=anc.research_id.astype(str)
    d=g.merge(isgdf,on='research_id',how='inner').merge(anc,on='research_id',how='left'); d['R85H_bq']=d.research_id.isin(r85h_bq).astype(int)
    S['R85H_PAV_vs_BQ']={'pav':int((d.R85H_d>0).sum()),'bq':int(d.R85H_bq.sum()),'concord':int(((d.R85H_d>0)==(d.R85H_bq>0)).sum())}
    d['R85H']=((d.R85H_d>0)|(d.R85H_bq>0)).astype(int); d['HAQ']=(d.HAQ_d>0).astype(int); d['COPImut']=((d.copi_dmg>0)|(d.R85H>0)).astype(int)
    a=d[d.ancestry_pred=='afr'].copy()
    S['n_ifn_join']=len(d); S['n_afr']=len(a); S['afr_R85H']=int(a.R85H.sum()); S['afr_HAQ']=int(a.HAQ.sum())
    print(f"   ISG-joined {len(d)} | AFR {len(a)} | R85H(AFR) {int(a.R85H.sum())} | HAQ(AFR) {int(a.HAQ.sum())} | R85H PAVvsBQ {S['R85H_PAV_vs_BQ']}")
    import statsmodels.formula.api as smf
    def ols(f,dd,term):
        try: r=smf.ols(f,data=dd.dropna(subset=['ifn'])).fit(); return {'beta':round(float(r.params[term]),4),'p':round(float(r.pvalues[term]),4),'n':int(r.nobs)}
        except Exception as e: return {'error':str(e)[:70]}
    print("== LADDER: STING -> ISG (AQ/HAQ hypomorph => NEGATIVE = validates) ==")
    for t in ['HAQ_d','AQ_d','R220H_d']:
        S[f'ladder_all_{t}']=ols(f'ifn ~ {t} + C(ancestry_pred)',d,t); S[f'ladder_afr_{t}']=ols(f'ifn ~ {t}',a,t)
        print(f"   ifn~{t}: ALL(+anc) {S[f'ladder_all_{t}']} | AFR {S[f'ladder_afr_{t}']}")
    print("== FLAGSHIP: R85H × HAQ -> ISG (HAQ protective => R85H+HAQ+ lowest; interaction NEGATIVE) ==")
    def cells(dd):
        out={}
        for r in (0,1):
            for h in (0,1):
                sub=dd[(dd.R85H==r)&(dd.HAQ==h)].ifn.dropna()
                out[f'R85H{r}_HAQ{h}']={'mean':round(float(sub.mean()),3) if len(sub) else None,'n':int(len(sub))}
        return out
    S['cells_afr']=cells(a); S['cells_all']=cells(d)
    print("   AFR 2x2:",S['cells_afr']); print("   ALL 2x2:",S['cells_all'])
    xa=a[a.R85H==1]
    S['haq_protection_afr']={'R85H+HAQ+_mean':round(float(xa[xa.HAQ==1].ifn.mean()),3) if (xa.HAQ==1).any() else None,'n+':int((xa.HAQ==1).sum()),
                             'R85H+HAQ-_mean':round(float(xa[xa.HAQ==0].ifn.mean()),3) if (xa.HAQ==0).any() else None,'n-':int((xa.HAQ==0).sum())}
    S['inter_afr']=ols('ifn ~ R85H*HAQ',a,'R85H:HAQ'); S['inter_all_adj']=ols('ifn ~ R85H*HAQ + C(ancestry_pred)',d,'R85H:HAQ')
    print("   HAQ-protection within R85H+ (AFR):",S['haq_protection_afr'])
    print("   interaction R85H:HAQ  AFR",S['inter_afr']," | ALL(+anc)",S['inter_all_adj'])
    print("== BROADER: any damaging COPI × HAQ -> ISG ==")
    S['inter_COPImut_all_adj']=ols('ifn ~ COPImut*HAQ + C(ancestry_pred)',d,'COPImut:HAQ')
    xc=d[d.COPImut==1]
    S['COPImut_HAQ_means']={'mut+HAQ+':round(float(xc[xc.HAQ==1].ifn.mean()),3) if (xc.HAQ==1).any() else None,'n+':int((xc.HAQ==1).sum()),
                            'mut+HAQ-':round(float(xc[xc.HAQ==0].ifn.mean()),3) if (xc.HAQ==0).any() else None,'n-':int((xc.HAQ==0).sum())}
    print("   COPImut×HAQ (ALL+anc):",S['inter_COPImut_all_adj']," means",S['COPImut_HAQ_means'])
    print("\n===== COPI×STING×ISG (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except FileNotFoundError:
    print(f"run failed: {CSV} not found — run copi_sting_extract_v9.py first")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1400:])
