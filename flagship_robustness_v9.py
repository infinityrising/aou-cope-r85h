#!/usr/bin/env python
"""AoU v9 — ROBUSTNESS of the R85H × AQ-strict ISG signal. The pre-reg PRIMARY (AQ-strict PROTECTIVE) INVERTED:
AQ-strict is significantly POSITIVE (RISK) in R85H carriers (AFR interaction +1.05 p0.006 / ALL-adj +0.89 p0.002),
OPPOSITE the protective HAQ. Before believing an n=11, pre-reg-contradicting result, check whether it's a BROAD shift
or 1-2 outliers: per-individual ISG (+ per-gene, to see if it's a coordinated IFN signature), median vs mean, a rank
test (outlier-robust), and leave-one-out interaction stability. Focus = AFR R85H carriers. STANDARD app.
Ends 'run complete' / 'run failed'.
"""
import os, subprocess, io, json
import numpy as np, pandas as pd
from scipy import stats
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
    for c in ['HAQ_d','AQ_d','R220H_d']: g[c]=pd.to_numeric(g[c],errors='coerce').fillna(0).astype(int)
    raw=sh(f'zcat "{TPM}" 2>/dev/null | awk -F"\\t" \'NR==1 || $1 ~ /{"|".join(IFN.values())}/\'')
    m=pd.read_csv(io.StringIO(raw),sep="\t"); m['ensg']=m.gene_id.str.split('.').str[0]
    m=m[m.ensg.isin(IFN.values())].set_index('ensg').drop(columns=['gene_id','transcript_id(s)'])
    ex=m.T.astype(float); ex.index=ex.index.astype(str); lg=np.log2(ex+1); zz=(lg-lg.mean())/lg.std()
    isg=zz.mean(axis=1); inv={v:k for k,v in IFN.items()}; zz.columns=[inv[c] for c in zz.columns]
    zz['ifn']=isg; zz['research_id']=zz.index
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    r85h=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid='{R85H_VID}'"))
    anc=pd.read_csv(ANC,sep="\t",usecols=['research_id','ancestry_pred']); anc['research_id']=anc.research_id.astype(str)
    d=g.merge(zz,on='research_id').merge(anc,on='research_id',how='left'); d['R85H']=d.research_id.isin(r85h).astype(int)
    d['AQ']=(d.AQ_d>0).astype(int); d['HAQ']=(d.HAQ_d>0).astype(int)
    a=d[d.ancestry_pred=='afr'].copy()
    genes=list(IFN.keys())
    def cls(r): return 'AQ' if r.AQ else ('HAQ' if r.HAQ else 'neither')
    rp=a[a.R85H==1].copy(); rp['cls']=rp.apply(cls,axis=1)
    print(f"AFR R85H+ carriers = {len(rp)} (AQ+ {int((rp.cls=='AQ').sum())}, HAQ+ {int((rp.cls=='HAQ').sum())}, neither {int((rp.cls=='neither').sum())})")
    print("== individual ISG by STING class (sorted desc; per-gene z shows if it's a coordinated IFN signature) ==")
    for c in ['AQ','HAQ','neither']:
        sub=rp[rp.cls==c].sort_values('ifn',ascending=False)
        S[f'{c}']={'n':len(sub),'mean':round(float(sub.ifn.mean()),3),'median':round(float(sub.ifn.median()),3),
                   'max':round(float(sub.ifn.max()),3),'min':round(float(sub.ifn.min()),3),
                   'trim10':round(float(stats.trim_mean(sub.ifn,0.1)),3) if len(sub)>=5 else None,
                   'n_gt1SD':int((sub.ifn>1).sum())}
        print(f"-- R85H+/{c}: {S[f'{c}']} --")
        print(sub[['research_id','ifn']+genes].round(2).to_string(index=False))
    print("== rank test (outlier-robust): R85H+AQ+ vs R85H+non-AQ (AFR) ==")
    aqp=rp[rp.cls=='AQ'].ifn; oth=rp[rp.cls!='AQ'].ifn
    u,pu=stats.mannwhitneyu(aqp,oth,alternative='greater')
    S['MWU']={'p_AQ_greater':round(float(pu),4),'median_AQ':round(float(aqp.median()),3),'median_nonAQ':round(float(oth.median()),3)}
    print("   ",S['MWU'])
    print("== leave-one-out: crude R85H×AQ interaction (AFR) dropping each R85H+AQ+ individual ==")
    import statsmodels.formula.api as smf
    full=smf.ols('ifn ~ R85H*AQ',data=a).fit(); S['full']={'beta':round(float(full.params['R85H:AQ']),3),'p':round(float(full.pvalues['R85H:AQ']),4)}
    loo=[]
    for rid in rp[rp.cls=='AQ'].research_id:
        r=smf.ols('ifn ~ R85H*AQ',data=a[a.research_id!=rid]).fit(); loo.append((rid,round(float(r.params['R85H:AQ']),3),round(float(r.pvalues['R85H:AQ']),4)))
    ps=[x[2] for x in loo]; S['full_result']=S['full']; S['LOO_p_range']=[min(ps),max(ps)]; S['LOO_still_sig_dropping_any1']=bool(max(ps)<0.05)
    print(f"   full {S['full']} | LOO p-range [{min(ps):.4f},{max(ps):.4f}] | stays<0.05 dropping ANY single AQ+: {S['LOO_still_sig_dropping_any1']}")
    for x in sorted(loo,key=lambda z:-z[2]): print(f"     drop {x[0]}: beta {x[1]} p {x[2]}")
    print("\n===== AQ ROBUSTNESS (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1000:])
