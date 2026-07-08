#!/usr/bin/env python
"""AoU v9 -- LEAVE-ONE-OUT / jackknife robustness of the phased R85H×cisAQ -> type-I ISG interaction (run 056: β+0.78,
p0.011, n=11). The appraisal flagged the crude precursor (run 025) as LOO-fragile (drop top contributor -> p0.058). This
drops each R85H+&cisAQ+ individual one at a time, refits the phased interaction, and reports the p-range (how many drops
keep p<0.05) + a RANK-transform (outlier-robust) refit. If p is fragile to any single individual, the centerpiece is
anecdote-driven and must be reported as such. PHASED (RNA∩PAV, cis). Standard app. Ends 'run complete'/'run failed'.
"""
import os, io, gzip, subprocess, json
import numpy as np, pandas as pd
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); RNADIR=f"{MNT}/v9/multiomics/rnaseq"; SNP=f"{MNT}/v9/wgs/short_read/snpindel"
TPM=f"{RNADIR}/rsem/aou_rnaseq_20260415.rsem_genes_tpm.txt.gz"; ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CELL=os.path.expanduser("~/cell_fractions_v9.csv"); PHPAV=os.path.expanduser("~/sting_phenome_pav_v9.csv")
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
R85H_VID='19-18911007-C-T'; STINGV={'R71H':'5-139481493-C-T','G230A':'5-139478340-C-G','R293Q':'5-139477397-C-T'}
ISG={'IFI27':'ENSG00000165949','IFI44L':'ENSG00000137959','IFIT1':'ENSG00000185745','ISG15':'ENSG00000187608','RSAD2':'ENSG00000134321','SIGLEC1':'ENSG00000088827'}
def sh(c): return subprocess.run(['bash','-lc',c],capture_output=True,text=True).stdout
def fnum(x):
    try: return float(x)
    except (TypeError,ValueError): return 0.0
S={}
try:
    pat="|".join(ISG.values()); raw=sh(f'zcat "{TPM}" 2>/dev/null | awk -F"\\t" \'NR==1 || $1 ~ /{pat}/\'')
    m=pd.read_csv(io.StringIO(raw),sep="\t"); m['ensg']=m.gene_id.str.split('.').str[0]
    m=m[m.ensg.isin(ISG.values())].drop_duplicates('ensg').set_index('ensg')
    dropc=[c for c in ['gene_id','transcript_id(s)'] if c in m.columns]; ex=m.drop(columns=dropc).T.astype(float)
    ex.index=ex.index.astype(str); lg=np.log2(ex+1); zz=(lg-lg.mean())/lg.std()
    d=pd.DataFrame(index=zz.index); d['typeI_ISG']=zz.mean(axis=1); d['research_id']=d.index.astype(str)
    pav=pd.read_csv(PHPAV); pav['research_id']=pav.research_id.astype(str)
    d=d.merge(pav[['research_id','AQ_d','HAQ_d']],on='research_id',how='inner')
    d['cisAQ']=(d.AQ_d>=1).astype(int); d['cisHAQ']=(d.HAQ_d>=1).astype(int); d.loc[d.cisHAQ==1,'cisAQ']=0
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    C_r85h=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid='{R85H_VID}'"))
    d['R85H']=d.research_id.isin(C_r85h).astype(int)
    anc=pd.read_csv(ANC,sep="\t"); anc['research_id']=anc.research_id.astype(str); COV=[]
    try:
        import ast
        P=anc['pca_features'].apply(lambda x: ast.literal_eval(x) if isinstance(x,str) else (x if isinstance(x,list) else []))
        k=min(min(map(len,P[P.map(len)>0])) if (P.map(len)>0).any() else 0,16)
        pcf=pd.DataFrame([row[:k] for row in P],columns=[f'PC{i+1}' for i in range(k)]); pcf['research_id']=anc.research_id.values
        d=d.merge(pcf,on='research_id',how='left'); COV+=[f'PC{i+1}' for i in range(k)]
    except Exception: pass
    try:
        ppl=bq.query(f"SELECT person_id, DATE_DIFF(DATE '2025-01-01', DATE(birth_datetime), YEAR) age, CAST(sex_at_birth_concept_id AS STRING) sexc FROM `{PROJ}.{DS}.person`").to_dataframe()
        ppl['research_id']=ppl.person_id.astype(str); d=d.merge(ppl[['research_id','age','sexc']],on='research_id',how='left'); COV+=['age','C(sexc)']
    except Exception: pass
    try:
        cf=pd.read_csv(CELL); cf['research_id']=cf.research_id.astype(str); cc=[c for c in cf.columns if c!='research_id']
        d=d.merge(cf,on='research_id',how='left'); COV+=cc
    except Exception: pass
    covs=(" + "+" + ".join(COV)) if COV else ""
    import statsmodels.formula.api as smf
    def inter_p(dd,ycol='typeI_ISG'):
        try:
            r=smf.ols(f'{ycol} ~ R85H*cisAQ + R85H*cisHAQ{covs}',data=dd,missing='drop').fit()
            return round(float(r.params.get('R85H:cisAQ',np.nan)),3),round(float(r.pvalues.get('R85H:cisAQ',np.nan)),4)
        except Exception: return None,None
    b0,p0=inter_p(d); S['full']={'beta':b0,'p':p0}
    carr=d[(d.R85H==1)&(d.cisAQ==1)]; ids=list(carr.research_id); print(f"RNA∩PAV {len(d)} | R85H&cisAQ carriers n={len(ids)} | FULL R85H:cisAQ β={b0} p={p0}")
    print("\n== LEAVE-ONE-OUT (drop each R85H+cisAQ+ individual, refit) ==")
    loo=[]
    for rid in ids:
        b,p=inter_p(d[d.research_id!=rid]); loo.append({'dropped':rid,'beta':b,'p':p})
        print(f"   drop {rid}: β={b} p={p}")
    ps=[x['p'] for x in loo if x['p'] is not None]; S['loo']=loo
    S['loo_summary']={'n_carriers':len(ids),'p_min':min(ps) if ps else None,'p_max':max(ps) if ps else None,'n_drops_p_lt_0.05':sum(1 for x in ps if x<0.05),'n_drops':len(ps)}
    print(f"\n   LOO p-range: {min(ps):.4f}–{max(ps):.4f} | drops keeping p<0.05: {sum(1 for x in ps if x<0.05)}/{len(ps)}")
    # rank-transform (outlier-robust)
    d['isg_rank']=d.typeI_ISG.rank()
    br,pr=inter_p(d,'isg_rank'); S['rank_transform']={'beta':br,'p':pr}
    print(f"\n== RANK-TRANSFORM (outlier-robust) R85H:cisAQ β={br} p={pr} ==")
    print("\n== VERDICT: if any single drop pushes p>0.05 OR rank-transform is ns -> the interaction is anecdote/outlier-fragile; report as suggestive only. ==")
    print("\n===== INTERACTION LOO (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
