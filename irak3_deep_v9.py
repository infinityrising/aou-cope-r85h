#!/usr/bin/env python
"""AoU v9 -- ELEVATE THE IRAK3 ARM to R85H-parallel rigor (ZS goal). IRAK3-LoF (IRAK-M truncating, the MyD88/NF-kB brake) is
the index proband's second hit. Test it the way we tested R85H: (1) INCIDENCE across the COPA/airway/autoimmune panel +
asthma (Balaci positive-control); (2) WITHIN-CASE SEVERITY (Firth, organ-failure hard outcome) -- is IRAK3-LoF also a
severity modifier?; (3) IRAK3-LoF × R85H incidence interaction (the digenic double-carriers). Covariate-adjusted; Firth for
small cells. Q: is IRAK3-LoF a co-headline modifier (severity and/or a distinct phenome) rather than a secondary route?
Standard app. Ends 'run complete' / 'run failed'.
"""
import os, gzip, ast, json
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, chi2
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); SNP=f"{MNT}/v9/wgs/short_read/snpindel"
ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
R85H_VID='19-18911007-C-T'
PANEL={'bronchiectasis':['J47','494'],'fibrotic_ILD':['J84.1','515','516.3'],'LIP':['J84.2'],'alveolar_hemorrhage':['R04.2','R04.89'],
 'infl_arthritis':['M05','M06','M08','714'],'ANCA_vasculitis':['M31.3','M31.7'],'glomerulonephritis':['N01','N03','N05','580','581','582','583'],
 'asthma':['J45','493'],'COPD':['J44','496'],'recurrent_pna':['J18','486']}
HARD=['J96','518.81','518.83','518.84','Z99.81','V46.2','Z99.2','N18.6','Z94','Z99.11']
MINCODES=int(os.environ.get("MINCODES","2"))
def fnum(x):
    try: return float(x)
    except (TypeError,ValueError): return 0.0
def firth(X,y,max_iter=200,tol=1e-8):
    n,p=X.shape; b=np.zeros(p)
    for _ in range(max_iter):
        eta=X@b; pi=np.clip(1/(1+np.exp(-eta)),1e-12,1-1e-12); W=pi*(1-pi); I=(X.T*W)@X
        try: Iinv=np.linalg.inv(I)
        except np.linalg.LinAlgError: Iinv=np.linalg.pinv(I)
        h=W*np.einsum('ij,jk,ik->i',X,Iinv,X); U=X.T@(y-pi+h*(0.5-pi)); step=Iinv@U; b=b+step
        if np.max(np.abs(step))<tol: break
    eta=X@b; pi=np.clip(1/(1+np.exp(-eta)),1e-12,1-1e-12)
    ll=np.sum(y*np.log(pi)+(1-y)*np.log(1-pi)); _,logdet=np.linalg.slogdet(I); return b,ll+0.5*logdet,np.sqrt(np.clip(np.diag(Iinv),0,None))
def firth_test(df,ex,cov):
    d=df.dropna(subset=[ex]+cov).copy(); cols=[ex]+[c for c in cov if d[c].nunique()>1]
    Xf=np.column_stack([np.ones(len(d))]+[d[c].values.astype(float) for c in cols]); y=d['__y'].values.astype(float)
    bf,pllf,sef=firth(Xf,y); Xn=np.column_stack([np.ones(len(d))]+[d[c].values.astype(float) for c in cols if c!=ex]); _,plln,_=firth(Xn,y)
    p=float(chi2.sf(2*(pllf-plln),1)); return round(float(np.exp(bf[1])),3),[round(float(np.exp(bf[1]-1.96*sef[1])),3),round(float(np.exp(bf[1]+1.96*sef[1])),3)],round(p,4),int(len(d))
S={}
try:
    import pysam
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    def q_ids(sql): return set(str(r[0]) for r in bq.query(sql))
    def caseset(codes):
        lk=" OR ".join([f"c.concept_code LIKE '{p}%'" for p in codes])
        return q_ids(f"SELECT co.person_id FROM `{PROJ}.{DS}.condition_occurrence` co JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id WHERE c.vocabulary_id LIKE 'ICD%' AND ({lk}) GROUP BY co.person_id HAVING COUNT(*)>={MINCODES}")
    def anyset(codes):
        lk=" OR ".join([f"c.concept_code LIKE '{p}%'" for p in codes])
        return q_ids(f"SELECT DISTINCT co.person_id FROM `{PROJ}.{DS}.condition_occurrence` co JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id WHERE c.vocabulary_id LIKE 'ICD%' AND ({lk})")
    with gzip.open(VAT,'rt') as fh: COLS=fh.readline().rstrip("\n").split("\t")
    IX={c:COLS.index(c) for c in ['vid','gene_symbol','gvs_afr_af','LoF']}
    tbx=pysam.TabixFile(VAT); irak3=[]
    for line in tbx.fetch('chr12',66183995,66259622):
        f=line.split("\t")
        if f[IX['gene_symbol']]=='IRAK3' and f[IX['LoF']]=='HC' and fnum(f[IX['gvs_afr_af']])<0.001: irak3.append(f[IX['vid']])
    def carr(vids):
        out=set()
        for i in range(0,len(vids),900):
            if not vids[i:i+900]: continue
            inl=",".join("'"+v+"'" for v in vids[i:i+900])
            out|=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid IN ({inl})"))
        return out
    C_irak3=carr(irak3); C_r85h=carr([R85H_VID]); hard=anyset(HARD)
    print(f"IRAK3-LoF vids {len(irak3)} -> carriers {len(C_irak3)} | R85H {len(C_r85h)} | double {len(C_irak3&C_r85h)}")
    anc=pd.read_csv(ANC,sep="\t"); anc['research_id']=anc.research_id.astype(str)
    d=anc[['research_id','ancestry_pred']].copy(); d['IRAK3']=d.research_id.isin(C_irak3).astype(float); d['R85H']=d.research_id.isin(C_r85h).astype(float); d['hard']=d.research_id.isin(hard).astype(int)
    PCS=[]
    try:
        P=anc['pca_features'].apply(lambda x: ast.literal_eval(x) if isinstance(x,str) else (x if isinstance(x,list) else []))
        k=min(min(map(len,P[P.map(len)>0])) if (P.map(len)>0).any() else 0,16)
        pcf=pd.DataFrame([row[:k] for row in P],columns=[f'PC{i+1}' for i in range(k)]); pcf['research_id']=anc.research_id.values
        d=d.merge(pcf,on='research_id',how='left'); PCS=[f'PC{i+1}' for i in range(k)]
    except Exception: pass
    try:
        ppl=bq.query(f"SELECT person_id, DATE_DIFF(DATE '2025-01-01', DATE(birth_datetime), YEAR) age, CAST(sex_at_birth_concept_id AS STRING) sexc FROM `{PROJ}.{DS}.person`").to_dataframe()
        ppl['research_id']=ppl.person_id.astype(str); d=d.merge(ppl[['research_id','age','sexc']],on='research_id',how='left')
        d['agez']=(d.age-d.age.mean())/d.age.std(); d['sex_m']=(d.sexc.astype(str)==d.sexc.astype(str).mode().iloc[0]).astype(float)
    except Exception: d['agez']=0.0; d['sex_m']=0.0
    covb=['agez','sex_m']+PCS
    import statsmodels.formula.api as smf
    def logit(fm,dd,t):
        try: r=smf.logit(fm,data=dd,missing='drop').fit(disp=0); return round(float(np.exp(r.params[t])),3),round(float(r.pvalues[t]),4)
        except Exception: return None,None
    pc=(" + "+" + ".join(PCS)) if PCS else ""
    print("\n== IRAK3-LoF: INCIDENCE (OR) + WITHIN-CASE SEVERITY (Firth) across the panel ==")
    print(f"{'phenotype':18s} {'ncase':>6} {'IRAK3_inc':>14} | {'IRAK3_SEV Firth OR[CI] p':>34}")
    S['panel']=[]
    for ph,cc in PANEL.items():
        cs=caseset(cc); d['Y']=d.research_id.isin(cs).astype(int)
        inc=logit(f'Y ~ IRAK3 + agez + C(sexc){pc}',d,'IRAK3')
        sub=d[d.Y==1].copy(); sub['__y']=sub.research_id.isin(hard).astype(float); ncar=int(sub.IRAK3.sum())
        sev=firth_test(sub,'IRAK3',covb) if ncar>=3 else (None,None,None,ncar)
        S['panel'].append({'phenotype':ph,'n_cases':len(cs),'IRAK3_inc_OR':inc[0],'IRAK3_inc_p':inc[1],'IRAK3_case_carr':ncar,'IRAK3_sev_firthOR':sev[0],'IRAK3_sev_CI':sev[1],'IRAK3_sev_p':sev[2]})
        sv=f"{sev[0]} {sev[1]} p{sev[2]} (n{ncar})" if sev[0] is not None else f"(n{ncar})"
        print(f"{ph:18s} {len(cs):>6} {str(inc[0])+'(p'+str(inc[1])+')':>14} | {sv:>34}")
    print("\n== IRAK3-LoF × R85H interaction (incidence): does the double-hit exceed additive? ==")
    for ph,cc in [('bronchiectasis',PANEL['bronchiectasis']),('fibrotic_ILD',PANEL['fibrotic_ILD']),('infl_arthritis',PANEL['infl_arthritis'])]:
        d['Y']=d.research_id.isin(caseset(cc)).astype(int)
        o,p=logit(f'Y ~ IRAK3*R85H + agez + C(sexc){pc}',d,'IRAK3:R85H')
        ndouble=int(((d.IRAK3==1)&(d.R85H==1)&(d.Y==1)).sum())
        print(f"   {ph:16s} IRAK3:R85H OR={o} p={p} | double-carrier cases={ndouble}")
        S.setdefault('interaction',{})[ph]={'IRAK3xR85H_OR':o,'p':p,'double_cases':ndouble}
    print("\n== READ: asthma = Balaci positive control (expect risk if IRAK-M-LoF acts like Balaci; protective = variant heterogeneity). Is IRAK3-LoF a SEVERITY modifier like R85H? ==")
    print("\n===== IRAK3 DEEP (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1300:])
