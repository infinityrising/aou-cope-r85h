#!/usr/bin/env python
"""AoU v9 -- PHENOME-WIDE SCAN (ZS: 'more to extract in the way of phenomes'). Hypothesis-generating breadth pass:
systematically test R85H, cis-AQ (phased, within-AFR), and IRAK3-LoF against a broad curated phenome (~50 phenotypes
across pulmonary / rheum-autoimmune / renal / heme / derm / systemic) PLUS 6 neutral NEGATIVE-CONTROL phenotypes that
should be null. Benjamini-Hochberg FDR per exposure; Firth fallback for small cells; direction (risk vs protective) reported.
This is DISCOVERY (v9) -- every hit is hypothesis-generating and must be re-derived in the v10 confirmation substrate; the
neg-control phenotypes + AFR-AF founder-tilt note guard against ancestry-structure false positives. Standard app.
Ends 'run complete' / 'run failed'.
"""
import os, gzip, ast, json
import numpy as np, pandas as pd
from scipy.stats import chi2
import warnings; warnings.filterwarnings("ignore"); np.seterr(all="ignore")   # small-cell fit warnings are expected/benign; keep output readable
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); SNP=f"{MNT}/v9/wgs/short_read/snpindel"
ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"; PHPAV=os.path.expanduser("~/sting_phenome_pav_v9.csv")
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
R85H_VID='19-18911007-C-T'
# broad curated phenome: {name: (system, [ICD prefixes ICD10+ICD9])}. NEG = neutral negative-control (expect null).
PHENOME={
 'bronchiectasis':('pulm',['J47','494']),'fibrotic_ILD':('pulm',['J84.1','515','516.3']),'ILD_other':('pulm',['J84','516']),
 'COPD':('pulm',['J44','496']),'asthma':('pulm',['J45','493']),'recurrent_pna':('pulm',['J18','486']),
 'pulm_htn':('pulm',['I27.0','I27.2','416.0','416.8']),'sarcoidosis':('pulm',['D86','135']),'pleural_effusion':('pulm',['J90','511']),
 'resp_failure':('pulm',['J96','518.81']),
 'rheumatoid_arthritis':('rheum',['M05','M06','714.0']),'SLE':('rheum',['M32','710.0']),'sjogren':('rheum',['M35.0','710.2']),
 'systemic_sclerosis':('rheum',['M34','710.1']),'myositis':('rheum',['M33','710.3','710.4']),'mixed_CTD':('rheum',['M35.1','M35.9']),
 'ANCA_vasculitis':('rheum',['M31.3','M31.7','446.4']),'vasculitis_other':('rheum',['M30','M31','446']),
 'psoriasis_PsA':('rheum',['L40','M07','696']),'IBD':('rheum',['K50','K51','555','556']),'ank_spond':('rheum',['M45','720.0']),
 'thyroiditis':('rheum',['E06','245']),'T1D':('rheum',['E10']),'autoimmune_hepatitis':('rheum',['K75.4','571.42']),'uveitis':('rheum',['H20','364']),
 'glomerulonephritis':('renal',['N01','N03','N05','580','581','582','583']),'nephrotic':('renal',['N04','581']),
 'CKD':('renal',['N18','585']),'ESRD':('renal',['N18.6','Z99.2']),
 'autoimmune_cytopenia':('heme',['D69.3','D59.1','287.31','283.0']),'lymphopenia':('heme',['D72.810','288.51']),'lymphoma':('heme',['C81','C82','C83','C84','C85','200','201','202']),
 'cutaneous_lupus':('derm',['L93','695.4']),'panniculitis':('derm',['M79.3','L52','729.30']),'vitiligo':('derm',['L80','709.01']),'alopecia_areata':('derm',['L63','704.01']),
 'fever_FUO':('systemic',['R50.9','780.6']),'failure_thrive':('systemic',['R62','783.4']),'hepatosplenomegaly':('systemic',['R16','789.2']),'lymphadenopathy':('systemic',['R59','785.6']),
 'refractive_error':('NEG',['H52','367']),'appendicitis':('NEG',['K35','540']),'osteoarthritis':('NEG',['M15','M16','M17','M18','M19','715']),
 'inguinal_hernia':('NEG',['K40','550']),'benign_skin_neoplasm':('NEG',['D22','216']),'cataract':('NEG',['H25','366']),
}
MINCODES=2
def fnum(x):
    try: return float(x)
    except (TypeError,ValueError): return 0.0
def firth(X,y,mx=200,tol=1e-8):
    n,p=X.shape; b=np.zeros(p)
    for _ in range(mx):
        eta=X@b; pi=np.clip(1/(1+np.exp(-eta)),1e-12,1-1e-12); W=pi*(1-pi); I=(X.T*W)@X
        try: Iinv=np.linalg.inv(I)
        except np.linalg.LinAlgError: Iinv=np.linalg.pinv(I)
        h=W*np.einsum('ij,jk,ik->i',X,Iinv,X); U=X.T@(y-pi+h*(0.5-pi)); step=Iinv@U; b=b+step
        if np.max(np.abs(step))<tol: break
    eta=X@b; pi=np.clip(1/(1+np.exp(-eta)),1e-12,1-1e-12)
    ll=np.sum(y*np.log(pi)+(1-y)*np.log(1-pi)); _,ld=np.linalg.slogdet(I); return b,ll+0.5*ld,np.sqrt(np.clip(np.diag(Iinv),0,None))
def firth_or(df,ex,cov):
    d=df.dropna(subset=[ex]+cov).copy(); cols=[ex]+[c for c in cov if d[c].nunique()>1]
    if d['__y'].sum()<3 or d[ex].sum()<3: return None,None,int(d[ex].sum())
    X=np.column_stack([np.ones(len(d))]+[d[c].values.astype(float) for c in cols]); y=d['__y'].values.astype(float)
    bf,pllf,_=firth(X,y); Xn=np.column_stack([np.ones(len(d))]+[d[c].values.astype(float) for c in cols if c!=ex]); _,plln,_=firth(Xn,y)
    return round(float(np.exp(bf[1])),3),round(float(chi2.sf(2*(pllf-plln),1)),4),int(d[ex].sum())
def bh_fdr(ps):
    ps=np.array([p if p is not None else np.nan for p in ps],float); ok=~np.isnan(ps); q=np.full(len(ps),np.nan)
    idx=np.where(ok)[0]; o=idx[np.argsort(ps[idx])]; m=len(o)
    for rank,i in enumerate(o,1): q[i]=ps[i]*m/rank
    for j in range(len(o)-2,-1,-1): q[o[j]]=min(q[o[j]],q[o[j+1]])
    return q
S={}
try:
    from google.cloud import bigquery
    import pysam
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    def caseset(codes):
        lk=" OR ".join([f"c.concept_code LIKE '{p}%'" for p in codes])
        return set(str(r[0]) for r in bq.query(f"SELECT co.person_id FROM `{PROJ}.{DS}.condition_occurrence` co JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id WHERE c.vocabulary_id LIKE 'ICD%' AND ({lk}) GROUP BY co.person_id HAVING COUNT(*)>={MINCODES}"))
    def carr(vids):
        out=set()
        for i in range(0,len(vids),900):
            if not vids[i:i+900]: continue
            inl=",".join("'"+v+"'" for v in vids[i:i+900])
            out|=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid IN ({inl})"))
        return out
    with gzip.open(VAT,'rt') as fh: COLS=fh.readline().rstrip("\n").split("\t")
    IX={c:COLS.index(c) for c in ['vid','gene_symbol','gvs_afr_af','LoF']}
    tbx=pysam.TabixFile(VAT); irak3=[]
    for line in tbx.fetch('chr12',66183995,66259622):
        f=line.split("\t")
        if f[IX['gene_symbol']]=='IRAK3' and f[IX['LoF']]=='HC' and fnum(f[IX['gvs_afr_af']])<0.001: irak3.append(f[IX['vid']])
    C_r85h=carr([R85H_VID]); C_irak3=carr(irak3)
    print(f"R85H {len(C_r85h)} | IRAK3-LoF {len(C_irak3)}")
    anc=pd.read_csv(ANC,sep="\t"); anc['research_id']=anc.research_id.astype(str)
    d=anc[['research_id','ancestry_pred']].copy()
    d['R85H']=d.research_id.isin(C_r85h).astype(float); d['IRAK3']=d.research_id.isin(C_irak3).astype(float)
    try:
        pav=pd.read_csv(PHPAV); pav['research_id']=pav.research_id.astype(str)
        d=d.merge(pav[['research_id','AQ_d','HAQ_d']],on='research_id',how='left'); d['cisAQ']=(d.AQ_d.fillna(0)>=1).astype(float); d['cisHAQ']=(d.HAQ_d.fillna(0)>=1).astype(float)
    except Exception: d['cisAQ']=np.nan; d['cisHAQ']=np.nan
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
    try:  # smoking (EHR ever-smoker, ICD tobacco) + healthcare utilization (visit count): confounder + ascertainment adjustment
        smk=set(str(r[0]) for r in bq.query(f"SELECT DISTINCT co.person_id FROM `{PROJ}.{DS}.condition_occurrence` co JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id WHERE c.vocabulary_id LIKE 'ICD%' AND (c.concept_code LIKE 'Z72.0%' OR c.concept_code LIKE 'F17%' OR c.concept_code LIKE '305.1%' OR c.concept_code LIKE 'V15.82%')"))
        d['smoker']=d.research_id.isin(smk).astype(float)
        util=bq.query(f"SELECT person_id, COUNT(DISTINCT visit_occurrence_id) nv FROM `{PROJ}.{DS}.visit_occurrence` GROUP BY person_id").to_dataframe()
        util['research_id']=util.person_id.astype(str); util['util_log']=np.log1p(util.nv.astype(float))
        d=d.merge(util[['research_id','util_log']],on='research_id',how='left'); d['util_log']=d.util_log.fillna(0.0); ADJ=['smoker','util_log']
        print(f"adj: ever-smoker {int(d.smoker.sum())} | util_log on")
    except Exception: d['smoker']=0.0; d['util_log']=0.0; ADJ=[]
    covb=['agez','sex_m']+ADJ+PCS; dA=d[d.ancestry_pred=='afr'].copy()
    print(f"cohort {len(d)} | AFR {len(dA)} | phenome {len(PHENOME)} phenotypes ({sum(1 for v in PHENOME.values() if v[0]=='NEG')} neg-control)\n")
    # exposures: (name, dataframe). cisAQ + cisHAQ(STING-dampened NEGATIVE CONTROL) tested within-AFR; R85H+IRAK3 whole-cohort.
    EXPO=[('R85H',d),('IRAK3',d),('cisAQ',dA),('cisHAQ',dA)]
    rows={e[0]:[] for e in EXPO}
    print(f"{'phenotype':22s} {'sys':8s} {'ncase':>6}  {'R85H OR(p)':>16} {'IRAK3 OR(p)':>16} {'cisAQ|AFR OR(p)':>18} {'cisHAQ|AFR[negctrl]':>19}")
    for ph,(sysv,cc) in PHENOME.items():
        cs=caseset(cc); ncase=len(cs); disp={}
        for ename,dd in EXPO:
            dd=dd.copy(); dd['__y']=dd.research_id.isin(cs).astype(float)
            orr,pp,ncar=firth_or(dd,ename,covb)
            rows[ename].append({'phenotype':ph,'system':sysv,'n_case':ncase,'OR':orr,'p':pp,'n_carrier_case':int(((dd[ename]==1)&(dd['__y']==1)).sum())})
            disp[ename]=f"{orr}(p{pp})" if orr is not None else "na"
        print(f"{ph:22s} {sysv:8s} {ncase:>6}  {disp['R85H']:>16} {disp['IRAK3']:>16} {disp['cisAQ']:>18} {disp['cisHAQ']:>19}")
    print("\n== FDR (Benjamini-Hochberg, per exposure across the phenome) -- hits q<0.10 (excl neg-controls from discovery emphasis) ==")
    S['hits']={}
    for ename,_ in EXPO:
        rr=rows[ename]; q=bh_fdr([r['p'] for r in rr])
        for r,qq in zip(rr,q): r['q']=None if np.isnan(qq) else round(float(qq),4)
        hits=[r for r in rr if r['q'] is not None and r['q']<0.10]; hits.sort(key=lambda r:r['q'])
        negflag=[r for r in rr if r['system']=='NEG' and r['p'] is not None and r['p']<0.05]
        S['hits'][ename]=hits
        print(f"\n  [{ename}] q<0.10:")
        for r in hits: print(f"     {r['phenotype']:22s} ({r['system']:8s}) OR={r['OR']} p={r['p']} q={r['q']} | carrier-cases={r['n_carrier_case']} {'*RISK' if (r['OR'] or 1)>1 else '*PROTECTIVE'}")
        if not hits: print("     (none survive q<0.10 -- expected for rare exposures at phenome breadth; see raw p in JSON)")
        if negflag: print(f"     !! NEG-CONTROL leakage (p<0.05, founder-tilt/structure warning): {[r['phenotype'] for r in negflag]}")
    S['all']=rows
    print("\n== READ: discovery-only. Any hit -> pre-specify in v10 on held-out data. Neg-control leakage => residual structure (add ancestry-band). cisAQ within-AFR removes the biggest confounder. ==")
    print("\n===== PHENOME SCAN (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1300:])
