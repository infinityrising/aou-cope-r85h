#!/usr/bin/env python
"""AoU v9 -- BY-ANCESTRY stratified + META-ANALYSIS (ZS: 'look at ALL backgrounds'). Tests the 4 core exposures
(R85H, AQ, HAQ, IRAK3-LoF) in EVERY ancestry (afr/eur/amr/eas/sas/mid) for the key endpoints, then meta-analyzes across
strata (fixed + random effects) with Cochran-Q / I^2 heterogeneity. Purpose: (1) principled de-confounder -- a signal that
replicates across >=2 ancestries with low I^2 is not a founder/structure artifact (the run-063 neg-control-leakage fix);
(2) biology -- is STING-gain->fibrosis (and R85H's effect) UNIVERSAL or ancestry-background-modulated?
Endpoints: fibrotic-ILD incidence (AQ primary), bronchiectasis WITHIN-CASE severity (R85H+AQ; RESPIRATORY-ONLY hard outcome,
NO renal), asthma (IRAK3), + neg-controls (osteoarthritis/appendicitis) as the per-ancestry leakage diagnostic. Firth per
stratum (stable under separation); full grid reported (underpowered cells flagged, never omitted). Allele-freq table = the
power map. Standard app. Ends 'run complete'/'run failed'.
"""
import os, gzip, ast, json
import numpy as np, pandas as pd
from scipy.stats import chi2, norm
import warnings; warnings.filterwarnings("ignore"); np.seterr(all="ignore")
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); SNP=f"{MNT}/v9/wgs/short_read/snpindel"
ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
R85H_VID='19-18911007-C-T'; STINGV={'R71H':'5-139481493-C-T','G230A':'5-139478340-C-G','R293Q':'5-139477397-C-T'}
# endpoints
INCID={'fibrotic_ILD':['J84.1','515','516.3'],'asthma':['J45','493'],'osteoarthritis':['M15','M16','M17','M18','M19','715'],'appendicitis':['K35','540']}
NEG={'osteoarthritis','appendicitis'}
BRONCH=['J47','494']; RESP_HARD=['J96','518.81','518.83','518.84','Z99.81','V46.2','Z94.2','J95.85','Z99.11']  # RESPIRATORY-ONLY (no renal) -- the clean severity outcome
ANCS=['afr','eur','amr','eas','sas','mid']
MINCC=3  # min carrier-cases per stratum to estimate
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
    return b,np.sqrt(np.clip(np.diag(Iinv),0,None))
def firth_fit(dd,ex,covs):
    d2=dd.dropna(subset=[ex,'__y']+covs).copy(); ncc=int(((d2[ex]==1)&(d2['__y']==1)).sum())
    if ncc<MINCC or d2[ex].sum()<MINCC or d2['__y'].sum()<MINCC: return None
    cols=[ex]+[c for c in covs if d2[c].nunique()>1]
    X=np.column_stack([np.ones(len(d2))]+[d2[c].values.astype(float) for c in cols]); y=d2['__y'].values.astype(float)
    try:
        b,se=firth(X,y); be,s=float(b[1]),float(se[1])
        if not np.isfinite(be) or not np.isfinite(s) or s>10: return None
        return {'beta':be,'se':s,'ncc':ncc,'n':int(len(d2)),'or':round(float(np.exp(be)),3),'ci':[round(float(np.exp(be-1.96*s)),3),round(float(np.exp(be+1.96*s)),3)]}
    except Exception: return None
def meta(strata):
    b=np.array([s['beta'] for s in strata]); se=np.array([s['se'] for s in strata]); w=1/se**2; k=len(b); df=k-1
    bFE=float((w*b).sum()/w.sum()); seFE=float((1/w.sum())**0.5); Q=float((w*(b-bFE)**2).sum())
    I2=max(0.0,(Q-df)/Q)*100 if (Q>0 and df>0) else 0.0; pQ=float(chi2.sf(Q,df)) if df>0 else None
    c=w.sum()-(w**2).sum()/w.sum(); tau2=max(0.0,(Q-df)/c) if (df>0 and c>0) else 0.0
    wR=1/(se**2+tau2); bRE=float((wR*b).sum()/wR.sum()); seRE=float((1/wR.sum())**0.5)
    orci=lambda bb,ss:{'OR':round(float(np.exp(bb)),3),'CI':[round(float(np.exp(bb-1.96*ss)),3),round(float(np.exp(bb+1.96*ss)),3)],'p':round(float(2*norm.sf(abs(bb/ss))),4)}
    return {'k':k,'FE':orci(bFE,seFE),'RE':orci(bRE,seRE),'I2':round(I2,1),'Q_p':(round(pQ,4) if pQ is not None else None)}
S={}
try:
    import pysam
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    def caseset(codes):
        lk=" OR ".join([f"c.concept_code LIKE '{p}%'" for p in codes])
        return set(str(r[0]) for r in bq.query(f"SELECT co.person_id FROM `{PROJ}.{DS}.condition_occurrence` co JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id WHERE c.vocabulary_id LIKE 'ICD%' AND ({lk}) GROUP BY co.person_id HAVING COUNT(*)>=2"))
    def anyset(codes):
        lk=" OR ".join([f"c.concept_code LIKE '{p}%'" for p in codes])
        return set(str(r[0]) for r in bq.query(f"SELECT DISTINCT co.person_id FROM `{PROJ}.{DS}.condition_occurrence` co JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id WHERE c.vocabulary_id LIKE 'ICD%' AND ({lk})"))
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
    st={k:carr([v]) for k,v in STINGV.items()}; C_aq=(st['G230A']&st['R293Q'])-st['R71H']; C_haq=st['R71H']&st['G230A']&st['R293Q']
    EXPO={'R85H':carr([R85H_VID]),'AQ':C_aq,'HAQ':C_haq,'IRAK3_LoF':carr(irak3)}
    print("exposure carriers:",{k:len(v) for k,v in EXPO.items()})
    anc=pd.read_csv(ANC,sep="\t"); anc['research_id']=anc.research_id.astype(str)
    d=anc[['research_id','ancestry_pred']].copy()
    for nm,cs in EXPO.items(): d[nm]=d.research_id.isin(cs).astype(float)
    PCS=[]
    try:
        P=anc['pca_features'].apply(lambda x: ast.literal_eval(x) if isinstance(x,str) else (x if isinstance(x,list) else []))
        k=min(min(map(len,P[P.map(len)>0])) if (P.map(len)>0).any() else 0,10)   # 10 PCs within-stratum (avoid overfit on small strata)
        pcf=pd.DataFrame([row[:k] for row in P],columns=[f'PC{i+1}' for i in range(k)]); pcf['research_id']=anc.research_id.values
        d=d.merge(pcf,on='research_id',how='left'); PCS=[f'PC{i+1}' for i in range(k)]
    except Exception: pass
    try:
        ppl=bq.query(f"SELECT person_id, DATE_DIFF(DATE '2025-01-01', DATE(birth_datetime), YEAR) age, CAST(sex_at_birth_concept_id AS STRING) sexc FROM `{PROJ}.{DS}.person`").to_dataframe()
        ppl['research_id']=ppl.person_id.astype(str); d=d.merge(ppl[['research_id','age','sexc']],on='research_id',how='left')
        d['agez']=(d.age-d.age.mean())/d.age.std(); d['sex_m']=(d.sexc.astype(str)==d.sexc.astype(str).mode().iloc[0]).astype(float)
    except Exception: d['agez']=0.0; d['sex_m']=0.0
    try:
        smk=anyset(['Z72.0','F17','305.1','V15.82']); d['smoker']=d.research_id.isin(smk).astype(float)
    except Exception: d['smoker']=0.0
    COVS=['agez','sex_m','smoker']+PCS
    present=[a for a in ANCS if (d.ancestry_pred==a).sum()>0]
    print("ancestry sizes:",{a:int((d.ancestry_pred==a).sum()) for a in present})
    # ---- allele-frequency / power map by ancestry ----
    print("\n== ALLELE-FREQUENCY / POWER MAP (carrier count [%] by ancestry) ==")
    S['freq']={}
    print(f"{'exposure':10s} "+" ".join(f"{a:>13s}" for a in present))
    for nm in EXPO:
        row=[]; S['freq'][nm]={}
        for a in present:
            sub=d[d.ancestry_pred==a]; nc=int(sub[nm].sum()); pct=100*nc/len(sub) if len(sub) else 0
            S['freq'][nm][a]={'carriers':nc,'pct':round(pct,3)}; row.append(f"{nc}({pct:.2f}%)")
        print(f"{nm:10s} "+" ".join(f"{x:>13s}" for x in row))
    # ---- INCIDENCE endpoints: per-ancestry + meta ----
    def run_grid(label, endpoint_case,restrict=None):
        print(f"\n== {label} : per-ancestry Firth OR[CI](carrier-cases) -> FE/RE meta, I^2 ==")
        S.setdefault('grid',{})[label]={}
        base=d if restrict is None else d[d.research_id.isin(restrict)]
        for nm in EXPO:
            strata=[]; cells={}
            for a in present:
                dd=base[base.ancestry_pred==a].copy(); dd['__y']=dd.research_id.isin(endpoint_case).astype(float)
                r=firth_fit(dd,nm,COVS); cells[a]=r
                if r: strata.append(r)
            m=meta(strata) if len(strata)>=2 else None
            S['grid'][label][nm]={'strata':cells,'meta':m}
            cellstr=" ".join(f"{a}:{(str(cells[a]['or'])+str(cells[a]['ci'])+'(n'+str(cells[a]['ncc'])+')') if cells[a] else '--'}" for a in present)
            mstr=(f"| FE {m['FE']['OR']}{m['FE']['CI']}(p{m['FE']['p']}) RE {m['RE']['OR']}(p{m['RE']['p']}) I2={m['I2']}% Qp={m['Q_p']}") if m else "| (meta n/a <2 strata)"
            flag='  <NEG' if label.split(':')[0].strip() in NEG else ''
            print(f"  {nm:10s} {cellstr} {mstr}{flag}")
    cases={ph:caseset(cc) for ph,cc in INCID.items()}
    for ph in INCID: run_grid(ph, cases[ph])
    # ---- WITHIN-CASE bronchiectasis SEVERITY (respiratory-only hard outcome, no renal) ----
    bronch=caseset(BRONCH); resp_hard=anyset(RESP_HARD)
    print(f"\n== bronchiectasis WITHIN-CASE SEVERITY (RESP-only hard outcome; among {len(bronch)} bronch cases) ==")
    S['grid']['bronch_severity']={}
    for nm in EXPO:
        strata=[]; cells={}
        for a in present:
            dd=d[(d.ancestry_pred==a)&(d.research_id.isin(bronch))].copy(); dd['__y']=dd.research_id.isin(resp_hard).astype(float)
            r=firth_fit(dd,nm,COVS); cells[a]=r
            if r: strata.append(r)
        m=meta(strata) if len(strata)>=2 else None
        S['grid']['bronch_severity'][nm]={'strata':cells,'meta':m}
        cellstr=" ".join(f"{a}:{(str(cells[a]['or'])+'(n'+str(cells[a]['ncc'])+')') if cells[a] else '--'}" for a in present)
        mstr=(f"| FE {m['FE']['OR']}{m['FE']['CI']}(p{m['FE']['p']}) RE {m['RE']['OR']}(p{m['RE']['p']}) I2={m['I2']}%") if m else "| (meta n/a)"
        print(f"  {nm:10s} {cellstr} {mstr}")
    print("\n== READ: does the effect REPLICATE across ancestries (consistent dir + I^2<50% + sig meta = real, de-confounded)? do NEG-controls stay null in every ancestry? is R85H benign-alone in ALL backgrounds? ==")
    print("\n===== ANCESTRY META (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1400:])
