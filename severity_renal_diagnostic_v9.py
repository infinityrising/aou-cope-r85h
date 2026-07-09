#!/usr/bin/env python
"""AoU v9 -- DECISIVE DIAGNOSTIC (per 3 adversarial critics): is the 'bronchiectasis SEVERITY' signal actually NEPHROLOGY?
The renal-contaminated severity outcome (Version A: RESP_HARD + N18.6/Z99.2/N17/Z94/V42.0) counts ESRD/dialysis/AKI/
transplant as 'severe bronchiectasis'. Because R85H/AQ/IRAK3 are SYSTEMIC interferonopathy genes with renal pleiotropy AND
R85H is ~82% AFR (APOL1 renal-risk confound), including renal codes can manufacture a false 'severity' signal.
TESTS, among bronchiectasis cases only: (1) what FRACTION of the contaminated 'severe' pool qualifies ONLY via renal/
transplant codes (no respiratory hard code)?; (2) are R85H/AQ/HAQ/IRAK3-LoF carriers ENRICHED in that PURE-RENAL stratum
(within-AFR + PCs)? -> if yes, the severity finding is nephrology; (3) side-by-side: exposure->severe under CONTAMINATED
(renal-incl) vs CLEAN (respiratory-only, >=2 codes) -- how much does the OR collapse? Also flags the in-sample SENTINEL
(1348881, R85H+AQ+COPA) the prereg required excluding. Standard app. Ends 'run complete'/'run failed'.
"""
import os, gzip, ast, json
import numpy as np, pandas as pd
import warnings; warnings.filterwarnings("ignore"); np.seterr(all="ignore")
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); SNP=f"{MNT}/v9/wgs/short_read/snpindel"
ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
R85H_VID='19-18911007-C-T'; STINGV={'R71H':'5-139481493-C-T','G230A':'5-139478340-C-G','R293Q':'5-139477397-C-T'}
SENTINEL='1348881'  # in-sample R85H+AQ+COPA; prereg required exclusion (deviation D4)
BRONCH=['J47','494']
RESP_HARD=['J96.1','J96.2','518.83','Z99.81','V46.2','Z94.2','Z99.11']   # CLEAN: chronic-respiratory-only (drop acute J96.0/518.81 + postproc J95.85 + renal)
RENAL_TX=['N18.6','Z99.2','N17','V42.0','Z94.0','Z94']                    # the contamination: ESRD/dialysis/AKI/kidney-tx/any-tx
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
def ffit(dd,ex,covs):
    d2=dd.dropna(subset=[ex,'__y']+covs).copy(); ncc=int(((d2[ex]==1)&(d2['__y']==1)).sum())
    if ncc<3 or d2[ex].sum()<3 or d2['__y'].sum()<3: return {'or':None,'p':None,'ncc':ncc}
    cols=[ex]+[c for c in covs if d2[c].nunique()>1]
    X=np.column_stack([np.ones(len(d2))]+[d2[c].values.astype(float) for c in cols]); y=d2['__y'].values.astype(float)
    try:
        b,se=firth(X,y); z=b[1]/se[1]
        from scipy.stats import norm
        return {'or':round(float(np.exp(b[1])),3),'ci':[round(float(np.exp(b[1]-1.96*se[1])),3),round(float(np.exp(b[1]+1.96*se[1])),3)],'p':round(float(2*norm.sf(abs(z))),4),'ncc':ncc}
    except Exception: return {'or':None,'p':None,'ncc':ncc}
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
    bronch=caseset(BRONCH); resp=anyset(RESP_HARD); renal=anyset(RENAL_TX)
    bc=bronch; resp_sev=bronch&resp; renal_pos=bronch&renal; pure_renal=(bronch&renal)-resp
    contaminated=resp_sev|renal_pos      # Version-A 'severe' (renal-incl)
    print(f"bronchiectasis cases {len(bc)}")
    print(f"  CONTAMINATED 'severe' (resp OR renal) = {len(contaminated)}")
    print(f"  respiratory-hard 'severe' (clean)      = {len(resp_sev)}")
    print(f"  renal/transplant-positive              = {len(renal_pos)}")
    print(f"  ★ PURE-RENAL 'severe' (renal+, resp-)  = {len(pure_renal)}  => {100*len(pure_renal)/max(len(contaminated),1):.1f}% of the contaminated 'severe' pool has NO respiratory hard code")
    S['counts']={'bronch':len(bc),'contaminated_severe':len(contaminated),'resp_severe':len(resp_sev),'renal_pos':len(renal_pos),'pure_renal_severe':len(pure_renal),'pct_pure_renal':round(100*len(pure_renal)/max(len(contaminated),1),1)}
    anc=pd.read_csv(ANC,sep="\t"); anc['research_id']=anc.research_id.astype(str)
    d=anc[['research_id','ancestry_pred']].copy()
    for nm,cs in EXPO.items(): d[nm]=d.research_id.isin(cs).astype(float)
    PCS=[]
    try:
        P=anc['pca_features'].apply(lambda x: ast.literal_eval(x) if isinstance(x,str) else (x if isinstance(x,list) else []))
        k=min(min(map(len,P[P.map(len)>0])) if (P.map(len)>0).any() else 0,10)
        pcf=pd.DataFrame([row[:k] for row in P],columns=[f'PC{i+1}' for i in range(k)]); pcf['research_id']=anc.research_id.values
        d=d.merge(pcf,on='research_id',how='left'); PCS=[f'PC{i+1}' for i in range(k)]
    except Exception: pass
    try:
        ppl=bq.query(f"SELECT person_id, DATE_DIFF(DATE '2025-01-01', DATE(birth_datetime), YEAR) age, CAST(sex_at_birth_concept_id AS STRING) sexc FROM `{PROJ}.{DS}.person`").to_dataframe()
        ppl['research_id']=ppl.person_id.astype(str); d=d.merge(ppl[['research_id','age','sexc']],on='research_id',how='left')
        d['agez']=(d.age-d.age.mean())/d.age.std(); d['sex_m']=(d.sexc.astype(str)==d.sexc.astype(str).mode().iloc[0]).astype(float)
    except Exception: d['agez']=0.0; d['sex_m']=0.0
    COVS=['agez','sex_m']+PCS
    dcase=d[d.research_id.isin(bc)].copy(); dA=dcase[dcase.ancestry_pred=='afr'].copy()
    print(f"\n== among bronchiectasis cases (within-AFR n={len(dA)}): is each exposure ENRICHED in PURE-RENAL 'severe' vs in CLEAN respiratory 'severe'? ==")
    print("   (if exposure -> PURE-RENAL is significant, the 'severity modifier' signal is NEPHROLOGY, not lung)")
    S['enrichment']={}
    for nm in EXPO:
        dA['__y']=dA.research_id.isin(pure_renal).astype(float); pr=ffit(dA,nm,COVS)
        dA['__y']=dA.research_id.isin(resp_sev).astype(float); rs=ffit(dA,nm,COVS)
        dA['__y']=dA.research_id.isin(contaminated).astype(float); ct=ffit(dA,nm,COVS)
        S['enrichment'][nm]={'pure_renal':pr,'resp_clean':rs,'contaminated':ct}
        print(f"   {nm:10s} PURE-RENAL OR={pr['or']}(p{pr['p']},n{pr['ncc']}) | CLEAN-resp OR={rs['or']}(p{rs['p']},n{rs['ncc']}) | CONTAMINATED OR={ct['or']}(p{ct['p']},n{ct['ncc']})")
    print("\n== SENTINEL check ==")
    insamp={nm:(SENTINEL in EXPO[nm]) for nm in EXPO}
    print(f"   in-sample sentinel {SENTINEL}: carrier of {[k for k,v in insamp.items() if v]} | bronchiectasis case: {SENTINEL in bc} | in contaminated-severe: {SENTINEL in contaminated}")
    S['sentinel']={'id':SENTINEL,'carrier_of':[k for k,v in insamp.items() if v],'is_bronch_case':SENTINEL in bc,'in_severe':SENTINEL in contaminated}
    print("\n== READ: (1) high %pure-renal => the 'severe' pool is contaminated by non-lung disease. (2) exposure->PURE-RENAL significant => that exposure's 'severity' signal is renal pleiotropy/APOL1, NOT bronchiectasis severity. (3) CLEAN-resp OR << CONTAMINATED OR => the headline was inflated by renal codes. ==")
    print("\n===== SEVERITY RENAL DIAGNOSTIC (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1400:])
