#!/usr/bin/env python
"""AoU v9 -- reconcile the AQ contradiction the appraisal flagged: UNPHASED AQ -> fibrotic-ILD incidence is OR 1.232
p=0.0015 (run 040), but CIS-PHASED AQ_d -> lung is NULL (OR 1.007, run 041). Which is real? Put the three views side by
side on the SAME fibrotic-ILD endpoint: (a) whole-cohort unphased AQ (reproduce), (b) WITHIN-AFR unphased AQ (does holding
ancestry constant kill it? -> if yes, the unphased signal was founder/ancestry tagging), (c) cis-phased AQ_d on the PAV
cohort (the gold standard, power-limited). Contrast with the SYSTEMIC endpoint (which WAS phase-confirmed, OR 1.365) as an
internal positive control. Verdict rule: if within-AFR unphased AQ attenuates to ~1, the fibrotic-ILD 'AQ incidence' claim
is ancestry-confounded and must be dropped; the phase-confirmed claim is AQ->SYSTEMIC only. Standard app. Ends 'run complete'/'run failed'.
"""
import os, ast, json
import numpy as np, pandas as pd
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); SNP=f"{MNT}/v9/wgs/short_read/snpindel"
ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
PHPAV=os.path.expanduser("~/sting_phenome_pav_v9.csv")
STINGV={'R71H':'5-139481493-C-T','G230A':'5-139478340-C-G','R293Q':'5-139477397-C-T'}
FIBROTIC=['J84.1','515','516.3']; SYSTEMIC=['M05','M06','M08','714','M30','M31','N01','N03','N05','580','581','582','583']
MINCODES=2
S={}
try:
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    def caseset(codes):
        lk=" OR ".join([f"c.concept_code LIKE '{p}%'" for p in codes])
        return set(str(r.person_id) for r in bq.query(f"SELECT co.person_id FROM `{PROJ}.{DS}.condition_occurrence` co JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id WHERE c.vocabulary_id LIKE 'ICD%' AND ({lk}) GROUP BY co.person_id HAVING COUNT(*)>={MINCODES}"))
    def carr(vids):
        out=set()
        for v in vids:
            out|=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid='{v}'"))
        return out
    st={k:carr([v]) for k,v in STINGV.items()}; C_aq=(st['G230A']&st['R293Q'])-st['R71H']; C_haq=st['R71H']&st['G230A']&st['R293Q']
    fib=caseset(FIBROTIC); sysd=caseset(SYSTEMIC)
    print(f"fibrotic-ILD cases {len(fib)} | systemic cases {len(sysd)} | AQ(unphased) {len(C_aq)} | HAQ(unphased) {len(C_haq)}")
    anc=pd.read_csv(ANC,sep="\t"); anc['research_id']=anc.research_id.astype(str)
    d=anc[['research_id','ancestry_pred']].copy(); d['AQ']=d.research_id.isin(C_aq).astype(int); d['HAQ']=d.research_id.isin(C_haq).astype(int)
    d['fib']=d.research_id.isin(fib).astype(int); d['sys']=d.research_id.isin(sysd).astype(int)
    try:
        P=anc['pca_features'].apply(lambda x: ast.literal_eval(x) if isinstance(x,str) else (x if isinstance(x,list) else []))
        k=min(min(map(len,P[P.map(len)>0])) if (P.map(len)>0).any() else 0,16)
        pcf=pd.DataFrame([row[:k] for row in P],columns=[f'PC{i+1}' for i in range(k)]); pcf['research_id']=anc.research_id.values
        d=d.merge(pcf,on='research_id',how='left'); PCS=[f'PC{i+1}' for i in range(k)]
    except Exception: PCS=[]
    try:
        ppl=bq.query(f"SELECT person_id, DATE_DIFF(DATE '2025-01-01', DATE(birth_datetime), YEAR) age, CAST(sex_at_birth_concept_id AS STRING) sexc FROM `{PROJ}.{DS}.person`").to_dataframe()
        ppl['research_id']=ppl.person_id.astype(str); d=d.merge(ppl[['research_id','age','sexc']],on='research_id',how='left')
    except Exception: pass
    pc16=(" + "+" + ".join(PCS)) if PCS else ""; pc5=(" + "+" + ".join(PCS[:5])) if PCS else ""
    import statsmodels.formula.api as smf
    def logit(fm,dd,t):
        try: r=smf.logit(fm,data=dd,missing='drop').fit(disp=0); return round(float(np.exp(r.params[t])),3),round(float(r.pvalues[t]),4),int(dd[t].sum())
        except Exception: return None,None,int(dd[t].sum())
    afr=d[d.ancestry_pred=='afr']
    print("\n== UNPHASED AQ vs HAQ(negative control) -> endpoint ==")
    S['unphased']={}
    for ep in ['fib','sys']:
        wc=logit(f'{ep} ~ AQ + age + C(sexc){pc16}',d,'AQ'); wa=logit(f'{ep} ~ AQ + age + C(sexc){pc5}',afr,'AQ')
        hwc=logit(f'{ep} ~ HAQ + age + C(sexc){pc16}',d,'HAQ'); hwa=logit(f'{ep} ~ HAQ + age + C(sexc){pc5}',afr,'HAQ')
        S['unphased'][ep]={'AQ_whole':wc,'AQ_afr':wa,'HAQ_whole':hwc,'HAQ_afr':hwa}
        print(f"   {ep}: AQ whole OR={wc[0]}(p{wc[1]}) AFR OR={wa[0]}(p{wa[1]}) | HAQ[neg-ctrl] whole OR={hwc[0]}(p{hwc[1]}) AFR OR={hwa[0]}(p{hwa[1]})")
    # ---- cis-phased AQ_d on PAV ----
    print("\n== CIS-PHASED AQ_d -> endpoint (PAV cohort) ==")
    try:
        pav=pd.read_csv(PHPAV); pav['research_id']=pav.research_id.astype(str)
        p=pav[['research_id','AQ_d','HAQ_d']].merge(d[['research_id','fib','sys','ancestry_pred']+PCS[:5]],on='research_id',how='left')
        p=p.merge(d[['research_id','age','sexc']],on='research_id',how='left')
        S['phased']={}
        for ep in ['fib','sys']:
            ph=logit(f'{ep} ~ AQ_d + HAQ_d + age + C(sexc){pc5}',p,'AQ_d')
            phh=logit(f'{ep} ~ AQ_d + HAQ_d + age + C(sexc){pc5}',p,'HAQ_d')
            S['phased'][ep]={'cis_AQ_d':ph,'cis_HAQ_d':phh,'n_cases_PAV':int(p[ep].sum())}
            print(f"   {ep}: cis-AQ_d OR={ph[0]} p={ph[1]} | cis-HAQ_d[neg-ctrl] OR={phh[0]} p={phh[1]} | PAV cases {int(p[ep].sum())}")
    except FileNotFoundError: print("   (sting_phenome_pav_v9.csv not found)"); S['phased']='(csv missing)'
    print("\n== VERDICT RULE ==")
    print("   If UNPHASED AQ->fib is significant whole-cohort but ATTENUATES within-AFR -> the fibrotic-ILD 'AQ incidence'")
    print("   signal is ancestry-confounded (founder tilt), NOT a real AQ effect -> drop it; keep only the phase-confirmed")
    print("   AQ->SYSTEMIC claim. If it HOLDS within-AFR, AQ->fibrotic-ILD is real and phasing is merely power-limited in PAV.")
    print("\n===== AQ PHASING RECONCILE (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
