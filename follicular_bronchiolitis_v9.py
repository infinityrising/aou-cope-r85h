#!/usr/bin/env python
"""AoU v9 -- FOLLICULAR-BRONCHIOLITIS / airway-lymphoid phenotype check. Addresses the phenotype-validity flag: COPA's
DOCUMENTED airway lesion is FOLLICULAR BRONCHIOLITIS + LIP (lymphoid small-airway), not classic bronchiectasis (Tsui 2018,
Vece 2016). The true lesion is ~uncodeable in EHR -- follicular bronchiolitis has NO dedicated ICD code, and LIP (J84.2) is
ultra-rare -- so J47-bronchiectasis has served as a pragmatic surrogate. Here we (1) test the closest CODEABLE proxies: LIP
(J84.2) alone, and an AIRWAY-LYMPHOID composite (LIP + bronchiectasis + lung-cyst + other/unspecified ILD), and (2) ask
whether the R85H severity signal is bronchiectasis-SPECIFIC or GENERALIZES to the broader airway-lymphoid phenotype closer
to the COPA lesion (generalizes -> strengthens COPA relevance; absent -> weakens it; bronchiectasis-only -> the surrogate
is doing real work but is narrower than the mechanism). Honest coding limitation reported. Standard app.
Ends 'run complete' / 'run failed'.
"""
import os, gzip, ast, json
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); SNP=f"{MNT}/v9/wgs/short_read/snpindel"
ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
R85H_VID='19-18911007-C-T'; STINGV={'R71H':'5-139481493-C-T','G230A':'5-139478340-C-G','R293Q':'5-139477397-C-T'}
# closest codeable proxies for the COPA airway-lymphoid lesion (follicular bronchiolitis is uncodeable)
PANEL={
 'LIP':['J84.2'],                                             # lymphoid interstitial pneumonia (nearest specific lesion)
 'airway_lymphoid':['J84.2','J47','494','J98.4','J84.89','J84.9'],   # LIP + bronchiectasis + cyst + other/unspec ILD
 'bronchiectasis':['J47','494'],                             # the surrogate we have used (comparator)
}
RESP_HARD=['J96','518.81','518.83','518.84','Z99.81','V46.2','Z94.2','J95.85','Z99.11']
MINCODES=int(os.environ.get("MINCODES","2"))
S={}
try:
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    def q_ids(sql): return set(str(r[0]) for r in bq.query(sql))
    def caseset(codes):
        lk=" OR ".join([f"c.concept_code LIKE '{p}%'" for p in codes])
        return q_ids(f"SELECT co.person_id FROM `{PROJ}.{DS}.condition_occurrence` co JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id WHERE c.vocabulary_id LIKE 'ICD%' AND ({lk}) GROUP BY co.person_id HAVING COUNT(*)>={MINCODES}")
    def anyset(codes):
        lk=" OR ".join([f"c.concept_code LIKE '{p}%'" for p in codes])
        return q_ids(f"SELECT DISTINCT co.person_id FROM `{PROJ}.{DS}.condition_occurrence` co JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id WHERE c.vocabulary_id LIKE 'ICD%' AND ({lk})")
    def carr(vids):
        out=set()
        for v in vids: out|=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid='{v}'"))
        return out
    C_r85h=carr([R85H_VID]); st={k:carr([v]) for k,v in STINGV.items()}
    C_aq=(st['G230A']&st['R293Q'])-st['R71H']; C_haq=st['R71H']&st['G230A']&st['R293Q']
    resp=anyset(RESP_HARD); cases={ph:caseset(cc) for ph,cc in PANEL.items()}
    print("== phenotype proxies for the COPA airway-lymphoid lesion ==")
    for ph in PANEL: print(f"   {ph:16s} cases {len(cases[ph]):>6}  respiratory-severe {len(cases[ph]&resp):>5}")
    anc=pd.read_csv(ANC,sep="\t"); anc['research_id']=anc.research_id.astype(str)
    d=anc[['research_id','ancestry_pred']].copy()
    d['R85H']=d.research_id.isin(C_r85h).astype(int); d['AQ']=d.research_id.isin(C_aq).astype(int); d['HAQ']=d.research_id.isin(C_haq).astype(int); d['sev']=d.research_id.isin(resp).astype(int)
    PCS=[]
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
        try: r=smf.logit(fm,data=dd,missing='drop').fit(disp=0); return round(float(np.exp(r.params[t])),3),round(float(r.pvalues[t]),4)
        except Exception: return None,None
    print("\n== R85H incidence + within-case severity, and AQ/HAQ incidence, across the airway-lymphoid proxies ==")
    print(f"{'phenotype':16s} {'ncase':>6} | {'R85H_inc':>13} | {'R85H_SEV(car/sev)':>22} | {'AQ_inc':>12} {'HAQ_inc':>12}")
    S['results']=[]
    for ph in PANEL:
        d['Y']=d.research_id.isin(cases[ph]).astype(int)
        inc=logit(f'Y ~ R85H + age + C(sexc){pc16}',d,'R85H'); aq=logit(f'Y ~ AQ + age + C(sexc){pc16}',d,'AQ'); haq=logit(f'Y ~ HAQ + age + C(sexc){pc16}',d,'HAQ')
        sub=d[d.Y==1].copy(); nc=int(sub.R85H.sum()); ns=int((sub.R85H*sub.sev).sum())
        if nc>=5: sev=logit(f'sev ~ R85H + age + C(sexc){pc5}',sub,'R85H')
        else:
            a=ns; b=nc-ns; c=int(((sub.R85H==0)&(sub.sev==1)).sum()); e=int(((sub.R85H==0)&(sub.sev==0)).sum())
            try: orr,pp=fisher_exact([[a,b],[c,e]]); sev=(round(float(orr),3),round(float(pp),4))
            except Exception: sev=(None,None)
        S['results'].append({'phenotype':ph,'n_cases':len(cases[ph]),'R85H_inc':inc,'R85H_sev':sev,'R85H_case_carr':nc,'R85H_severe':ns,'AQ_inc':aq,'HAQ_inc':haq})
        sv=f"{sev[0]}(p{sev[1]}) {nc}/{ns}" if sev[0] is not None else f"(n={nc})"
        print(f"{ph:16s} {len(cases[ph]):>6} | {str(inc[0]):>6}(p{str(inc[1]):>5}) | {sv:>22} | {str(aq[0]):>5}(p{str(aq[1]):>4}) {str(haq[0]):>5}(p{str(haq[1]):>4})")
    print("\n== READ ==")
    print("   LIP (J84.2) = nearest specific lesion but ultra-rare (n small -> Fisher). airway_lymphoid = broader proxy.")
    print("   If R85H severity holds on airway_lymphoid ~ as on bronchiectasis -> signal is not an artifact of the J47 code")
    print("   and tracks the COPA-type airway phenotype. If only on bronchiectasis -> the surrogate is narrower than the lesion.")
    print("   TRUE follicular bronchiolitis is uncodeable in EHR; this is the closest achievable test.")
    print("\n===== FOLLICULAR BRONCHIOLITIS / AIRWAY-LYMPHOID (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
