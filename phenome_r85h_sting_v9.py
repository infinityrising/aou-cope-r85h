#!/usr/bin/env python
"""AoU v9 — PHENOME arm: R85H × STING-haplotype -> COPA-spectrum PHENOTYPE (full long-read cohort; EHR, not RNA).
The molecular arm capped at n=46 R85H∩RNA. Phenome uses the FULL long-read PAV cohort (~13,252; R85H∩ ~138) with
phased STING (sting_only_extract PHENOME=1 -> ~/sting_phenome_pav_v9.csv) + EHR — 3x the R85H carriers, phenotype not
RNA-gated. Q: does AQ-strict (the molecular risk-permissive STING) also raise COPA-spectrum DISEASE in R85H carriers?
Phenotype = standard ICD COPA-spectrum (ILD/arthritis/vasculitis; proband-BLIND, NOT sentinel-derived). Logistic,
PC/age/sex-adjusted; interaction + within-R85H-carrier test + prevalence table. STANDARD app. Ends 'run complete'/'run failed'.
"""
import os, ast, json
import numpy as np, pandas as pd
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); SNP=f"{MNT}/v9/wgs/short_read/snpindel"
ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1); R85H_VID='19-18911007-C-T'
CSV=os.path.expanduser("~/sting_phenome_pav_v9.csv")
PHENO={'ILD':['J84%','515%','516%','J98.4','D86%'],'arthritis':['M05%','M06%','M08%','714%'],'vasculitis':['M31.3%','M31.7%','M30%','M31%','446%','I77.6']}
S={}
def sc_of(r):
    h={r.hap1,r.hap2}
    if 'HAQ' in h and 'AQ' in h: return 'HAQ/AQ'
    if 'HAQ' in h: return 'HAQ'
    if 'AQ' in h: return 'AQ'
    return 'WT/other'
try:
    g=pd.read_csv(CSV); g['research_id']=g.research_id.astype(str)
    for c in ['HAQ_d','AQ_d','R220H_d']: g[c]=pd.to_numeric(g[c],errors='coerce').fillna(0).astype(int)
    print(f"loaded {len(g)} phenome cohort (all PAV) | hap freq {dict(pd.Series(list(g.hap1)+list(g.hap2)).value_counts())}")
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    r85h=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid='{R85H_VID}'"))
    def cases(codes):
        lk=" OR ".join([f"c.concept_code LIKE '{p}'" for p in codes])
        q=f"""SELECT DISTINCT co.person_id FROM `{PROJ}.{DS}.condition_occurrence` co
              JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id
              WHERE c.vocabulary_id LIKE 'ICD%' AND ({lk})"""
        return set(str(r.person_id) for r in bq.query(q))
    domsets={k:cases(v) for k,v in PHENO.items()}; copa=set().union(*domsets.values())
    S['cohortwide_cases']={k:len(v) for k,v in domsets.items()}|{'COPA_any':len(copa)}
    print("   cohort-wide EHR cases:",S['cohortwide_cases'])
    anc=pd.read_csv(ANC); anc['research_id']=anc.research_id.astype(str)
    d=g.merge(anc[['research_id','ancestry_pred']],on='research_id',how='left')
    d['R85H']=d.research_id.isin(r85h).astype(int); d['AQ']=(d.AQ_d>0).astype(int); d['HAQ']=(d.HAQ_d>0).astype(int)
    d['COPA']=d.research_id.isin(copa).astype(int)
    for k,s in domsets.items(): d[k]=d.research_id.isin(s).astype(int)
    S['n_cohort']=len(d); S['n_R85H']=int(d.R85H.sum()); S['n_COPA']=int(d.COPA.sum())
    print(f"   cohort {len(d)} | R85H {int(d.R85H.sum())} | COPA-spectrum {int(d.COPA.sum())} ({100*d.COPA.mean():.2f}%)")
    COVcols=[]
    try:
        if 'pca_features' in anc.columns:
            P=anc['pca_features'].apply(lambda x: ast.literal_eval(x) if isinstance(x,str) else (x if isinstance(x,list) else []))
            k=min(min(map(len,P[P.map(len)>0])) if (P.map(len)>0).any() else 0,16)
            pcs=pd.DataFrame([row[:k] for row in P],columns=[f'PC{i+1}' for i in range(k)]); pcs['research_id']=anc.research_id.values
            d=d.merge(pcs,on='research_id',how='left'); COVcols+=[f'PC{i+1}' for i in range(k)]
    except Exception as ex: print("   (PCs skipped:",str(ex)[:40],")")
    try:
        ppl=bq.query(f"SELECT person_id, DATE_DIFF(DATE '2025-01-01', DATE(birth_datetime), YEAR) age, CAST(sex_at_birth_concept_id AS STRING) sexc FROM `{PROJ}.{DS}.person`").to_dataframe()
        ppl['research_id']=ppl.person_id.astype(str); d=d.merge(ppl[['research_id','age','sexc']],on='research_id',how='left'); COVcols+=['age','C(sexc)']
    except Exception as ex: print("   (age/sex skipped:",str(ex)[:40],")")
    COV=(" + "+" + ".join(COVcols)) if COVcols else ""; print(f"   covariates: {COVcols}")
    import statsmodels.formula.api as smf
    def logit(f,dd,term):
        try:
            r=smf.logit(f,data=dd,missing='drop').fit(disp=0)
            return {'OR':round(float(np.exp(r.params[term])),3),'p':round(float(r.pvalues[term]),4),'n':int(r.nobs)}
        except Exception as e: return {'error':str(e)[:70]}
    print("== COPA-spectrum prevalence by R85H × STING class (both hits required?) ==")
    def prev(dd):
        dd=dd.copy(); dd['sc']=dd.apply(sc_of,axis=1)
        return {c:{'n':int((dd.sc==c).sum()),'cases':int(dd[dd.sc==c].COPA.sum()),'copa_pct':round(float(dd[dd.sc==c].COPA.mean()*100),2) if (dd.sc==c).any() else None} for c in ['WT/other','HAQ','AQ']}
    S['prev_R85Hneg']=prev(d[d.R85H==0]); S['prev_R85Hpos']=prev(d[d.R85H==1])
    print("   R85H- :",S['prev_R85Hneg']); print("   R85H+ :",S['prev_R85Hpos'])
    print("== R85H × STING -> phenotype (logistic, PC/age/sex-adj) ==")
    for out in ['COPA','ILD','arthritis','vasculitis']:
        S[f'{out}_R85HxAQ']=logit(f'{out} ~ R85H*AQ + C(ancestry_pred){COV}',d,'R85H:AQ')
        S[f'{out}_R85HxHAQ']=logit(f'{out} ~ R85H*HAQ + C(ancestry_pred){COV}',d,'R85H:HAQ')
        print(f"   {out}: R85H×AQ {S[f'{out}_R85HxAQ']} | R85H×HAQ {S[f'{out}_R85HxHAQ']}")
    print("== within R85H carriers (n≈138): AQ -> phenotype ==")
    rp=d[d.R85H==1]
    for out in ['COPA','ILD','arthritis','vasculitis']:
        S[f'within_AQ_{out}']=logit(f'{out} ~ AQ + C(ancestry_pred)',rp,'AQ')
        print(f"   R85H+ {out} ~ AQ: {S[f'within_AQ_{out}']}")
    print("\n===== PHENOME R85H×STING (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except FileNotFoundError:
    print("run failed: ~/sting_phenome_pav_v9.csv not found — run sting_only_extract_v9.py with PHENOME=1 first")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
