#!/usr/bin/env python
"""AoU v9 — SENTINEL EXCLUSION AUDIT. The preregistration requires the sentinel (AoU person_id 1348881 — the
COPA-phenocopy R85H carrier the analysis 'converges on') be EXCLUDED from association tests (circular ascertainment).
This checks whether 1348881 was INCLUDED in either test cohort — molecular (~/copi_sting_pav_v9.csv) or phenome
(~/sting_phenome_pav_v9.csv) — plus their STING haplotype and COPA-spectrum EHR status, so we can exclude + re-verify.
OA37 is external to AoU (not checkable here). STANDARD app. Ends 'run complete' / 'run failed'.
"""
import os, json
import pandas as pd
SENT='1348881'
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1); R85H_VID='19-18911007-C-T'
RNA=os.path.expanduser("~/copi_sting_pav_v9.csv"); PH=os.path.expanduser("~/sting_phenome_pav_v9.csv")
PHENO={'ILD':['J84%','515%','516%','J98.4','D86%'],'infl_arthritis':['M05%','M06%','M07%','M08%','714%'],'vasculitis':['M31.3%','M31.7%','M30%','M31%','446%','I77.6']}
S={'sentinel':SENT}
def hapinfo(path,label):
    try:
        g=pd.read_csv(path); g['research_id']=g.research_id.astype(str); row=g[g.research_id==SENT]
        if len(row):
            r=row.iloc[0]; S[label]={'present':True,'hap1':r.get('hap1'),'hap2':r.get('hap2'),'AQ_d':int(r.get('AQ_d',0)),'HAQ_d':int(r.get('HAQ_d',0)),'R85H_d':int(r.get('R85H_d',0)) if 'R85H_d' in g.columns else 'n/a'}
        else: S[label]={'present':False,'cohort_n':len(g)}
    except FileNotFoundError: S[label]='(csv not found)'
    print(f"   {label}: {S[label]}")
try:
    print("== is the sentinel in the test cohorts? ==")
    hapinfo(RNA,'in_molecular_RNA_PAV'); hapinfo(PH,'in_phenome_allPAV')
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    r85h=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid='{R85H_VID}'"))
    S['sentinel_is_R85H_carrier']=SENT in r85h; print(f"   R85H carrier (BQ): {SENT in r85h}")
    codes=[c for v in PHENO.values() for c in v]; lk=" OR ".join([f"c.concept_code LIKE '{p}'" for p in codes])
    q=f"""SELECT DISTINCT co.person_id FROM `{PROJ}.{DS}.condition_occurrence` co
          JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id
          WHERE c.vocabulary_id LIKE 'ICD%' AND ({lk}) AND co.person_id={SENT}"""
    is_copa=len([r for r in bq.query(q)])>0; S['sentinel_is_COPA_case']=is_copa
    print(f"   sentinel is a COPA-spectrum EHR case: {is_copa}")
    print("\n== VERDICT ==")
    inmol=isinstance(S.get('in_molecular_RNA_PAV'),dict) and S['in_molecular_RNA_PAV'].get('present')
    inph=isinstance(S.get('in_phenome_allPAV'),dict) and S['in_phenome_allPAV'].get('present')
    if inph and is_copa: print("   ⚠️ SENTINEL IS IN THE PHENOME COHORT AND IS A COPA CASE -> circular; EXCLUDE + re-run phenome (conservative).")
    elif inph: print("   sentinel in phenome cohort but not a COPA case -> low impact; exclude for compliance.")
    elif inmol: print("   sentinel in molecular cohort -> exclude for compliance (molecular arm is discovery anyway).")
    else: print("   sentinel NOT in either test cohort -> no exclusion needed; pre-reg exclusion satisfied de facto. Document it.")
    print("\n===== SENTINEL AUDIT (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-900:])
