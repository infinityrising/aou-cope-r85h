#!/usr/bin/env python
"""AoU v9 -- the run-039 confound battery GENERALIZED to OA37's TRIAD. Run 040 found R85H amplifies within-case SEVERITY
on bronchiectasis (OR 2.77), inflammatory arthritis (1.72) and vasculitis (11.2, n=6) but NOT fibrotic-ILD / renal.
Bronchiectasis survived ancestry+utilization (run 039); the arthritis + vasculitis hits have NOT been controlled. Apply
the identical battery to each phenotype in the triad:
 (1) ANCESTRY -- WITHIN-AFR 2x2 (R85H+ vs R85H- severe, ancestry held constant) + 16-PC adjustment.
 (2) CODING INTENSITY -- total condition-code count (utilization) covariate.
 HARD severity = organ failure (resp failure / O2-dep / dialysis-ESRD / transplant / mech-vent). Among CASES of each
 phenotype, does R85H->SEVERE survive? bronchiectasis = built-in positive control (must reproduce run 039 ~3x).
 PRE-SPECIFIED (locked): R85H->severe survives ancestry+utilization for OA37's triad. DISCOVERY. Standard app.
 Ends 'run complete' / 'run failed'.
"""
import os, gzip, ast, json
import numpy as np, pandas as pd
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); SNP=f"{MNT}/v9/wgs/short_read/snpindel"
ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
R85H_VID='19-18911007-C-T'
PANEL={
 'bronchiectasis':['J47','494'],          # positive control (survived run 039)
 'infl_arthritis':['M05','M06','M08','714'],
 'vasculitis_broad':['M30','M31'],
 'ANCA_vasculitis':['M31.3','M31.7'],      # OA37's exact Dx (small)
}
HARD=['J96','518.81','518.83','518.84','Z99.81','V46.2','Z99.2','N18.6','Z94','V42.6','V42.0','Z99.11']
MINCODES=int(os.environ.get("MINCODES","2"))
def fnum(x):
    try: return float(x)
    except (TypeError,ValueError): return 0.0
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
    hard=anyset(HARD); cases={ph:caseset(cc) for ph,cc in PANEL.items()}
    def carr(vids):
        out=set()
        for i in range(0,len(vids),900):
            if not vids[i:i+900]: continue
            inl=",".join("'"+v+"'" for v in vids[i:i+900])
            out|=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid IN ({inl})"))
        return out
    C_r85h=carr([R85H_VID])
    # utilization (total condition codes/person) for the union of all triad cases
    allc=set().union(*cases.values()); util={}; bl=[b for b in allc]
    for i in range(0,len(bl),5000):
        inl=",".join(bl[i:i+5000])
        if not inl: continue
        for r in bq.query(f"SELECT person_id, COUNT(*) n FROM `{PROJ}.{DS}.condition_occurrence` WHERE person_id IN ({inl}) GROUP BY person_id"): util[str(r.person_id)]=r.n
    # cohort frame
    anc=pd.read_csv(ANC,sep="\t"); anc['research_id']=anc.research_id.astype(str)
    d=anc[['research_id','ancestry_pred']].copy()
    d['R85H']=d.research_id.isin(C_r85h).astype(int); d['hard']=d.research_id.isin(hard).astype(int)
    d['log_util']=np.log1p(d.research_id.map(lambda r: util.get(r,0)))
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
    import statsmodels.formula.api as smf
    def logit(fm,dd,t):
        try: r=smf.logit(fm,data=dd,missing='drop').fit(disp=0); return round(float(np.exp(r.params[t])),3),float(r.pvalues[t])
        except Exception: return None,None
    pc16=(" + "+" + ".join(PCS)) if PCS else ""; pc5=(" + "+" + ".join(PCS[:5])) if PCS else ""
    print(f"cohort {len(d)} | PCs {len(PCS)} | R85H {int(d.R85H.sum())} | hard-outcome {int(d.hard.sum())}")
    print("\n== R85H -> SEVERE, confound-controlled, across OA37's triad (bronchiectasis = run-039 positive control) ==")
    for ph in PANEL:
        sub=d[d.research_id.isin(cases[ph])].copy(); sub['SEV']=sub.research_id.isin(hard).astype(int)
        ncar=int(sub.R85H.sum()); nsev=int((sub.R85H*sub.SEV).sum())
        or16,p16=logit(f'SEV ~ R85H + age + C(sexc){pc16}',sub,'R85H') if ncar>=3 else (None,None)
        oru,pu=logit(f'SEV ~ R85H + age + C(sexc) + log_util{pc16}',sub,'R85H') if ncar>=3 else (None,None)
        afr=sub[sub.ancestry_pred=='afr']; pos=afr[afr.R85H==1]; neg=afr[afr.R85H==0]
        fp=float(pos.SEV.mean()) if len(pos) else float('nan'); fn=float(neg.SEV.mean()) if len(neg) else float('nan')
        ora,pa=logit(f'SEV ~ R85H + age + C(sexc){pc5}',afr,'R85H') if int(pos.R85H.sum() if len(pos) else 0)>=3 else (None,None)
        S[ph]={'cases':len(cases[ph]),'R85H_cases':ncar,'R85H_severe':nsev,
               'OR_16PC':or16,'p_16PC':p16,'OR_util':oru,'p_util':pu,
               'AFR_R85Hpos':f"{int(pos.SEV.sum())}/{len(pos)}" if len(pos) else "0/0",'AFR_R85Hpos_frac':round(fp,3) if len(pos) else None,'AFR_R85Hneg_frac':round(fn,3) if len(neg) else None,'AFR_adjOR':ora,'AFR_p':pa}
        pf=f"{100*fp:.0f}%" if len(pos) else "n/a"; nf=f"{100*fn:.0f}%" if len(neg) else "n/a"
        print(f"\n  [{ph}]  cases {len(cases[ph])} | R85H-cases {ncar} severe {nsev}")
        print(f"     all-cases +16PC   OR={or16} p={p16}")
        print(f"     all-cases +util   OR={oru} p={pu}")
        print(f"     WITHIN-AFR        R85H+ {S[ph]['AFR_R85Hpos']} ({pf}) severe vs R85H- {nf} | adj OR={ora} p={pa}")
    print("\n== READ: does R85H->severe survive within-AFR + 16PC + utilization for arthritis & vasculitis (as it did for bronchiectasis)? ==")
    print("\n===== SEVERITY SPECTRUM CONTROL (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
