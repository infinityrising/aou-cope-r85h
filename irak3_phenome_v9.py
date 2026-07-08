#!/usr/bin/env python
"""AoU v9 — IRAK3 (IRAK-M) phenome scan on the FULL WGS cohort. IRAK3 has almost NO human disease literature (only
asthma: Balaci/OMIM 611064; LoF-tolerant; healthy hom-nonsense adults exist). ★ NOVEL QUESTION: is IRAK3-LoF tied to
ILD / CYSTIC LUNG DISEASE? = OA37's phenotype (R85H + IRAK3-L210*, cystic lung dz). Powered: ~2,500 IRAK3-damaging
carriers, full ~535k EHR cohort (no phasing, not RNA-gated). Tests (logistic, 16 PC + age + sex adj):
  (1) IRAK3 -> asthma        = POSITIVE CONTROL (established RISK; validates the calls + pathway).
  (2) IRAK3 -> ILD / cystic_ILD = NOVEL discovery, direction OPEN (mouse: IRAK-M KO PROTECTS bleomycin fibrosis but
      WORSENS airway inflammation -> bidirectional; report honestly whichever way).
  (3) IRAK3 -> infl_arthritis / vasculitis (COPA-spectrum).
  (4) R85H x IRAK3 -> lung disease (the two-hit; double-carriers rare -> expect underpowered, report cell counts).
DISCOVERY / exploratory (a positive novel hit -> pre-register a replication). STANDARD app. Ends 'run complete'/'run failed'.
"""
import os, gzip, ast, json
import numpy as np, pandas as pd
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); SNP=f"{MNT}/v9/wgs/short_read/snpindel"
ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1); R85H_VID='19-18911007-C-T'
IRAK3=('12',66183995,66259622)
PHE={'asthma':['J45%','493%'],
     'ILD':['J84%','515%','516%','J98.4'],
     'cystic_ILD':['J84%','515%','516%','J98.4','J47%','494%'],   # ILD + bronchiectasis (cystic airway) = OA37-style cystic/ILD
     'infl_arthritis':['M05%','M06%','M07%','M08%','714%'],
     'vasculitis':['M31.3%','M31.7%','M30%','M31%','446%','I77.6']}
S={}
def fnum(x):
    try: return float(x)
    except (TypeError,ValueError): return 0.0
try:
    print("== IRAK3 damaging vids (VAT: rare LoF-HC or missense REVEL>0.5) ==")
    import pysam
    with gzip.open(VAT,'rt') as fh: COLS=fh.readline().rstrip("\n").split("\t")
    IX={c:COLS.index(c) for c in ['vid','gene_symbol','consequence','gvs_afr_af','revel','LoF']}
    tbx=pysam.TabixFile(VAT); lof=[]; dmg=[]
    for line in tbx.fetch('chr'+IRAK3[0],IRAK3[1],IRAK3[2]):
        f=line.split("\t")
        if f[IX['gene_symbol']]!='IRAK3': continue
        if fnum(f[IX['gvs_afr_af']])>=0.001: continue
        if f[IX['LoF']]=='HC': lof.append(f[IX['vid']]); dmg.append(f[IX['vid']])
        elif 'missense' in f[IX['consequence']] and fnum(f[IX['revel']])>0.5: dmg.append(f[IX['vid']])
    S['IRAK3_LoF_vids']=len(lof); S['IRAK3_dmg_vids']=len(dmg); print(f"   IRAK3 LoF-HC {len(lof)} | damaging {len(dmg)}")
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    def carriers(vids):
        out=set()
        for i in range(0,len(vids),900):
            if not vids[i:i+900]: continue
            inl=",".join("'"+v+"'" for v in vids[i:i+900])
            out|=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid IN ({inl})"))
        return out
    irak3_lof=carriers(lof); irak3_dmg=carriers(dmg); r85h=carriers([R85H_VID])
    print(f"   carriers: IRAK3-LoF {len(irak3_lof)} | IRAK3-damaging {len(irak3_dmg)} | R85H {len(r85h)}")
    def cases(codes):
        lk=" OR ".join([f"c.concept_code LIKE '{p}'" for p in codes])
        q=f"""SELECT DISTINCT co.person_id FROM `{PROJ}.{DS}.condition_occurrence` co JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id WHERE c.vocabulary_id LIKE 'ICD%' AND ({lk})"""
        return set(str(r.person_id) for r in bq.query(q))
    caseset={k:cases(v) for k,v in PHE.items()}; S['cohortwide_cases']={k:len(v) for k,v in caseset.items()}
    smoke=cases(['F17%','Z72.0%','305.1%','V15.82'])   # ever-smoker proxy (EHR) — the lung-specific confound to rule out
    S['smokers']=len(smoke); print("   cohort-wide cases:",S['cohortwide_cases'],"| smokers:",len(smoke))
    print("== cohort (WGS w/ ancestry+PCs) + covariates ==")
    anc=pd.read_csv(ANC,sep="\t"); anc['research_id']=anc.research_id.astype(str)
    d=anc[['research_id','ancestry_pred']].copy()
    d['IRAK3']=d.research_id.isin(irak3_dmg).astype(int); d['IRAK3_lof']=d.research_id.isin(irak3_lof).astype(int); d['R85H']=d.research_id.isin(r85h).astype(int)
    for k,s in caseset.items(): d[k]=d.research_id.isin(s).astype(int)
    d['smoking']=d.research_id.isin(smoke).astype(int)
    COV=[]
    try:
        P=anc['pca_features'].apply(lambda x: ast.literal_eval(x) if isinstance(x,str) else (x if isinstance(x,list) else []))
        k=min(min(map(len,P[P.map(len)>0])) if (P.map(len)>0).any() else 0,16)
        pcs=pd.DataFrame([row[:k] for row in P],columns=[f'PC{i+1}' for i in range(k)]); pcs['research_id']=anc.research_id.values
        d=d.merge(pcs,on='research_id',how='left'); COV+=[f'PC{i+1}' for i in range(k)]
    except Exception as ex: print("   (PCs skipped:",str(ex)[:40],")")
    try:
        ppl=bq.query(f"SELECT person_id, DATE_DIFF(DATE '2025-01-01', DATE(birth_datetime), YEAR) age, CAST(sex_at_birth_concept_id AS STRING) sexc FROM `{PROJ}.{DS}.person`").to_dataframe()
        ppl['research_id']=ppl.person_id.astype(str); d=d.merge(ppl[['research_id','age','sexc']],on='research_id',how='left'); COV+=['age','C(sexc)']
    except Exception as ex: print("   (age/sex skipped:",str(ex)[:40],")")
    covs=(" + "+" + ".join(COV)) if COV else ""; S['n_cohort']=len(d); print(f"   cohort {len(d)} | covariates {len(COV)}")
    import statsmodels.formula.api as smf
    def logit(f,dd,term):
        try:
            r=smf.logit(f,data=dd,missing='drop').fit(disp=0)
            return {'OR':round(float(np.exp(r.params[term])),3),'p':round(float(r.pvalues[term]),5),'n':int(r.nobs)}
        except Exception as e: return {'error':str(e)[:70]}
    print("== IRAK3 -> phenotypes (PC+age+sex adj, then +SMOKING). asthma=positive control; ILD/cystic=NOVEL ==")
    for ph in ['asthma','ILD','cystic_ILD','infl_arthritis','vasculitis']:
        ncar=int(d.IRAK3.sum()); ncc=int((d.IRAK3*d[ph]).sum())
        S[f'{ph}_IRAK3']=logit(f'{ph} ~ IRAK3{covs}',d,'IRAK3'); S[f'{ph}_IRAK3lof']=logit(f'{ph} ~ IRAK3_lof{covs}',d,'IRAK3_lof')
        S[f'{ph}_IRAK3_smk']=logit(f'{ph} ~ IRAK3{covs} + smoking',d,'IRAK3')   # smoking-adjusted (does the lung signal survive?)
        print(f"   {ph}: carriers {ncar} cases {ncc} ({100*ncc/max(ncar,1):.1f}%) -> dmg {S[f'{ph}_IRAK3']} | +smk {S[f'{ph}_IRAK3_smk']} | LoF {S[f'{ph}_IRAK3lof']}")
    print("== R85H × IRAK3 -> lung disease (two-hit; report double-carrier counts) ==")
    for ph in ['cystic_ILD','ILD']:
        dd=d.copy(); dd['inter']=dd.R85H*dd.IRAK3
        ndbl=int((dd.R85H*dd.IRAK3).sum()); ndc=int((dd.R85H*dd.IRAK3*dd[ph]).sum())
        S[f'2hit_{ph}']=logit(f'{ph} ~ R85H + IRAK3 + inter{covs}',dd,'inter')
        print(f"   R85H×IRAK3 -> {ph}: double-carriers {ndbl}, double-carrier cases {ndc} -> {S[f'2hit_{ph}']}")
    print("\n===== IRAK3 PHENOME (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
