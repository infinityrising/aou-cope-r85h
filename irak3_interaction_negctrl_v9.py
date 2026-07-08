#!/usr/bin/env python
"""AoU v9 -- IRAK3 co-equal package (completes the IRAK3 arm per the lit agent). Three tests:
 (1) CONVERGENCE/REDUNDANCY: IRAK3-LoF × STING-AQ / × HAQ / × R85H interactions on fibrotic-ILD + bronchiectasis +
     autoinflammatory -- if both second hits feed one axis, the double is sub-additive/redundant (interaction ~0 or <0),
     not synergistic.
 (2) DUAL SIGNATURE (H5c): IRAK3-LoF should be PROTECTIVE for fibrotic-ILD (Ballinger M1-skew) but ELEVATED for
     autoinflammatory/exacerbation phenotypes -- the M1 phenocopy.
 (3) LIPSITCH SPECIFICITY: negative-control OUTCOMES (refractive error / appendicitis / back pain / benign neoplasm --
     expect NULL) + negative-control GENOTYPE (synonymous-IRAK3 burden -> fibrotic-ILD, expect NULL -> confirms the signal
     is LoF-specific, not a general IRAK3-locus effect). Covariate-adjusted. Standard app. Ends 'run complete'/'run failed'.
"""
import os, gzip, ast, json
import numpy as np, pandas as pd
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); SNP=f"{MNT}/v9/wgs/short_read/snpindel"
ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
R85H_VID='19-18911007-C-T'; STINGV={'R71H':'5-139481493-C-T','G230A':'5-139478340-C-G','R293Q':'5-139477397-C-T'}
IRAK3=('12',66183995,66259622)
PHENO={'fibrotic_ILD':['J84.1','515','516.3'],'bronchiectasis':['J47','494'],
 'autoinflammatory':['M05','M06','M08','714','M30','M31','N01','N03','N05','L40','K50','K51'],  # RA/vasculitis/GN/psoriasis/IBD
 'asthma':['J45','493']}
NEGCTRL_OUT={'refractive_error':['H52','367'],'appendicitis':['K35','540'],'back_pain':['M54','724'],'benign_neoplasm':['D10','D12','D18','D36','216']}
MINCODES=2
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
    def carr(vids):
        out=set()
        for i in range(0,len(vids),900):
            if not vids[i:i+900]: continue
            inl=",".join("'"+v+"'" for v in vids[i:i+900])
            out|=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid IN ({inl})"))
        return out
    with gzip.open(VAT,'rt') as fh: COLS=fh.readline().rstrip("\n").split("\t")
    IX={c:COLS.index(c) for c in ['vid','gene_symbol','consequence','gvs_afr_af','LoF']}
    tbx=pysam.TabixFile(VAT); irak3_lof=[]; irak3_syn=[]
    for line in tbx.fetch('chr'+IRAK3[0],IRAK3[1],IRAK3[2]):
        f=line.split("\t")
        if f[IX['gene_symbol']]!='IRAK3' or fnum(f[IX['gvs_afr_af']])>=0.001: continue
        if f[IX['LoF']]=='HC': irak3_lof.append(f[IX['vid']])
        elif 'synonymous' in f[IX['consequence']]: irak3_syn.append(f[IX['vid']])
    C_lof=carr(irak3_lof); C_syn=carr(irak3_syn); C_r85h=carr([R85H_VID]); st={k:carr([v]) for k,v in STINGV.items()}
    C_aq=(st['G230A']&st['R293Q'])-st['R71H']; C_haq=st['R71H']&st['G230A']&st['R293Q']
    print(f"IRAK3-LoF {len(C_lof)} | IRAK3-synonymous {len(C_syn)} | AQ {len(C_aq)} | HAQ {len(C_haq)} | R85H {len(C_r85h)}")
    anc=pd.read_csv(ANC,sep="\t"); anc['research_id']=anc.research_id.astype(str)
    d=anc[['research_id','ancestry_pred']].copy()
    for nm,cs in [('IRAK3',C_lof),('IRAK3syn',C_syn),('AQ',C_aq),('HAQ',C_haq),('R85H',C_r85h)]: d[nm]=d.research_id.isin(cs).astype(int)
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
    except Exception: pass
    pc=(" + "+" + ".join(PCS)) if PCS else ""
    import statsmodels.formula.api as smf
    def logit(fm,dd,t):
        try: r=smf.logit(fm,data=dd,missing='drop').fit(disp=0); return round(float(np.exp(r.params[t])),3),round(float(r.pvalues[t]),4)
        except Exception: return None,None
    print("\n== (2) IRAK3-LoF INCIDENCE across phenotypes + Lipsitch NEG-CONTROL OUTCOMES + synonymous-IRAK3 neg-control genotype ==")
    S['incidence']={}
    for ph,cc in {**PHENO,**NEGCTRL_OUT}.items():
        d['Y']=d.research_id.isin(caseset(cc)).astype(int)
        lof=logit(f'Y ~ IRAK3 + age + C(sexc){pc}',d,'IRAK3'); syn=logit(f'Y ~ IRAK3syn + age + C(sexc){pc}',d,'IRAK3syn') if ph=='fibrotic_ILD' else (None,None)
        S['incidence'][ph]={'IRAK3_LoF':lof,'IRAK3_syn':syn}
        tag='  <NEG-CTRL' if ph in NEGCTRL_OUT else ''
        print(f"   {ph:18s} IRAK3-LoF OR={lof[0]}(p{lof[1]})"+(f" | synonymous OR={syn[0]}(p{syn[1]})" if syn[0] is not None else "")+tag)
    print("\n== (1) CONVERGENCE/REDUNDANCY: IRAK3-LoF × (AQ / HAQ / R85H) interaction (predict ~0/sub-additive if one axis) ==")
    S['interaction']={}
    for ph,cc in [('fibrotic_ILD',PHENO['fibrotic_ILD']),('bronchiectasis',PHENO['bronchiectasis']),('autoinflammatory',PHENO['autoinflammatory'])]:
        d['Y']=d.research_id.isin(caseset(cc)).astype(int); S['interaction'][ph]={}
        for other in ['AQ','HAQ','R85H']:
            o,p=logit(f'Y ~ IRAK3*{other} + age + C(sexc){pc}',d,f'IRAK3:{other}')
            nd=int(((d.IRAK3==1)&(d[other]==1)&(d.Y==1)).sum()); S['interaction'][ph][other]={'OR':o,'p':p,'double_cases':nd}
            print(f"   {ph:16s} IRAK3×{other:4s} interaction OR={o}(p{p}) | double-carrier cases={nd}")
    print("\n== READ: IRAK3-LoF protective for fibrotic-ILD but risk for autoinflammatory (dual M1 signature)? flat across neg-ctrl outcomes + synonymous (specificity)? interaction ~0 (convergence, not synergy)? ==")
    print("\n===== IRAK3 INTERACTION + NEG-CONTROL (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
