#!/usr/bin/env python
"""AoU v9 — PHENOTYPE-FIRST reverse burden: start from parenchymal ILD, ask WHICH genes are enriched for rare damaging
variants. Tightened ILD cohort (J84/515/516, requiring >=2 occurrences = chronic, to cut the noisy single-code cases).
Curated 24-gene panel: ILD POSITIVE-CONTROLS (TERT/RTEL1/MUC5B/PARN/ABCA3/SFTPC) + interferonopathy (TREX1/IFIH1/
SAMHD1/DNASE2) + our COPI/STING/IRAK pathway. Per-gene rare-damaging burden (VAT LoF-HC or missense REVEL>0.5, afr_af
<0.1%) -> logistic ILD ~ burden + 16 PC + age + sex on the full ~535k cohort; FDR across the panel.
VALIDATION: do the known ILD genes light up (OR>1)? POSITIONING: does our pathway? DISCOVERY. Ends 'run complete'/'run failed'.
"""
import os, gzip, ast, json
import numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); SNP=f"{MNT}/v9/wgs/short_read/snpindel"
ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
GENES={  # GRCh38, Ensembl-verified
 'TERT':('5',1253147,1295086),'RTEL1':('20',63657810,63696253),'MUC5B':('11',1223066,1262172),'PARN':('16',14435699,14632728),'ABCA3':('16',2275878,2340788),'SFTPC':('8',22156913,22164479),
 'TREX1':('3',48465811,48467645),'IFIH1':('2',162267074,162318684),'SAMHD1':('20',36890216,36951893),'DNASE2':('19',12875209,12881595),
 'COPE':('19',18899511,18919407),'COPA':('1',160288580,160343566),'COPB1':('11',14443357,14500027),'COPB2':('3',139353942,139389736),'COPG1':('3',129249575,129278068),'COPG2':('7',130505553,130668755),'COPZ1':('12',54301202,54351848),'COPZ2':('17',48026142,48038030),'ARCN1':('11',118572384,118603033),
 'STING1':('5',139475528,139482935),'IRAK3':('12',66183995,66259622),'IRAK1':('X',154005501,154025650),'IRAK4':('12',43753938,43803307),'MYD88':('3',38133552,38148024)}
CAT={**{g:'ILD_posctrl' for g in ['TERT','RTEL1','MUC5B','PARN','ABCA3','SFTPC']},**{g:'interferonopathy' for g in ['TREX1','IFIH1','SAMHD1','DNASE2']}}
for g in GENES:
    CAT.setdefault(g,'pathway')
ILD_CODES=['J84','515','516']
S={}
def fnum(x):
    try: return float(x)
    except (TypeError,ValueError): return 0.0
try:
    import pysam
    with gzip.open(VAT,'rt') as fh: COLS=fh.readline().rstrip("\n").split("\t")
    IX={c:COLS.index(c) for c in ['vid','gene_symbol','consequence','gvs_afr_af','revel','LoF']}
    tbx=pysam.TabixFile(VAT); gv={}
    for g,(ch,s,e) in GENES.items():
        vids=[]
        for line in tbx.fetch('chr'+ch,s,e):
            f=line.split("\t")
            if f[IX['gene_symbol']]!=g: continue
            if fnum(f[IX['gvs_afr_af']])>=0.001: continue
            if f[IX['LoF']]=='HC' or ('missense' in f[IX['consequence']] and fnum(f[IX['revel']])>0.5): vids.append(f[IX['vid']])
        gv[g]=vids
    print("damaging vids/gene:",{g:len(v) for g,v in gv.items()})
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    def carr(vids):
        out=set()
        for i in range(0,len(vids),900):
            if not vids[i:i+900]: continue
            inl=",".join("'"+v+"'" for v in vids[i:i+900])
            out|=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid IN ({inl})"))
        return out
    lk=" OR ".join([f"c.concept_code LIKE '{p}%'" for p in ILD_CODES])
    ild=set(str(r.person_id) for r in bq.query(f"SELECT co.person_id FROM `{PROJ}.{DS}.condition_occurrence` co JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id WHERE c.vocabulary_id LIKE 'ICD%' AND ({lk}) GROUP BY co.person_id HAVING COUNT(*)>=2"))
    S['ILD_cases_chronic']=len(ild); print(f"parenchymal-ILD cases (>=2 codes): {len(ild)}")
    anc=pd.read_csv(ANC,sep="\t"); anc['research_id']=anc.research_id.astype(str)
    d=anc[['research_id','ancestry_pred']].copy(); d['ILD']=d.research_id.isin(ild).astype(int)
    COV=[]
    try:
        P=anc['pca_features'].apply(lambda x: ast.literal_eval(x) if isinstance(x,str) else (x if isinstance(x,list) else []))
        k=min(min(map(len,P[P.map(len)>0])) if (P.map(len)>0).any() else 0,16)
        p=pd.DataFrame([row[:k] for row in P],columns=[f'PC{i+1}' for i in range(k)]); p['research_id']=anc.research_id.values; d=d.merge(p,on='research_id',how='left'); COV+=[f'PC{i+1}' for i in range(k)]
    except Exception: pass
    try:
        ppl=bq.query(f"SELECT person_id, DATE_DIFF(DATE '2025-01-01', DATE(birth_datetime), YEAR) age, CAST(sex_at_birth_concept_id AS STRING) sexc FROM `{PROJ}.{DS}.person`").to_dataframe()
        ppl['research_id']=ppl.person_id.astype(str); d=d.merge(ppl[['research_id','age','sexc']],on='research_id',how='left'); COV+=['age','C(sexc)']
    except Exception: pass
    covs=(" + "+" + ".join(COV)) if COV else ""; print(f"cohort {len(d)} | ILD prevalence {100*d.ILD.mean():.2f}% | covariates {len(COV)}")
    import statsmodels.formula.api as smf
    def logit(f,dd,t):
        try: r=smf.logit(f,data=dd,missing='drop').fit(disp=0); return {'OR':round(float(np.exp(r.params[t])),3),'p':float(r.pvalues[t])}
        except Exception as e: return {'OR':None,'p':None,'err':str(e)[:40]}
    res=[]
    for g in GENES:
        d[g]=d.research_id.isin(carr(gv[g])).astype(int); ncar=int(d[g].sum()); ncc=int((d[g]*d.ILD).sum())
        r=logit(f'ILD ~ {g}{covs}',d,g); res.append({'gene':g,'cat':CAT[g],'carriers':ncar,'ILD_carriers':ncc,'OR':r.get('OR'),'p':r.get('p')})
    ps=[x['p'] for x in res if x['p'] is not None]
    fdr=dict(zip([x['gene'] for x in res if x['p'] is not None],multipletests(ps,method='fdr_bh')[1])) if ps else {}
    for x in res: x['fdr']=round(float(fdr.get(x['gene'],np.nan)),4) if x['gene'] in fdr else None; x['p']=round(x['p'],5) if x['p'] is not None else None
    res=sorted(res,key=lambda x:(x['p'] if x['p'] is not None else 9))
    S['results']=res
    print("\n== reverse burden: ILD ~ gene rare-damaging (sorted by p) ==")
    print(f"{'gene':8s} {'cat':16s} {'carr':>6s} {'ILDcarr':>7s} {'OR':>6s} {'p':>9s} {'FDR':>7s}")
    for x in res: print(f"{x['gene']:8s} {x['cat']:16s} {x['carriers']:>6d} {x['ILD_carriers']:>7d} {str(x['OR']):>6s} {str(x['p']):>9s} {str(x['fdr']):>7s}")
    pc=[x for x in res if x['cat']=='ILD_posctrl']
    print(f"\n== VALIDATION: ILD positive-controls with OR>1: {sum(1 for x in pc if x['OR'] and x['OR']>1)}/{len(pc)} ({[x['gene'] for x in pc if x['OR'] and x['OR']>1 and x['p'] and x['p']<0.05]} sig) ==")
    print("\n===== ILD REVERSE BURDEN (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
