#!/usr/bin/env python
"""AoU v9 -- R85H as a SEVERITY MODIFIER across the FULL COPA/R85H SPECTRUM (generalizes the bronchiectasis-only runs
038/039, where R85H was incidence-neutral but ~3x hard-severe). Test the SAME logic on every COPA-syndrome phenotype:
  LUNG   : bronchiectasis, fibrotic-ILD, LIP (lymphoid interstitial pneumonia -- COPA-characteristic), alveolar-hemorrhage
  SYSTEMIC (autoimmune -- OA37 = RA + ANCA): inflammatory arthritis, ANCA-vasculitis (GPA/MPA), broad vasculitis, glomerulonephritis
For EACH phenotype: (i) R85H INCIDENCE (expect ~null -- incidence-neutral); (ii) R85H WITHIN-CASE SEVERITY -> organ-failure
HARD outcome (resp failure / O2 / dialysis-ESRD / transplant / mech-vent) = the generalizable severity-modifier test;
(iii) STING second-hit route: AQ / HAQ incidence per phenotype. Covariate-adjusted; FDR across phenotypes on the severity
column. Q: is R85H a severity amplifier of the WHOLE spectrum, or bronchiectasis-specific? DISCOVERY. Standard app.
Ends 'run complete' / 'run failed'.
"""
import os, gzip, ast, json
import numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); SNP=f"{MNT}/v9/wgs/short_read/snpindel"
ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
R85H_VID='19-18911007-C-T'; STINGV={'R71H':'5-139481493-C-T','G230A':'5-139478340-C-G','R293Q':'5-139477397-C-T'}
PANEL={
 'bronchiectasis':['J47','494'],
 'fibrotic_ILD':['J84.1','515','516.3'],
 'LIP':['J84.2'],                                     # lymphoid interstitial pneumonia (COPA-characteristic)
 'alveolar_hemorrhage':['R04.2','R04.89'],            # DAH proxy (hemoptysis/pulmonary hemorrhage) -- COPA-classic; noisy
 'infl_arthritis':['M05','M06','M08','714'],          # RA/JIA (OA37)
 'ANCA_vasculitis':['M31.3','M31.7'],                 # GPA/MPA (OA37 ANCA+)
 'vasculitis_broad':['M30','M31'],
 'glomerulonephritis':['N01','N03','N05','580','581','582','583'],
}
LUNGP=['bronchiectasis','fibrotic_ILD','LIP','alveolar_hemorrhage']
# organ-failure HARD outcome (uniform severity proxy across the spectrum)
HARD=['J96','518.81','518.83','518.84','Z99.81','V46.2','Z99.2','N18.6','Z94','V42.6','V42.0','Z99.11']
MINCODES=int(os.environ.get("MINCODES","2"))
def fnum(x):
    try: return float(x)
    except (TypeError,ValueError): return 0.0
def f2(x): return f"{x:.2f}" if isinstance(x,(int,float)) else "na"
def f3(x): return f"{x:.3f}" if isinstance(x,(int,float)) else "na"
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
    print("== phenotype panel: cases (>=%d codes) + hard-severe ==" % MINCODES)
    hard=anyset(HARD); cases={ph:caseset(cc) for ph,cc in PANEL.items()}
    for ph in PANEL: print(f"   {ph:20s} cases {len(cases[ph]):>6}  hard-severe {len(cases[ph]&hard):>5}")
    S['n_cases']={ph:len(cases[ph]) for ph in PANEL}; S['n_hard']=len(hard)
    # ---- exposures ----
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
    C_r85h=carr([R85H_VID]); C_irak3=carr(irak3); st={k:carr([v]) for k,v in STINGV.items()}
    C_aq=(st['G230A']&st['R293Q'])-st['R71H']; C_haq=st['R71H']&st['G230A']&st['R293Q']
    print(f"   R85H {len(C_r85h)} | AQ(unphased) {len(C_aq)} | HAQ(unphased) {len(C_haq)} | IRAK3-LoF {len(C_irak3)}")
    # ---- cohort + covariates ----
    anc=pd.read_csv(ANC,sep="\t"); anc['research_id']=anc.research_id.astype(str)
    d=anc[['research_id','ancestry_pred']].copy()
    d['R85H']=d.research_id.isin(C_r85h).astype(int); d['AQ']=d.research_id.isin(C_aq).astype(int)
    d['HAQ']=d.research_id.isin(C_haq).astype(int); d['IRAK3']=d.research_id.isin(C_irak3).astype(int)
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
    covfull=" + age + C(sexc)"+((" + "+" + ".join(PCS)) if PCS else "")
    covlite=" + age + C(sexc)"+((" + "+" + ".join(PCS[:5])) if PCS else "")
    import statsmodels.formula.api as smf
    def logit(fm,dd,t):
        try: r=smf.logit(fm,data=dd,missing='drop').fit(disp=0); return round(float(np.exp(r.params[t])),3),float(r.pvalues[t])
        except Exception: return None,None
    # ---- per-phenotype: R85H incidence, R85H within-case severity, AQ/HAQ incidence ----
    rows=[]
    for ph in PANEL:
        d['Y']=d.research_id.isin(cases[ph]).astype(int)
        inc_or,inc_p=logit(f'Y ~ R85H{covfull}',d,'R85H')
        aq_or,aq_p=logit(f'Y ~ AQ{covfull}',d,'AQ'); haq_or,haq_p=logit(f'Y ~ HAQ{covfull}',d,'HAQ')
        sub=d[d.Y==1].copy(); sub['SEV']=sub.research_id.isin(hard).astype(int)
        ncar=int(sub.R85H.sum()); nsev=int((sub.R85H*sub.SEV).sum())
        sev_or,sev_p=logit(f'SEV ~ R85H{covlite}',sub,'R85H') if ncar>=3 else (None,None)
        rows.append({'phenotype':ph,'n_cases':len(cases[ph]),'R85H_inc_OR':inc_or,'R85H_inc_p':inc_p,
                     'R85H_case_carr':ncar,'R85H_severe':nsev,'R85H_sev_OR':sev_or,'R85H_sev_p':sev_p,
                     'AQ_inc_OR':aq_or,'AQ_inc_p':aq_p,'HAQ_inc_OR':haq_or,'HAQ_inc_p':haq_p})
    sp=[r['R85H_sev_p'] for r in rows if r['R85H_sev_p'] is not None]
    fdr=dict(zip([r['phenotype'] for r in rows if r['R85H_sev_p'] is not None],multipletests(sp,method='fdr_bh')[1])) if sp else {}
    for r in rows: r['R85H_sev_fdr']=round(float(fdr.get(r['phenotype'],np.nan)),4) if r['phenotype'] in fdr else None
    S['spectrum']=rows
    print("\n== R85H across the COPA/R85H spectrum: INCIDENCE (expect ~null) vs WITHIN-CASE SEVERITY (the modifier test) ==")
    print(f"{'phenotype':20s} {'ncase':>6} | {'R85H_inc OR(p)':>16} | {'R85H_SEV OR(p) car/sev':>26} {'FDR':>6} | {'AQ_inc':>10} {'HAQ_inc':>10}")
    for r in rows:
        inc=f"{f3(r['R85H_inc_OR'])}({f2(r['R85H_inc_p'])})"
        sev=(f"{f3(r['R85H_sev_OR'])}({f3(r['R85H_sev_p'])}) {r['R85H_case_carr']}/{r['R85H_severe']}") if r['R85H_sev_OR'] is not None else f"(n={r['R85H_case_carr']})"
        aq=f"{f3(r['AQ_inc_OR'])}({f2(r['AQ_inc_p'])})"; haq=f"{f3(r['HAQ_inc_OR'])}({f2(r['HAQ_inc_p'])})"
        print(f"{r['phenotype']:20s} {r['n_cases']:>6} | {inc:>16} | {sev:>26} {f3(r['R85H_sev_fdr']):>6} | {aq:>10} {haq:>10}")
    print("\n   LUNG:",LUNGP," | SYSTEMIC: infl_arthritis, ANCA_vasculitis, vasculitis_broad, glomerulonephritis")
    print("   READ: R85H incidence ~1.0 everywhere but severity OR>1 across the spectrum => R85H = SPECTRUM-WIDE severity amplifier (not bronchiectasis-specific).")
    print("\n===== SPECTRUM SEVERITY (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
