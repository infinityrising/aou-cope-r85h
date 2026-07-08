#!/usr/bin/env python
"""AoU v9 -- SEVERITY axis, confound-controlled. Run 038 found R85H converts bronchiectasis mild->severe WITHIN cases
(15/19 severe, OR 3.72, p=0.02) while leaving incidence untouched (OR 1.00) -- i.e. R85H is a SEVERITY modifier, not an
incidence driver. That within-cases model was age+sex only; before believing it, close two confounds:
 (1) ANCESTRY. R85H is AFR-enriched (1.3%); severe-coding could tag the AFR background, not R85H. Controls: (a) add 16
     PCs to the within-cases severity model; (b) repeat WITHIN AFR bronchiectasis cases only (R85H+ vs R85H-, ancestry
     held constant -- the clean 2x2).
 (2) CODING INTENSITY. More encounters -> more incidental severity codes. Controls: adjust for total condition-code
     count (utilization), and use a HARD severity tier (respiratory failure / O2-dependence / lung transplant) that is
     encounter-robust, alongside the broad any-marker tier.
PRE-SPECIFIED direction (locked): R85H (and R85H+2nd-hit) -> SEVERE (OR>1) survives ancestry + utilization adjustment
and holds within AFR. DISCOVERY / exploratory. Standard app. Ends 'run complete' / 'run failed'.
"""
import os, gzip, ast, json
import numpy as np, pandas as pd
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); SNP=f"{MNT}/v9/wgs/short_read/snpindel"
ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
R85H_VID='19-18911007-C-T'; STINGV={'R71H':'5-139481493-C-T','G230A':'5-139478340-C-G','R293Q':'5-139477397-C-T'}
BRONCH=['J47','494']
SEV_ANY=['J96','518.81','518.83','518.84','Z99.81','V46.2','B96.5','R04.2','786.3','I27','416','Z94.2']
SEV_HARD=['J96','518.81','518.83','518.84','Z99.81','V46.2','Z94.2','V42.6']   # resp failure / O2-dep / transplant (encounter-robust)
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
    bronch=caseset(BRONCH); sev_any=bronch & anyset(SEV_ANY); sev_hard=bronch & anyset(SEV_HARD)
    S['n_bronchiect']=len(bronch); S['n_severe_any']=len(sev_any); S['n_severe_hard']=len(sev_hard)
    print(f"bronchiectasis {len(bronch)} | severe_any {len(sev_any)} ({100*len(sev_any)/max(len(bronch),1):.0f}%) | severe_hard {len(sev_hard)} ({100*len(sev_hard)/max(len(bronch),1):.0f}%)")
    # ---- exposures ----
    def carr(vids):
        out=set()
        for i in range(0,len(vids),900):
            if not vids[i:i+900]: continue
            inl=",".join("'"+v+"'" for v in vids[i:i+900])
            out|=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid IN ({inl})"))
        return out
    with gzip.open(VAT,'rt') as fh: COLS=fh.readline().rstrip("\n").split("\t")
    IX={c:COLS.index(c) for c in ['vid','gene_symbol','consequence','gvs_afr_af','LoF']}
    tbx=pysam.TabixFile(VAT); irak3=[]
    for line in tbx.fetch('chr12',66183995,66259622):
        f=line.split("\t")
        if f[IX['gene_symbol']]=='IRAK3' and f[IX['LoF']]=='HC' and fnum(f[IX['gvs_afr_af']])<0.001: irak3.append(f[IX['vid']])
    C_r85h=carr([R85H_VID]); C_irak3=carr(irak3); stcar={k:carr([v]) for k,v in STINGV.items()}
    C_aq=(stcar['G230A']&stcar['R293Q'])-stcar['R71H']
    # ---- utilization (total condition codes / person) among cases ----
    util={}; bl=[b for b in bronch]
    for i in range(0,len(bl),5000):
        inl=",".join(bl[i:i+5000])
        if not inl: continue
        for r in bq.query(f"SELECT person_id, COUNT(*) n FROM `{PROJ}.{DS}.condition_occurrence` WHERE person_id IN ({inl}) GROUP BY person_id"): util[str(r.person_id)]=r.n
    # ---- case-level frame ----
    anc=pd.read_csv(ANC,sep="\t"); anc['research_id']=anc.research_id.astype(str)
    d=anc[anc.research_id.isin(bronch)][['research_id','ancestry_pred']].copy()
    d['sev_any']=d.research_id.isin(sev_any).astype(int); d['sev_hard']=d.research_id.isin(sev_hard).astype(int)
    d['R85H']=d.research_id.isin(C_r85h).astype(int)
    d['R85H_x_AQ']=d.research_id.isin(C_r85h&C_aq).astype(int); d['IRAK3_LoF']=d.research_id.isin(C_irak3).astype(int)
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
    print(f"case frame {len(d)} | PCs {len(PCS)} | R85H cases {int(d.R85H.sum())} | AFR cases {int((d.ancestry_pred=='afr').sum())}")
    import statsmodels.formula.api as smf
    def logit(fm,dd,t):
        try:
            r=smf.logit(fm,data=dd,missing='drop').fit(disp=0); return round(float(np.exp(r.params[t])),3),float(r.pvalues[t]),int(r.nobs)
        except Exception: return None,None,None
    # ---- (1) WITHIN-AFR: R85H+ vs R85H- severe fraction, ancestry held constant ----
    print("\n== (1) WITHIN-AFR bronchiectasis cases: R85H+ vs R85H- severe (the clean 2x2) ==")
    afr=d[d.ancestry_pred=='afr'].copy(); afrpc=" + ".join(PCS[:5]); afrpc=(" + "+afrpc) if afrpc else ""
    S['within_afr']={}
    for tier in ['sev_any','sev_hard']:
        pos=afr[afr.R85H==1]; neg=afr[afr.R85H==0]
        fp=float(pos[tier].mean()) if len(pos) else float('nan'); fn=float(neg[tier].mean()) if len(neg) else float('nan')
        OR,p,n=logit(f'{tier} ~ R85H + age + C(sexc){afrpc}',afr,'R85H')
        S['within_afr'][tier]={'R85Hpos':f"{int(pos[tier].sum())}/{len(pos)}",'R85Hpos_frac':round(fp,3),'R85Hneg_frac':round(fn,3),'adjOR':OR,'p':p}
        pf=f"{100*fp:.0f}%" if len(pos) else "n/a"; nf=f"{100*fn:.0f}%" if len(neg) else "n/a"
        print(f"   [{tier:8s}] AFR R85H+ {int(pos[tier].sum())}/{len(pos)} ({pf}) severe  vs  R85H- {nf}  | adj OR={OR} p={p}")
    # ---- (2) ALL cases: ancestry(16PC)- and utilization-adjusted ----
    print("\n== (2) ALL bronchiectasis cases: severe ~ exposure + 16PC + age + sex  (+utilization) ==")
    pcterm=(" + "+" + ".join(PCS)) if PCS else ""; S['adjusted']=[]
    for tier in ['sev_any','sev_hard']:
        for ex in ['R85H','R85H_x_AQ','IRAK3_LoF']:
            n_ex=int(d[ex].sum()); n_sev=int((d[ex]*d[tier]).sum())
            OR,p,_=logit(f'{tier} ~ {ex} + age + C(sexc){pcterm}',d,ex)
            ORu,pu,_=logit(f'{tier} ~ {ex} + age + C(sexc) + log_util{pcterm}',d,ex)
            S['adjusted'].append({'tier':tier,'exposure':ex,'case_carriers':n_ex,'severe':n_sev,'OR_16PC':OR,'p_16PC':p,'OR_util':ORu,'p_util':pu})
            print(f"   [{tier:8s}] {ex:15s} n={n_ex:>4} sev={n_sev:>3} | +16PC OR={str(OR):>7} p={str(round(p,4) if p is not None else p):>8} | +util OR={str(ORu):>7} p={str(round(pu,4) if pu is not None else pu):>8}")
    print("\n== READ: does R85H->severe survive ancestry(16PC/within-AFR)+utilization? if yes, R85H is a SEVERITY modifier ==")
    print("\n===== SEVERITY ANCESTRY CONTROL (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
