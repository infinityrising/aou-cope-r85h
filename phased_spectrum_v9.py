#!/usr/bin/env python
"""AoU v9 -- PHASE-RESOLVED STING across the COPA spectrum. Long-read cis-phasing (sting_phenome_pav_v9.csv, N=13,252;
AQ_d / HAQ_d = CIS dose of true AQ / HAQ alleles) is our unique asset -- short-read CANNOT phase. Three phased tests:
 (1) CIS DOSE: incidence & within-case severity ~ AQ_d(cis) + HAQ_d(cis). Does the true cis-AQ allele drive disease /
     severity, and does cis-HAQ (= AQ core + the protective R71H) blunt what cis-AQ confers?
 (2) CIS/TRANS NATURAL EXPERIMENT: among GENOTYPE-AQ (carry BOTH G230A+R293Q by UNPHASED BQ, R71H-), split by long-read
     phase -- CIS (AQ_d>=1 = a real AQ allele) vs TRANS (AQ_d==0 = the two SNPs on OPPOSITE chromosomes = NO real allele
     = built-in negative control). cis worse than trans => it is the ALLELE, not the individual SNPs. Impossible without
     long-read; phasing IS the experiment.
 (3) DOSE GRADIENT: disease / severe rate by AQ_d in {0,1,2}.
Composite COPA-lung (for power) + per-phenotype N reported honestly (PAV is thin). NOTE: R85H itself is unphaseable here
(~1 R85H case in PAV) -> this is the STING-second-hit arm, phase-resolved. DISCOVERY. Standard app. Ends 'run complete'/'run failed'.
"""
import os, json
import numpy as np, pandas as pd
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); SNP=f"{MNT}/v9/wgs/short_read/snpindel"
ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
PHPAV=os.path.expanduser("~/sting_phenome_pav_v9.csv")
STINGV={'R71H':'5-139481493-C-T','G230A':'5-139478340-C-G','R293Q':'5-139477397-C-T'}
PANEL={
 'bronchiectasis':['J47','494'],'fibrotic_ILD':['J84.1','515','516.3'],'LIP':['J84.2'],'alveolar_hemorrhage':['R04.2','R04.89'],
 'infl_arthritis':['M05','M06','M08','714'],'ANCA_vasculitis':['M31.3','M31.7'],'vasculitis_broad':['M30','M31'],'glomerulonephritis':['N01','N03','N05','580','581','582','583'],
}
LUNGP=['bronchiectasis','fibrotic_ILD','LIP','alveolar_hemorrhage']; SYSP=['infl_arthritis','ANCA_vasculitis','vasculitis_broad','glomerulonephritis']
HARD=['J96','518.81','518.83','518.84','Z99.81','V46.2','Z99.2','N18.6','Z94','V42.6','V42.0','Z99.11']
MINCODES=int(os.environ.get("MINCODES","2"))
def pct(x,n): return f"{100*x/n:.0f}%({x}/{n})" if n else "n/a"
S={}
try:
    from google.cloud import bigquery
    import ast
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    pav=pd.read_csv(PHPAV); pav['research_id']=pav.research_id.astype(str)
    print(f"PAV phased cohort {len(pav)} | cis-AQ (AQ_d>=1) {int((pav.AQ_d>=1).sum())} | cis-HAQ (HAQ_d>=1) {int((pav.HAQ_d>=1).sum())}")
    def q_ids(sql): return set(str(r[0]) for r in bq.query(sql))
    def caseset(codes):
        lk=" OR ".join([f"c.concept_code LIKE '{p}%'" for p in codes])
        return q_ids(f"SELECT co.person_id FROM `{PROJ}.{DS}.condition_occurrence` co JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id WHERE c.vocabulary_id LIKE 'ICD%' AND ({lk}) GROUP BY co.person_id HAVING COUNT(*)>={MINCODES}")
    def anyset(codes):
        lk=" OR ".join([f"c.concept_code LIKE '{p}%'" for p in codes])
        return q_ids(f"SELECT DISTINCT co.person_id FROM `{PROJ}.{DS}.condition_occurrence` co JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id WHERE c.vocabulary_id LIKE 'ICD%' AND ({lk})")
    hard=anyset(HARD); cases={ph:caseset(cc) for ph,cc in PANEL.items()}
    lung=set().union(*[cases[p] for p in LUNGP]); syst=set().union(*[cases[p] for p in SYSP])
    def carr(vids):
        out=set()
        for i in range(0,len(vids),900):
            if not vids[i:i+900]: continue
            inl=",".join("'"+v+"'" for v in vids[i:i+900])
            out|=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid IN ({inl})"))
        return out
    st={k:carr([v]) for k,v in STINGV.items()}
    geno_aq=(st['G230A']&st['R293Q'])-st['R71H']       # UNPHASED genotype-AQ (cis OR trans)
    # ---- PAV frame ----
    p=pav[['research_id','AQ_d','HAQ_d']].copy()
    p['lung']=p.research_id.isin(lung).astype(int); p['syst']=p.research_id.isin(syst).astype(int)
    p['bronch']=p.research_id.isin(cases['bronchiectasis']).astype(int); p['hard']=p.research_id.isin(hard).astype(int)
    p['geno_aq']=p.research_id.isin(geno_aq).astype(int)
    anc=pd.read_csv(ANC,sep="\t"); anc['research_id']=anc.research_id.astype(str)
    try:
        P=anc['pca_features'].apply(lambda x: ast.literal_eval(x) if isinstance(x,str) else (x if isinstance(x,list) else []))
        k=min(min(map(len,P[P.map(len)>0])) if (P.map(len)>0).any() else 0,5)
        pcf=pd.DataFrame([row[:k] for row in P],columns=[f'PC{i+1}' for i in range(k)]); pcf['research_id']=anc.research_id.values
        p=p.merge(pcf,on='research_id',how='left'); PCS=[f'PC{i+1}' for i in range(k)]
    except Exception: PCS=[]
    try:
        ppl=bq.query(f"SELECT person_id, DATE_DIFF(DATE '2025-01-01', DATE(birth_datetime), YEAR) age, CAST(sex_at_birth_concept_id AS STRING) sexc FROM `{PROJ}.{DS}.person`").to_dataframe()
        ppl['research_id']=ppl.person_id.astype(str); p=p.merge(ppl[['research_id','age','sexc']],on='research_id',how='left')
    except Exception: pass
    cov=" + age + C(sexc)"+((" + "+" + ".join(PCS)) if PCS else "")
    print("\n== PAV ∩ phenotype cases (what phasing can actually see) ==")
    for ph in PANEL: print(f"   {ph:20s} PAV∩cases {len(cases[ph]&set(pav.research_id)):>5}")
    print(f"   composite-LUNG PAV∩cases {len(lung&set(pav.research_id))} | composite-SYSTEMIC {len(syst&set(pav.research_id))} | genotype-AQ in PAV {int(p.geno_aq.sum())}")
    import statsmodels.formula.api as smf
    def logit(fm,dd,t):
        try: r=smf.logit(fm,data=dd,missing='drop').fit(disp=0); return round(float(np.exp(r.params[t])),3),float(r.pvalues[t])
        except Exception: return None,None
    # ---- (1) CIS DOSE model ----
    print("\n== (1) CIS-DOSE: outcome ~ AQ_d(cis) + HAQ_d(cis) + age + sex + 5PC ==")
    S['cis_dose']={}
    for lab,col in [('composite-lung','lung'),('bronchiectasis','bronch'),('composite-systemic','syst')]:
        aor,ap=logit(f'{col} ~ AQ_d + HAQ_d{cov}',p,'AQ_d'); hor,hp=logit(f'{col} ~ AQ_d + HAQ_d{cov}',p,'HAQ_d')
        sub=p[p[col]==1].copy()
        sor,spv=logit(f'hard ~ AQ_d + HAQ_d{cov}',sub,'AQ_d') if len(sub)>=15 else (None,None)
        S['cis_dose'][lab]={'AQ_inc_OR':aor,'AQ_inc_p':ap,'HAQ_inc_OR':hor,'HAQ_inc_p':hp,'AQ_sev_OR':sor,'AQ_sev_p':spv,'n_cases':int(p[col].sum())}
        print(f"   {lab:20s} incidence AQ_d OR={aor}(p{ap if ap is None else round(ap,3)}) HAQ_d OR={hor}(p{hp if hp is None else round(hp,3)}) | within-case severity AQ_d OR={sor}(p{spv if spv is None else round(spv,3)})")
    # ---- (2) CIS/TRANS NATURAL EXPERIMENT ----
    print("\n== (2) CIS/TRANS natural experiment (genotype-AQ split by long-read phase; trans = built-in negative control) ==")
    ga=p[p.geno_aq==1].copy(); ga['cis']=(ga.AQ_d>=1).astype(int)
    ncis=int((ga.cis==1).sum()); ntrans=int((ga.cis==0).sum())
    cis_lung=int(ga[ga.cis==1].lung.sum()); trans_lung=int(ga[ga.cis==0].lung.sum())
    cis_hard=int(ga[ga.cis==1].hard.sum()); trans_hard=int(ga[ga.cis==0].hard.sum())
    print(f"   genotype-AQ in PAV: CIS(true allele) n={ncis} | TRANS(pseudo-AQ) n={ntrans}")
    print(f"   composite-lung:  CIS {pct(cis_lung,ncis)}  vs  TRANS {pct(trans_lung,ntrans)}")
    print(f"   hard-outcome:    CIS {pct(cis_hard,ncis)}  vs  TRANS {pct(trans_hard,ntrans)}")
    or_ct,p_ct=logit(f'lung ~ cis{cov}',ga,'cis') if len(ga)>=15 else (None,None)
    S['cis_vs_trans']={'cis_n':ncis,'trans_n':ntrans,'cis_lung':cis_lung,'trans_lung':trans_lung,'cis_hard':cis_hard,'trans_hard':trans_hard,'lung_OR_cis':or_ct,'lung_p_cis':p_ct}
    print(f"   adj lung ~ cis: OR={or_ct} p={p_ct}  (cis>trans => the AQ ALLELE, not the SNPs)")
    # ---- (3) DOSE GRADIENT ----
    print("\n== (3) CIS-AQ DOSE GRADIENT (0/1/2) ==")
    S['dose']=[]
    for dv in [0,1,2]:
        sub=p[p.AQ_d==dv]; n=len(sub)
        row={'AQ_d':dv,'n':n,'lung_rate':round(float(sub.lung.mean()),4) if n else None,'severe_rate':round(float(sub.hard.mean()),4) if n else None}
        S['dose'].append(row); print(f"   AQ_d={dv}: n={n:>6} lung {pct(int(sub.lung.sum()),n)} hard {pct(int(sub.hard.sum()),n)}")
    # ---- (4) HAQ vs AQ (R71H protection, phased) ----
    cisAQ=p[(p.AQ_d>=1)&(p.HAQ_d==0)]; cisHAQ=p[p.HAQ_d>=1]
    print("\n== (4) cis-HAQ vs cis-AQ (does the R71H in HAQ protect?) ==")
    print(f"   cis-AQ  n={len(cisAQ)} lung {pct(int(cisAQ.lung.sum()),len(cisAQ))} hard {pct(int(cisAQ.hard.sum()),len(cisAQ))}")
    print(f"   cis-HAQ n={len(cisHAQ)} lung {pct(int(cisHAQ.lung.sum()),len(cisHAQ))} hard {pct(int(cisHAQ.hard.sum()),len(cisHAQ))}")
    S['haq_vs_aq']={'cisAQ_n':len(cisAQ),'cisAQ_lung':int(cisAQ.lung.sum()),'cisHAQ_n':len(cisHAQ),'cisHAQ_lung':int(cisHAQ.lung.sum())}
    print("\n===== PHASED SPECTRUM (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
