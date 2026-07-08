#!/usr/bin/env python
"""AoU v9 -- THE DECISIVE de-confounder: the negative-control-variant BAND (pre-registered v5 §5.6/§6/§8 as the decisive
test; flagged PENDING; never built until now). Question: is R85H's bronchiectasis-SEVERITY effect SPECIFIC, or a founder-
tilt / AFR-ascertainment artifact? Calibrate R85H against a NULL BAND of ~30 AFR-allele-frequency-MATCHED, functionally
NEUTRAL (synonymous/intron/intergenic) variants pushed through the IDENTICAL within-case severity pipeline.
Also fixes two appraisal findings simultaneously: severity is restricted to RESPIRATORY organ-failure (resp failure /
O2-dependence / LUNG transplant / mech-vent -- NO renal codes, so the APOL1 confound is removed), and this LOCKS ONE
severity definition. VERDICT: R85H severity OR beyond the band's 95th percentile -> SPECIFIC; inside the band -> ARTIFACT
(declare null per the lock). PRE-SPECIFIED. Standard app. Ends 'run complete' / 'run failed'.
"""
import os, gzip, ast, json
import numpy as np, pandas as pd
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); SNP=f"{MNT}/v9/wgs/short_read/snpindel"
ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
R85H_VID='19-18911007-C-T'
BRONCH=['J47','494']
# LOCKED, organ-appropriate severity: RESPIRATORY organ failure ONLY (no renal ESRD/dialysis -> removes APOL1 confound)
RESP_HARD=['J96','518.81','518.83','518.84','Z99.81','V46.2','Z94.2','J95.85','Z99.11']
AFR_LO,AFR_HI=0.008,0.02      # match R85H AFR AF ~1.3%
EUR_MAX=0.003                 # AFR-enriched / founder-like (as R85H is)
NEUTRAL=('synonymous','intron','intergenic','non_coding','upstream','downstream','5_prime_utr','3_prime_utr','5_prime_UTR','3_prime_UTR')
NONNEUTRAL=('missense','stop','frameshift','splice','start_lost','stop_lost','inframe')
WINDOWS=[('1',60000000,65000000),('2',60000000,65000000),('3',70000000,75000000),('4',60000000,65000000),('5',60000000,65000000),('7',60000000,65000000),('8',60000000,65000000),('9',80000000,85000000),('10',60000000,65000000),('11',70000000,75000000),('13',50000000,55000000),('18',30000000,35000000),('20',30000000,35000000)]
EXCLUDE_GENES={'COPA','COPB1','COPB2','COPG1','COPG2','COPZ1','COPZ2','COPE','ARCN1','STING1','TMEM173','IRAK1','IRAK3','IRAK4','MYD88'}
PER_WINDOW=int(os.environ.get("PER_WINDOW","3")); MINCODES=int(os.environ.get("MINCODES","2"))
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
    def carr(vids):
        out=set()
        for i in range(0,len(vids),900):
            if not vids[i:i+900]: continue
            inl=",".join("'"+v+"'" for v in vids[i:i+900])
            out|=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid IN ({inl})"))
        return out
    bronch=caseset(BRONCH); resp=anyset(RESP_HARD)
    print(f"bronchiectasis {len(bronch)} | respiratory-hard-severe {len(bronch&resp)} ({100*len(bronch&resp)/max(len(bronch),1):.0f}%)")
    # ---- VAT setup + sample AFR-AF-matched NEUTRAL negative-control variants ----
    with gzip.open(VAT,'rt') as fh: COLS=fh.readline().rstrip("\n").split("\t")
    has_eur='gvs_eur_af' in COLS
    IX={c:COLS.index(c) for c in (['vid','gene_symbol','consequence','gvs_afr_af']+(['gvs_eur_af'] if has_eur else []))}
    tbx=pysam.TabixFile(VAT); ctrl=[]
    for (ch,s,e) in WINDOWS:
        found=0
        try: it=tbx.fetch('chr'+ch,s,e)
        except Exception: continue
        for line in it:
            f=line.split("\t"); cons=f[IX['consequence']]
            if any(x in cons for x in NONNEUTRAL): continue
            if not any(x in cons for x in NEUTRAL): continue
            if f[IX['gene_symbol']] in EXCLUDE_GENES: continue
            af=fnum(f[IX['gvs_afr_af']])
            if not (AFR_LO<=af<=AFR_HI): continue
            if has_eur and fnum(f[IX['gvs_eur_af']])>EUR_MAX: continue
            ctrl.append(f[IX['vid']]); found+=1
            if found>=PER_WINDOW: break
    print(f"AFR-AF-matched neutral negative-control variants sampled: {len(ctrl)} (target ~{PER_WINDOW*len(WINDOWS)}; eur-filter={has_eur})")
    S['n_controls_sampled']=len(ctrl)
    # ---- case frame (bronchiectasis) + covariates ----
    anc=pd.read_csv(ANC,sep="\t"); anc['research_id']=anc.research_id.astype(str)
    d=anc[anc.research_id.isin(bronch)][['research_id','ancestry_pred']].copy(); d['SEV']=d.research_id.isin(resp).astype(int)
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
    afr=d[d.ancestry_pred=='afr'].copy(); pc5=(" + "+" + ".join(PCS[:5])) if PCS else ""; pc16=(" + "+" + ".join(PCS)) if PCS else ""
    print(f"bronchiectasis case-frame {len(d)} | AFR cases {len(afr)} | PCs {len(PCS)}")
    import statsmodels.formula.api as smf
    def logit(fm,dd,t):
        try: r=smf.logit(fm,data=dd,missing='drop').fit(disp=0); return round(float(np.exp(r.params[t])),3),float(r.pvalues[t])
        except Exception: return None,None
    def sev_within_afr(cset):
        a=afr.copy(); a['V']=a.research_id.isin(cset).astype(int); n=int(a.V.sum())
        if n<3: return None,None,n,int((a.V*a.SEV).sum())
        OR,p=logit(f'SEV ~ V + age + C(sexc){pc5}',a,'V'); return OR,p,n,int((a.V*a.SEV).sum())
    def sev_allcases(cset):
        a=d.copy(); a['V']=a.research_id.isin(cset).astype(int)
        if int(a.V.sum())<3: return None,None
        return logit(f'SEV ~ V + age + C(sexc){pc16}',a,'V')
    # ---- R85H reference ----
    r_or,r_p,r_n,r_sev=sev_within_afr(carr([R85H_VID])); ra_or,ra_p=sev_allcases(carr([R85H_VID]))
    S['R85H']={'withinAFR_OR':r_or,'withinAFR_p':r_p,'withinAFR_n':r_n,'withinAFR_severe':r_sev,'allcases16PC_OR':ra_or,'allcases16PC_p':ra_p}
    print(f"\n★ R85H (respiratory-hard severity): within-AFR OR={r_or} p={r_p} (n={r_n}, severe={r_sev}) | all-cases+16PC OR={ra_or} p={ra_p}")
    # ---- NEGATIVE-CONTROL BAND ----
    print(f"\n== NULL BAND: {len(ctrl)} AFR-AF-matched neutral variants through the identical pipeline ==")
    band=[]
    for v in ctrl:
        cs=carr([v]); o,p,n,sev=sev_within_afr(cs); ao,ap=sev_allcases(cs)
        if o is not None: band.append({'vid':v,'withinAFR_OR':o,'allcases_OR':ao,'n':n})
    S['band']=band; ors=sorted([b['withinAFR_OR'] for b in band]); aors=sorted([b['allcases_OR'] for b in band if b['allcases_OR'] is not None])
    def pctile_of(val,arr): return round(100*sum(1 for x in arr if x<val)/len(arr),1) if arr else None
    if ors:
        q=lambda a,p: round(float(np.percentile(a,p)),3)
        S['band_summary']={'n':len(ors),'median':q(ors,50),'p90':q(ors,90),'p95':q(ors,95),'max':round(max(ors),3),'R85H_percentile':pctile_of(r_or,ors) if r_or else None}
        print(f"   within-AFR band (n={len(ors)}): median OR={q(ors,50)} | 90th={q(ors,90)} | 95th={q(ors,95)} | max={round(max(ors),3)}")
        print(f"   ★ R85H within-AFR OR={r_or} sits at the {pctile_of(r_or,ors)}th percentile of the null band")
        if aors and ra_or: print(f"   all-cases band (n={len(aors)}): median={round(float(np.median(aors)),3)} 95th={round(float(np.percentile(aors,95)),3)} | R85H OR={ra_or} = {pctile_of(ra_or,aors)}th pct")
        verdict=("SPECIFIC (R85H beyond 95th pct of the founder-matched null -> not an ancestry artifact)" if (r_or and pctile_of(r_or,ors) is not None and pctile_of(r_or,ors)>=95) else "INSIDE THE BAND -> R85H severity NOT distinguishable from founder-matched neutral variants -> ARTIFACT / declare null per lock")
        S['verdict']=verdict; print(f"\n== VERDICT: {verdict} ==")
    else: print("   band empty (no controls converged) — widen AF window or PER_WINDOW"); S['verdict']='band empty'
    print("\n===== NEGATIVE-CONTROL SEVERITY BAND (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
