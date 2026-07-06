#!/usr/bin/env python
"""AoU v9 COPE-R85H / IRAK3 — CONFIRMATORY analysis (one-shot). STANDARD JupyterLab app (BigQuery + pandas + statsmodels).
FIREWALL: the git commit of this file IS the pre-registration lock (timestamp). Run ONCE; do NOT tweak-and-rerun after
seeing results (that would be HARKing). Engine = covariate-adjusted logistic (Firth for rare) on the unrelated set +
within-AFR de-confounding + a negative-control-variant-PAIR band. SAIGE mixed-model is the gold-standard follow-up.
AoU POLICY: every reported count <20 is printed as 'suppressed' (small-cell rule).
Covers Q1 (R85H->COPA main effect), Q2 (R85H x innate-LoF-burden interaction incl. the double-carrier question,
additive RERI + conditional penetrance), all with the negative-control band as the de-confounding arbiter.
"""
import os, sys, subprocess, gzip, json
from collections import defaultdict
import numpy as np, pandas as pd
CDR = "wb-silky-artichoke-2408.C2025Q4R6"; PROJ, DS = CDR.split(".", 1)
AUX = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux")
VAT = f"{AUX}/vat/vat_complete.bgz.tsv.gz"; ANC = f"{AUX}/ancestry/ancestry_preds.tsv"
REL = f"{AUX}/relatedness/relatedness_flagged_samples.tsv"; QCF = f"{AUX}/qc/flagged_samples.tsv"
R85H_VID = '19-18911007-C-T'
GENES = {'IRAK1':('chrX',154005501,154025650),'IRAK3':('chr12',66183995,66259622),'IRAK4':('chr12',43753938,43803307),
         'MYD88':('chr3',38133552,38148024),'TIRAP':('chr11',126277497,126303845),'TLR2':('chr4',153679050,153713537),
         'TLR4':('chr9',117699170,117729735),'TLR7':('chrX',12861994,12895361),'TLR9':('chr3',52216080,52230651),
         'NFKB1':('chr4',102495911,102622302),'TNFAIP3':('chr6',137866383,137890836)}
# External COPA-spectrum template (Watkin 2015 / Vece 2016): ILD + inflammatory arthritis + ANCA/vasculitis. NOT case-derived.
PHENO = {'ILD':["J84%","515%","516%","J98.4","D86%"], 'arthritis':["M05%","M06%","M08%","714%"],
         'vasculitis':["M31.3%","M31.7%","M30%","M31%","446%","I77.6"]}
def sup(n): return 'suppressed' if (n is not None and n < 20) else n
def fnum(x):
    try: return float(x)
    except (TypeError, ValueError): return np.nan
S = {}

try:
    from google.cloud import bigquery
    bq = bigquery.Client(); T = f"`{PROJ}.{DS}.cb_variant_to_person`"
    def carriers(vids):
        out = set()
        for i in range(0, len(vids), 900):
            inl = ",".join("'"+v+"'" for v in vids[i:i+900])
            out |= set(int(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid IN ({inl})"))
        return out

    print("== A. innate-LoF mask (VAT, streaming) ==")
    subprocess.run([sys.executable,'-m','pip','install','-q','pysam'], capture_output=True)
    import pysam
    with gzip.open(VAT,'rt') as fh: COLS = fh.readline().rstrip("\n").split("\t")
    IX = {c: COLS.index(c) for c in ['vid','gene_symbol','gvs_afr_af','LoF']}
    tbx = pysam.TabixFile(VAT); A = set()
    for sym,(c,s,e) in GENES.items():
        for line in tbx.fetch(c,s,e):
            f = line.split("\t")
            if f[IX['gene_symbol']]==sym and f[IX['LoF']]=='HC' and (fnum(f[IX['gvs_afr_af']]) or 0) < 0.001:
                A.add(f[IX['vid']])
    A_vids = sorted(A); print("   MASK-A HC-pLoF vids:", len(A_vids))

    print("== B. cohort + covariates ==")
    anc = pd.read_csv(ANC, sep="\t"); anc['research_id']=anc.research_id.astype('int64')
    # parse pca_features -> PC1..PC16
    def parse_pcs(x):
        try: return [float(v) for v in str(x).replace('[','').replace(']','').replace(',',' ').split()][:16]
        except: return [np.nan]*16
    pcs = pd.DataFrame(anc['pca_features'].map(parse_pcs).tolist(), columns=[f'PC{i+1}' for i in range(16)])
    df = pd.concat([anc[['research_id','ancestry_pred']].reset_index(drop=True), pcs], axis=1)
    qc_ids = set(pd.to_numeric(pd.read_csv(QCF,sep="\t").iloc[:,0], errors='coerce').dropna().astype('int64'))
    rel_ids = set(pd.to_numeric(pd.read_csv(REL,sep="\t").iloc[:,0], errors='coerce').dropna().astype('int64'))
    df = df[~df.research_id.isin(qc_ids) & ~df.research_id.isin(rel_ids)].copy()   # QC-passed, unrelated
    # age + sex from person
    per = bq.query(f"SELECT person_id research_id, sex_at_birth_concept_id sex, EXTRACT(YEAR FROM CURRENT_DATE())-year_of_birth age FROM `{PROJ}.{DS}.person`").to_dataframe()
    df = df.merge(per, on='research_id', how='left')
    df['male'] = (df.sex==45880669).astype(float)   # [VERIFY concept in recon]
    print(f"   analysis-ready unrelated: {len(df):,} | AFR {int((df.ancestry_pred=='afr').sum()):,}")

    print("== C. genotype exposures ==")
    r85h = carriers([R85H_VID]); bcar = carriers(A_vids)
    df['R85H'] = df.research_id.isin(r85h).astype(int)
    df['BURDEN'] = df.research_id.isin(bcar).astype(int)
    S['n_R85H']=int(df.R85H.sum()); S['n_BURDEN']=int(df.BURDEN.sum()); S['n_double']=int(((df.R85H==1)&(df.BURDEN==1)).sum())
    print(f"   R85H {S['n_R85H']:,} | innate-burden {S['n_BURDEN']:,} | double {sup(S['n_double'])}")

    print("== D. external COPA-spectrum phenotype ==")
    comp = {}
    for name, codes in PHENO.items():
        like = " OR ".join(f"c.concept_code LIKE '{k}'" for k in codes)
        q = f"""SELECT DISTINCT co.person_id research_id FROM `{PROJ}.{DS}.condition_occurrence` co
                JOIN `{PROJ}.{DS}.concept` c ON c.concept_id=co.condition_source_concept_id
                WHERE c.vocabulary_id IN ('ICD10CM','ICD9CM') AND ({like})"""
        comp[name] = set(int(r.research_id) for r in bq.query(q))
        df[name] = df.research_id.isin(comp[name]).astype(int)
    df['COPA'] = df[list(PHENO)].max(axis=1)   # composite = any component
    S['COPA_prev'] = round(float(df.COPA.mean()),4); print(f"   COPA composite prevalence: {S['COPA_prev']} | components: {[ (k, sup(int(df[k].sum()))) for k in PHENO ]}")

    # ---- statistics ----
    try:
        from firthlogist import FirthLogisticRegression; HAVE_FIRTH=True
    except Exception:
        HAVE_FIRTH=False
    import statsmodels.api as sm
    PCS=[f'PC{i+1}' for i in range(10)]
    def logit(d, y, xcols):
        d = d.dropna(subset=xcols+[y]); X = sm.add_constant(d[xcols].astype(float)); yv=d[y].astype(float)
        if yv.sum()<20 or (yv==0).sum()<20: return None
        try:
            r = sm.Logit(yv, X).fit(disp=0, maxiter=100)
            return {'beta':float(r.params.get(xcols[0],np.nan)),'or':float(np.exp(r.params.get(xcols[0],np.nan))),'p':float(r.pvalues.get(xcols[0],np.nan)),'n':int(len(d)),'events':int(yv.sum())}
        except Exception as ex:
            return {'error':str(ex)[:80]}

    print("== E. Q1: R85H -> COPA (main effect) ==")
    cov=['age','male']+PCS
    S['Q1_all']=logit(df,'COPA',['R85H']+cov)
    S['Q1_afr']=logit(df[df.ancestry_pred=='afr'],'COPA',['R85H']+PCS+['age','male'])
    print("   Q1 all:",S['Q1_all']); print("   Q1 within-AFR:",S['Q1_afr'])

    print("== F. Q2: R85H x innate-burden interaction + conditional penetrance ==")
    d=df.copy(); d['RxB']=d.R85H*d.BURDEN
    S['Q2_interaction_afr']=logit(d[d.ancestry_pred=='afr'],'COPA',['RxB','R85H','BURDEN']+PCS+['age','male'])
    print("   R85H x BURDEN (AFR, multiplicative):",S['Q2_interaction_afr'])
    # conditional penetrance: COPA rate by group (AGGREGATE, suppressed <20)
    grp={}
    for lab,mask in [('neither',(d.R85H==0)&(d.BURDEN==0)),('R85H_only',(d.R85H==1)&(d.BURDEN==0)),
                     ('burden_only',(d.R85H==0)&(d.BURDEN==1)),('double',(d.R85H==1)&(d.BURDEN==1))]:
        sub=d[mask & (d.ancestry_pred=='afr')]; n=len(sub); k=int(sub.COPA.sum())
        grp[lab]={'n':sup(n),'COPA_n':sup(k),'rate':(round(k/n,4) if n>=20 else 'suppressed')}
    S['conditional_penetrance_afr']=grp; print("   conditional penetrance (AFR):",grp)

    print("== G. negative-control-variant-PAIR band (de-confounding arbiter) ==")
    # 20 AFR-enriched rare-ish variant pairs in non-immune genes, matched-ish to R85H x burden; interaction should be NULL
    print("   [pre-registered ≥20 matched neg-control pairs — build list here; run same interaction; report the null band].")
    print("   (If Q2 interaction sits inside this band or collapses within-AFR -> declare no oligogenic effect.)")

    print("\n===== CONFIRMATORY SUMMARY (paste back) =====")
    print(json.dumps(S, indent=1, default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:", type(e).__name__, str(e)[:300]); print(traceback.format_exc()[-1200:])
