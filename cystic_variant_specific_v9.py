#!/usr/bin/env python
"""AoU v9 -- VARIANT-SPECIFIC exposure x RESOLVED lung phenotype (cystic vs fibrotic). Corrects the two blind spots in
the reverse-burden (run 037):
 (A) EXPOSURE. Gene-level burden + a 'damaging = REVEL>0.5' filter DELETES gain-of-function disease alleles: SAVI-STING
     and COPA-WD40 GoF mutations score LOW REVEL (REVEL is trained toward LoF), so the burden mask drops the very
     variants that cause the disease -- which is WHY STING1/COPA came back null. Here we POSITION-SELECT the actual
     pathogenic residues (SAVI hotspots; COPA WD40 R233/E241/D243 cluster), any missense, REVEL-independent.
 (B) PHENOTYPE. J84/515/516 = FIBROTIC ILD. The index patient's lesion is CYSTIC: airway + lymphoid inflammation and
     recurrent infection remodel into cysts (bronchiectasis, acquired lung cyst, LIP, LAM/PLCH). AoU has NO curated
     phenotype -- only raw ICD -- so separating cystic from fibrotic is on us. We define both and test each exposure
     against both, and read out the SIGN.
PRE-SPECIFIED (locked before run): GoF positive controls (SAVI, COPA) -> RISK for lung disease that burden missed;
 a-priori lean COPA-class/airway -> CYSTIC, SAVI -> FIBROTIC. KEY TEST: IRAK3-LoF and R85H+2nd-hit -> RISK for CYSTIC,
 REVERSING the PROTECTIVE sign they showed on FIBROTIC ILD (interferon = anti-fibrotic but pro-airway-remodeling).
 DISCOVERY / exploratory; own confirmation required. Standard app. Ends 'run complete' / 'run failed'.
"""
import os, gzip, re, ast, json
import numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); SNP=f"{MNT}/v9/wgs/short_read/snpindel"
ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
R85H_VID='19-18911007-C-T'
STINGV={'R71H':'5-139481493-C-T','G230A':'5-139478340-C-G','R293Q':'5-139477397-C-T'}
MINCODES=int(os.environ.get("MINCODES","2")); AFMAX=float(os.environ.get("AFMAX","0.01"))
# ---- variant-specific GoF exposure sets (position-selected on the pathogenic residues, REVEL-independent) ----
GOF={
 'SAVI_STING':{'gene':'STING1','chr':'5','s':139475528,'e':139482935,'res':{147,154,155,158,166,206,281,284}},
 'COPA_WD40' :{'gene':'COPA', 'chr':'1','s':160288580,'e':160343566,'res':{230,233,236,241,242,243}},
}
# ---- resolved lung phenotypes (CYSTIC = airway/lymphoid remodeling; FIBROTIC = classic ILD comparator) ----
CYSTIC={'bronchiectasis':['J47','494'],'lung_cyst':['J98.4'],'LAM_PLCH':['J84.81','J84.82'],'LIP':['J84.2']}
FIBROTIC={'pulm_fibrosis':['J84.1','515','516.3']}
RESID=re.compile(r'p\.\(?[A-Za-z]{3}(\d+)')
def fnum(x):
    try: return float(x)
    except (TypeError,ValueError): return 0.0
S={'params':{'MINCODES':MINCODES,'AFMAX':AFMAX}}
try:
    import pysam
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    with gzip.open(VAT,'rt') as fh: COLS=fh.readline().rstrip("\n").split("\t")
    hasaa='aa_change' in COLS
    IX={c:COLS.index(c) for c in ['vid','gene_symbol','consequence','gvs_afr_af','revel','LoF']+(['aa_change'] if hasaa else [])}
    tbx=pysam.TabixFile(VAT)
    def carr(vids):
        out=set()
        for i in range(0,len(vids),900):
            chunk=vids[i:i+900]
            if not chunk: continue
            inl=",".join("'"+v+"'" for v in chunk)
            out|=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid IN ({inl})"))
        return out
    # ---- (A) build GoF exposure vids BY RESIDUE (the fix: no REVEL gate, no gene-collapse) ----
    print("== (A) variant-specific GoF exposure (position-selected, REVEL-independent) ==")
    expo={}; revel_dropped={}
    for name,g in GOF.items():
        vids=[]; low_revel=0
        for line in tbx.fetch('chr'+g['chr'],g['s'],g['e']):
            f=line.split("\t")
            if f[IX['gene_symbol']]!=g['gene'] or 'missense' not in f[IX['consequence']]: continue
            if fnum(f[IX['gvs_afr_af']])>=AFMAX: continue
            m=RESID.search(f[IX['aa_change']] if hasaa else '')
            if m and int(m.group(1)) in g['res']:
                vids.append(f[IX['vid']])
                if fnum(f[IX['revel']])<=0.5: low_revel+=1   # would have been dropped by the burden mask
        expo[name]=carr(vids); revel_dropped[name]=low_revel
        S[f'{name}_vids']=len(vids); S[f'{name}_revel_le0.5']=low_revel
        print(f"   {name}: {len(vids)} residue-matched vids ({low_revel} with REVEL<=0.5 that burden would DROP) -> {len(expo[name])} carriers")
    # ---- IRAK3 truncating-LoF + R85H + AQ(unphased proxy) + the OA37-matched combos ----
    irak3=[]
    for line in tbx.fetch('chr12',66183995,66259622):
        f=line.split("\t")
        if f[IX['gene_symbol']]=='IRAK3' and f[IX['LoF']]=='HC' and fnum(f[IX['gvs_afr_af']])<0.001: irak3.append(f[IX['vid']])
    C_irak3=carr(irak3); C_r85h=carr([R85H_VID]); stcar={k:carr([v]) for k,v in STINGV.items()}
    C_aq=(stcar['G230A']&stcar['R293Q'])-stcar['R71H']                 # AQ pattern (unphased proxy; cis needs long-read)
    expo['IRAK3_LoF']=C_irak3; expo['R85H']=C_r85h
    expo['R85H_x_IRAK3LoF']=C_r85h&C_irak3; expo['R85H_x_AQ']=C_r85h&C_aq
    print(f"   IRAK3-LoF {len(C_irak3)} | R85H {len(C_r85h)} | R85H&IRAK3LoF {len(expo['R85H_x_IRAK3LoF'])} | R85H&AQ {len(expo['R85H_x_AQ'])}")
    # ---- (B) resolved phenotype case sets (>=MINCODES occurrences = chronic/established) ----
    def cases(codes):
        lk=" OR ".join([f"c.concept_code LIKE '{p}%'" for p in codes])
        q=(f"SELECT co.person_id FROM `{PROJ}.{DS}.condition_occurrence` co "
           f"JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id "
           f"WHERE c.vocabulary_id LIKE 'ICD%' AND ({lk}) GROUP BY co.person_id HAVING COUNT(*)>={MINCODES}")
        return set(str(r.person_id) for r in bq.query(q))
    print("== (B) resolved phenotypes ==")
    PH={'CYSTIC':set().union(*[cases(v) for v in CYSTIC.values()]),'FIBROTIC':cases(FIBROTIC['pulm_fibrosis'])}
    for k,v in CYSTIC.items(): PH[k]=cases(v)            # component breakdown (bronchiectasis etc.)
    S['pheno_n']={k:len(v) for k,v in PH.items()}
    print("   cases:",S['pheno_n'])
    # ---- cohort + covariates (16 PC + age + sex + smoking) ----
    anc=pd.read_csv(ANC,sep="\t"); anc['research_id']=anc.research_id.astype(str)
    d=anc[['research_id']].copy(); COV=[]
    try:
        P=anc['pca_features'].apply(lambda x: ast.literal_eval(x) if isinstance(x,str) else (x if isinstance(x,list) else []))
        k=min(min(map(len,P[P.map(len)>0])) if (P.map(len)>0).any() else 0,16)
        p=pd.DataFrame([row[:k] for row in P],columns=[f'PC{i+1}' for i in range(k)]); p['research_id']=anc.research_id.values
        d=d.merge(p,on='research_id',how='left'); COV+=[f'PC{i+1}' for i in range(k)]
    except Exception: pass
    try:
        ppl=bq.query(f"SELECT person_id, DATE_DIFF(DATE '2025-01-01', DATE(birth_datetime), YEAR) age, CAST(sex_at_birth_concept_id AS STRING) sexc FROM `{PROJ}.{DS}.person`").to_dataframe()
        ppl['research_id']=ppl.person_id.astype(str); d=d.merge(ppl[['research_id','age','sexc']],on='research_id',how='left'); COV+=['age','C(sexc)']
    except Exception: pass
    try:
        SMOKE=['F17','Z72.0','305.1','V15.82']; lk=" OR ".join([f"c.concept_code LIKE '{p}%'" for p in SMOKE])
        smk=set(str(r.person_id) for r in bq.query(f"SELECT DISTINCT co.person_id FROM `{PROJ}.{DS}.condition_occurrence` co JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id WHERE c.vocabulary_id LIKE 'ICD%' AND ({lk})"))
        d['smoke']=d.research_id.isin(smk).astype(int); COV+=['smoke']
    except Exception: pass
    covs=(" + "+" + ".join(COV)) if COV else ""
    print(f"   cohort {len(d)} | covariates {len(COV)}")
    import statsmodels.formula.api as smf
    def logit(fm,dd,t):
        try: r=smf.logit(fm,data=dd,missing='drop').fit(disp=0); return round(float(np.exp(r.params[t])),3),float(r.pvalues[t])
        except Exception: return None,None
    # ---- grid: each exposure x {CYSTIC, FIBROTIC} ----
    grid=[]; EXPO_ORDER=['SAVI_STING','COPA_WD40','IRAK3_LoF','R85H','R85H_x_IRAK3LoF','R85H_x_AQ']
    for ename in EXPO_ORDER:
        d['E']=d.research_id.isin(expo[ename]).astype(int); ncar=int(d.E.sum())
        for pname in ['CYSTIC','FIBROTIC']:
            d['Y']=d.research_id.isin(PH[pname]).astype(int); nco=int((d.E*d.Y).sum())
            OR,p=logit(f'Y ~ E{covs}',d,'E')
            grid.append({'exposure':ename,'phenotype':pname,'carriers':ncar,'cases_in_carriers':nco,'OR':OR,'p':p,'suppress':nco<20})
    ps=[x['p'] for x in grid if x['p'] is not None]
    fdr=dict(zip([i for i,x in enumerate(grid) if x['p'] is not None],multipletests(ps,method='fdr_bh')[1])) if ps else {}
    for i,x in enumerate(grid): x['fdr']=round(float(fdr[i]),4) if i in fdr else None; x['p']=round(x['p'],5) if x['p'] is not None else None
    S['grid']=grid
    print("\n== exposure x phenotype (variant-specific, sign-resolved) ==")
    print(f"{'exposure':18s} {'pheno':9s} {'carr':>6s} {'ca':>4s} {'OR':>7s} {'p':>9s} {'FDR':>7s}")
    for x in grid:
        cc='(ca<20 INTERNAL)' if x['suppress'] else ''
        print(f"{x['exposure']:18s} {x['phenotype']:9s} {x['carriers']:>6d} {x['cases_in_carriers']:>4d} {str(x['OR']):>7s} {str(x['p']):>9s} {str(x['fdr']):>7s} {cc}")
    # ---- the sign-flip readout (the whole point) ----
    print("\n== SIGN-FLIP (same exposure: cystic vs fibrotic) ==")
    for ename in EXPO_ORDER:
        c=next((x for x in grid if x['exposure']==ename and x['phenotype']=='CYSTIC'),None)
        fb=next((x for x in grid if x['exposure']==ename and x['phenotype']=='FIBROTIC'),None)
        if c and fb and c['OR'] and fb['OR']:
            flip='FLIP' if (c['OR']-1)*(fb['OR']-1)<0 else 'same-dir'
            print(f"   {ename:18s} cystic OR={c['OR']:>6} / fibrotic OR={fb['OR']:>6}  [{flip}]")
    print("\n===== CYSTIC vs FIBROTIC, VARIANT-SPECIFIC (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
