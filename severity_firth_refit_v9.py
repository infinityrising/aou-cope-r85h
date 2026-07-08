#!/usr/bin/env python
"""AoU v9 -- FIRTH penalized-logistic refit of the small-cell / quasi-separation ORs (appraisal fix #2). The reported
vasculitis OR 11-12 (n=6) and COPA_WD40 OR 99-224 (n=5) are maximum-likelihood separation artifacts with meaningless
CIs. Refit with Firth's bias-reduced logistic (penalized-likelihood-ratio p, valid under separation) + a crude Fisher's
exact 2x2 as an assumption-free anchor. Also adds APOL1-G1 (S342G) as a covariate for the renal-containing systemic
outcomes (appraisal: APOL1 confound on organ-failure endpoints in an AFR-enriched exposure). Self-contained Firth (numpy).
Organ-appropriate severity: lung -> respiratory-hard; systemic -> organ-failure-hard. Standard app. Ends 'run complete'/'run failed'.
"""
import os, gzip, ast, json, re
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, chi2
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); SNP=f"{MNT}/v9/wgs/short_read/snpindel"
ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
R85H_VID='19-18911007-C-T'
PANEL={'bronchiectasis':['J47','494'],'infl_arthritis':['M05','M06','M08','714'],'vasculitis_broad':['M30','M31']}
RESP_HARD=['J96','518.81','518.83','518.84','Z99.81','V46.2','Z94.2','J95.85','Z99.11']              # lung severity (no renal)
ORGAN_HARD=RESP_HARD+['Z99.2','N18.6','V42.0','N17']                                                  # systemic severity (renal appropriate)
SEVMAP={'bronchiectasis':RESP_HARD,'infl_arthritis':ORGAN_HARD,'vasculitis_broad':ORGAN_HARD}
COPA=('1',160288580,160343566); COPA_RES={230,233,236,241,242,243}
APOL1=('22',36253071,36267530); RESID=re.compile(r'p\.\(?[A-Za-z]{3}(\d+)')
MINCODES=2
def fnum(x):
    try: return float(x)
    except (TypeError,ValueError): return 0.0
def firth(X,y,max_iter=200,tol=1e-8):
    n,p=X.shape; b=np.zeros(p)
    for _ in range(max_iter):
        eta=X@b; pi=np.clip(1/(1+np.exp(-eta)),1e-12,1-1e-12); W=pi*(1-pi)
        I=(X.T*W)@X
        try: Iinv=np.linalg.inv(I)
        except np.linalg.LinAlgError: Iinv=np.linalg.pinv(I)
        h=W*np.einsum('ij,jk,ik->i',X,Iinv,X)
        U=X.T@(y-pi+h*(0.5-pi)); step=Iinv@U; b=b+step
        if np.max(np.abs(step))<tol: break
    eta=X@b; pi=np.clip(1/(1+np.exp(-eta)),1e-12,1-1e-12)
    ll=np.sum(y*np.log(pi)+(1-y)*np.log(1-pi)); _,logdet=np.linalg.slogdet(I); pll=ll+0.5*logdet
    se=np.sqrt(np.clip(np.diag(Iinv),0,None)); return b,se,pll
def firth_test(df,exposure,covcols):
    d=df.dropna(subset=[exposure]+covcols).copy()
    cols=[exposure]+[c for c in covcols if d[c].nunique()>1]
    Xf=np.column_stack([np.ones(len(d))]+[d[c].values.astype(float) for c in cols]); y=d[exposure+'__y'].values.astype(float) if False else d['__y'].values.astype(float)
    bf,sef,pllf=firth(Xf,y)
    Xn=np.column_stack([np.ones(len(d))]+[d[c].values.astype(float) for c in cols if c!=exposure]); _,_,plln=firth(Xn,y)
    lr=2*(pllf-plln); p=float(chi2.sf(lr,1)); j=1
    OR=float(np.exp(bf[j])); lo=float(np.exp(bf[j]-1.96*sef[j])); hi=float(np.exp(bf[j]+1.96*sef[j]))
    return {'OR':round(OR,3),'CI':[round(lo,3),round(hi,3)],'p_LRT':round(p,4),'n':int(len(d))}
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
    with gzip.open(VAT,'rt') as fh: COLS=fh.readline().rstrip("\n").split("\t")
    hasaa='aa_change' in COLS; IX={c:COLS.index(c) for c in ['vid','gene_symbol','consequence','gvs_afr_af','revel','LoF']+(['aa_change'] if hasaa else [])}
    tbx=pysam.TabixFile(VAT)
    copa_v=[]; apol1_v=[]
    for line in tbx.fetch('chr'+COPA[0],COPA[1],COPA[2]):
        f=line.split("\t")
        if f[IX['gene_symbol']]=='COPA' and 'missense' in f[IX['consequence']] and fnum(f[IX['gvs_afr_af']])<0.01:
            mm=RESID.search(f[IX['aa_change']] if hasaa else '')
            if mm and int(mm.group(1)) in COPA_RES: copa_v.append(f[IX['vid']])
    for line in tbx.fetch('chr'+APOL1[0],APOL1[1],APOL1[2]):
        f=line.split("\t")
        if f[IX['gene_symbol']]=='APOL1' and hasaa:
            mm=RESID.search(f[IX['aa_change']])
            if mm and int(mm.group(1)) in {342,384}: apol1_v.append(f[IX['vid']])   # G1 (S342G / I384M)
    C_r85h=carr([R85H_VID]); C_copa=carr(copa_v); C_apol1=carr(apol1_v)
    print(f"R85H {len(C_r85h)} | COPA-WD40 vids {len(copa_v)}->carr {len(C_copa)} | APOL1-G1 vids {len(apol1_v)}->carr {len(C_apol1)}")
    hard={ph:anyset(SEVMAP[ph]) for ph in PANEL}; cases={ph:caseset(cc) for ph,cc in PANEL.items()}
    # cohort covariates
    anc=pd.read_csv(ANC,sep="\t"); anc['research_id']=anc.research_id.astype(str)
    base=anc[['research_id','ancestry_pred']].copy()
    try:
        P=anc['pca_features'].apply(lambda x: ast.literal_eval(x) if isinstance(x,str) else (x if isinstance(x,list) else []))
        k=min(min(map(len,P[P.map(len)>0])) if (P.map(len)>0).any() else 0,3)
        pcf=pd.DataFrame([row[:k] for row in P],columns=[f'PC{i+1}' for i in range(k)]); pcf['research_id']=anc.research_id.values
        base=base.merge(pcf,on='research_id',how='left'); PCS=[f'PC{i+1}' for i in range(k)]
    except Exception: PCS=[]
    try:
        ppl=bq.query(f"SELECT person_id, DATE_DIFF(DATE '2025-01-01', DATE(birth_datetime), YEAR) age, CAST(sex_at_birth_concept_id AS STRING) sexc FROM `{PROJ}.{DS}.person`").to_dataframe()
        ppl['research_id']=ppl.person_id.astype(str); base=base.merge(ppl[['research_id','age','sexc']],on='research_id',how='left')
        base['agez']=(base.age-base.age.mean())/base.age.std(); base['sex_m']=(base.sexc.astype(str)==base.sexc.astype(str).mode().iloc[0]).astype(float)
    except Exception: base['agez']=0.0; base['sex_m']=0.0
    base['apol1']=base.research_id.isin(C_apol1).astype(float)
    covbase=['agez','sex_m']+PCS
    # ---- (1) R85H within-case SEVERITY, Firth + Fisher ----
    print("\n== (1) R85H -> SEVERE within cases: Firth (penalized-LRT) + Fisher exact (crude 2x2) ==")
    S['severity']={}
    for ph in PANEL:
        sub=base[base.research_id.isin(cases[ph])].copy(); sub['__y']=sub.research_id.isin(hard[ph]).astype(float); sub['R85H']=sub.research_id.isin(C_r85h).astype(float)
        a=int(((sub.R85H==1)&(sub.__y==1)).sum()); b=int(((sub.R85H==1)&(sub.__y==0)).sum()); c=int(((sub.R85H==0)&(sub.__y==1)).sum()); dd=int(((sub.R85H==0)&(sub.__y==0)).sum())
        fo,fp=fisher_exact([[a,b],[c,dd]])
        cov=covbase+(['apol1'] if ph!='bronchiectasis' else [])
        ft=firth_test(sub,'R85H',cov)
        S['severity'][ph]={'2x2':[a,b,c,dd],'fisher_OR':round(float(fo),3),'fisher_p':round(float(fp),4),'firth':ft}
        print(f"   {ph:16s} 2x2[sev+/sev-|R85H+ {a}/{b}, R85H- {c}/{dd}] Fisher OR={round(float(fo),2)} p={round(float(fp),3)} | Firth OR={ft['OR']} CI{ft['CI']} p_LRT={ft['p_LRT']} (n={ft['n']})")
    # ---- (2) COPA-WD40 INCIDENCE, Firth (fix the OR 99-224 separation) ----
    print("\n== (2) COPA-WD40 incidence: Firth (replaces the ML separation OR 99-224) ==")
    S['copa_incidence']={}
    for ph,cc in [('bronchiectasis',cases['bronchiectasis']),('fibrotic_ILD',caseset(['J84.1','515','516.3']))]:
        sub=base.copy(); sub['__y']=sub.research_id.isin(cc).astype(float); sub['COPA']=sub.research_id.isin(C_copa).astype(float)
        a=int(((sub.COPA==1)&(sub.__y==1)).sum()); nC=int(sub.COPA.sum())
        ft=firth_test(sub,'COPA',covbase)
        S['copa_incidence'][ph]={'copa_carriers':nC,'copa_cases':a,'firth':ft}
        print(f"   COPA-WD40 -> {ph:14s}: carriers {nC}, cases {a} | Firth OR={ft['OR']} CI{ft['CI']} p_LRT={ft['p_LRT']}")
    print("\n== READ: Firth CIs are the honest ones; if they span 1, the ML point estimate (11/224) was an artifact. ==")
    print("\n===== FIRTH REFIT (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1400:])
