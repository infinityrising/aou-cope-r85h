#!/usr/bin/env python
"""AoU v9 — IRAK3 TRUNCATING-LoF ONLY (LoF='HC': stop-gain / frameshift / splice). ZS: OA37 + relative NC48 carry a
TRUNCATION (IRAK3 p.L210*); a damaging-MISSENSE is NOT mechanistically equivalent for a LoF-tolerant brake gene, and
mixing it in dilutes the signal. So restrict to true LoF = the OA37-matched class. Re-runs the key analyses LoF-ONLY:
 (1) IRAK3-LoF -> ISG (RNA, covariate-adjusted) — does our strongest molecular signal hold for the truncation class?
 (2) R85H × IRAK3-LoF double-carriers -> who + phenotype (the truncation-matched OA37 analogs; subset of the earlier 4).
 (3) IRAK3-LoF -> asthma/ILD/cystic (full 535k) — re-confirm the population direction, LoF-only.
Individual-level part (2) is <20 -> INTERNAL ONLY. STANDARD app. Ends 'run complete' / 'run failed'.
"""
import os, gzip, io, ast, json, subprocess
import numpy as np, pandas as pd
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); SNP=f"{MNT}/v9/wgs/short_read/snpindel"
RNADIR=f"{MNT}/v9/multiomics/rnaseq"; ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
TPM=f"{RNADIR}/rsem/aou_rnaseq_20260415.rsem_genes_tpm.txt.gz"; CELLS=os.path.expanduser("~/cell_fractions_v9.csv")
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1); R85H_VID='19-18911007-C-T'; IRAK3=('12',66183995,66259622)
IFN={'IFI27':'ENSG00000165949','IFI44L':'ENSG00000137959','IFIT1':'ENSG00000185745','ISG15':'ENSG00000187608','RSAD2':'ENSG00000134321','SIGLEC1':'ENSG00000088827'}
PHE={'asthma':['J45','493'],'ILD':['J84','515','516','J98.4'],'cystic_ILD':['J84','515','516','J98.4','J47','494']}
PANEL={'asthma':['J45','493'],'ILD':['J84','515','516','J98.4'],'bronchiectasis':['J47','494'],'RA':['M05','M06','714'],'vasculitis':['M31','M30','I77.6','446'],'SLE':['M32','710.0'],'sarcoid':['D86']}
S={}
def sh(c): return subprocess.run(['bash','-lc',c],capture_output=True,text=True).stdout
def fnum(x):
    try: return float(x)
    except (TypeError,ValueError): return 0.0
try:
    import pysam
    with gzip.open(VAT,'rt') as fh: COLS=fh.readline().rstrip("\n").split("\t")
    hasaa='aa_change' in COLS; IX={c:COLS.index(c) for c in ['vid','gene_symbol','consequence','gvs_afr_af','LoF']+(['aa_change'] if hasaa else [])}
    tbx=pysam.TabixFile(VAT); lof={}
    for line in tbx.fetch('chr'+IRAK3[0],IRAK3[1],IRAK3[2]):
        f=line.split("\t")
        if f[IX['gene_symbol']]!='IRAK3' or f[IX['LoF']]!='HC': continue
        if fnum(f[IX['gvs_afr_af']])>=0.001: continue
        lof[f[IX['vid']]]={'aa':f[IX['aa_change']] if hasaa else '','cons':f[IX['consequence']]}
    S['IRAK3_truncating_LoF_vids']=len(lof); print(f"IRAK3 TRUNCATING-LoF (HC) vids {len(lof)}")
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    def carr(vids):
        m={}
        for i in range(0,len(vids),900):
            inl=",".join("'"+v+"'" for v in vids[i:i+900])
            for r in bq.query(f"SELECT vid, e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid IN ({inl})"): m.setdefault(str(r.pid),[]).append(r.vid)
        return m
    lofcar=carr(list(lof.keys())); lofset=set(lofcar); r85h=set(carr([R85H_VID]).keys())
    print(f"IRAK3-LoF carriers {len(lofset)} | R85H {len(r85h)}")
    anc=pd.read_csv(ANC,sep="\t"); anc['research_id']=anc.research_id.astype(str); amap=dict(zip(anc.research_id,anc.ancestry_pred))
    def pcs_of(df):
        try:
            P=anc['pca_features'].apply(lambda x: ast.literal_eval(x) if isinstance(x,str) else (x if isinstance(x,list) else []))
            k=min(min(map(len,P[P.map(len)>0])) if (P.map(len)>0).any() else 0,16)
            p=pd.DataFrame([row[:k] for row in P],columns=[f'PC{i+1}' for i in range(k)]); p['research_id']=anc.research_id.values
            return df.merge(p,on='research_id',how='left'),[f'PC{i+1}' for i in range(k)]
        except Exception: return df,[]
    import statsmodels.formula.api as smf
    # ---- (1) ISG ----
    print("== (1) IRAK3-LoF -> ISG (RNA, covariate-adjusted) ==")
    raw=sh(f'zcat "{TPM}" 2>/dev/null | awk -F"\\t" \'NR==1 || $1 ~ /{"|".join(IFN.values())}/\'')
    m=pd.read_csv(io.StringIO(raw),sep="\t"); m['ensg']=m.gene_id.str.split('.').str[0]
    m=m[m.ensg.isin(IFN.values())].set_index('ensg').drop(columns=['gene_id','transcript_id(s)'])
    ex=m.T.astype(float); ex.index=ex.index.astype(str); lg=np.log2(ex+1); isg=((lg-lg.mean())/lg.std()).mean(axis=1)
    d=pd.DataFrame({'research_id':isg.index,'ifn':isg.values}).merge(anc[['research_id','ancestry_pred']],on='research_id',how='left')
    d['IRAK3_lof']=d.research_id.isin(lofset).astype(int); d,pc=pcs_of(d); COV=pc[:]
    try:
        cf=pd.read_csv(CELLS); cf['research_id']=cf.research_id.astype(str); ctc=[c for c in cf.columns if c!='research_id']; d=d.merge(cf,on='research_id',how='left'); COV+=ctc
    except Exception: pass
    covs=(" + "+" + ".join(COV)) if COV else ""
    def ols(f,dd,t):
        try: r=smf.ols(f,data=dd,missing='drop').fit(); return {'beta':round(float(r.params[t]),4),'p':round(float(r.pvalues[t]),4),'n':int(r.nobs),'n_carr':int(dd[t].sum())}
        except Exception as e: return {'error':str(e)[:60]}
    S['ISG_IRAK3lof_crude']=ols('ifn ~ IRAK3_lof + C(ancestry_pred)',d,'IRAK3_lof'); S['ISG_IRAK3lof_adj']=ols(f'ifn ~ IRAK3_lof + C(ancestry_pred){covs}',d,'IRAK3_lof')
    print(f"   ISG ~ IRAK3-LoF: crude(+anc) {S['ISG_IRAK3lof_crude']} | +PC+cell {S['ISG_IRAK3lof_adj']}  (vs combined-damaging was +0.31 p0.006, n56)")
    # ---- (2) doubles ----
    print("== (2) R85H × IRAK3-LoF double-carriers (truncation-matched to OA37) ==")
    doubles=sorted(lofset & r85h); S['n_doubles_LoF']=len(doubles); print(f"   R85H × IRAK3-LoF doubles: {len(doubles)}")
    if doubles:
        ids=",".join(doubles)
        demo={str(r.person_id):{'age':r.age,'sex':int(r.sex)} for r in bq.query(f"SELECT person_id, DATE_DIFF(DATE '2025-01-01', DATE(birth_datetime), YEAR) age, sex_at_birth_concept_id sex FROM `{PROJ}.{DS}.person` WHERE person_id IN ({ids})")}
        allc=[c for v in PANEL.values() for c in v]; lk=" OR ".join([f"c.concept_code LIKE '{p}%'" for p in allc]); phe={x:set() for x in doubles}
        for r in bq.query(f"SELECT DISTINCT co.person_id, c.concept_code FROM `{PROJ}.{DS}.condition_occurrence` co JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id WHERE co.person_id IN ({ids}) AND c.vocabulary_id LIKE 'ICD%' AND ({lk})"):
            for lab,codes in PANEL.items():
                if any(str(r.concept_code).startswith(p) for p in codes): phe[str(r.person_id)].add(lab)
        rows=[]
        for x in doubles:
            rec={'id':x,'ancestry':amap.get(x),'age':demo.get(x,{}).get('age'),'sex':demo.get(x,{}).get('sex'),'IRAK3_LoF':[(v,lof[v]['aa']) for v in lofcar[x]],'disease':sorted(phe[x]) if phe[x] else 'NONE'}
            rows.append(rec); print("   ",rec)
        S['doubles_LoF']=rows; S['n_affected']=sum(1 for r in rows if r['disease']!='NONE'); print(f"   -> {S['n_affected']}/{len(doubles)} affected")
    # ---- (3) full-cohort lung ----
    print("== (3) IRAK3-LoF -> lung (full 535k) ==")
    full=anc[['research_id','ancestry_pred']].copy(); full['IRAK3_lof']=full.research_id.isin(lofset).astype(int); full,fpc=pcs_of(full); fcovs=(" + "+" + ".join(fpc)) if fpc else ""
    def cases(codes):
        lk=" OR ".join([f"c.concept_code LIKE '{p}%'" for p in codes])
        return set(str(r.person_id) for r in bq.query(f"SELECT DISTINCT co.person_id FROM `{PROJ}.{DS}.condition_occurrence` co JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id WHERE c.vocabulary_id LIKE 'ICD%' AND ({lk})"))
    def logit(f,dd,t):
        try: r=smf.logit(f,data=dd,missing='drop').fit(disp=0); return {'OR':round(float(np.exp(r.params[t])),3),'p':round(float(r.pvalues[t]),4),'n':int(r.nobs)}
        except Exception as e: return {'error':str(e)[:50]}
    for ph,codes in PHE.items():
        full[ph]=full.research_id.isin(cases(codes)).astype(int); ncc=int((full.IRAK3_lof*full[ph]).sum())
        S[f'{ph}_IRAK3lof']=logit(f'{ph} ~ IRAK3_lof{fcovs}',full,'IRAK3_lof'); print(f"   {ph}: LoF carrier-cases {ncc} -> {S[f'{ph}_IRAK3lof']}")
    print("\n===== IRAK3-LoF-ONLY (paste back; part 2 INTERNAL) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
