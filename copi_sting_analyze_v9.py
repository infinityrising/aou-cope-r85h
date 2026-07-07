#!/usr/bin/env python
"""AoU v9 — PHASE 2 (schema-robust + covariate-adjusted): STING phased haplotypes -> ISG stratified by COPI-mut × STING.
Auto-loads STING CSV (prefers ~/copi_sting_pav_v9.csv, else ~/sting_only_pav_v9.csv). R85H + damaging-COPI from BigQuery.
Outcome = per-person 6-gene ISG z (Rice/Crow). MULTIVARIABLE OLS with covariates loaded robustly (each optional, degrades
gracefully): fine-scale ancestry PCs (founder-tilt confound), age, sex, and cell-composition scores (~/cell_fractions_v9.csv
from cell_deconv_v9.py — bulk-blood ISG confound). Models: LADDER ISG~HAQ/AQ/R220H (+cov); FLAGSHIP ISG~R85H*HAQ+cov
(HAQ protective => R85H:HAQ<0); BROADER ISG~COPImut*HAQ+cov. Also a 5-gene ISG (drop monocyte-marker SIGLEC1) sensitivity.
STANDARD app. Ends 'run complete' / 'run failed'.
"""
import os, subprocess, io, json, gzip, ast
from collections import Counter
import numpy as np, pandas as pd
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9")
RNADIR=f"{MNT}/v9/multiomics/rnaseq"; SNP=f"{MNT}/v9/wgs/short_read/snpindel"
TPM=f"{RNADIR}/rsem/aou_rnaseq_20260415.rsem_genes_tpm.txt.gz"; ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"
VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1); R85H_VID='19-18911007-C-T'
CSV_FULL=os.path.expanduser("~/copi_sting_pav_v9.csv"); CSV_STING=os.path.expanduser("~/sting_only_pav_v9.csv")
CELLS=os.path.expanduser("~/cell_fractions_v9.csv")
COPI={'COPA':('1',160288580,160343566),'COPB1':('11',14443357,14500027),'COPB2':('3',139353942,139389736),
      'COPG1':('3',129249575,129278068),'COPG2':('7',130505553,130668755),'COPZ1':('12',54301202,54351848),
      'COPZ2':('17',48026142,48038030),'COPE':('19',18899511,18919407),'ARCN1':('11',118572384,118603033)}
IFN={'IFI27':'ENSG00000165949','IFI44L':'ENSG00000137959','IFIT1':'ENSG00000185745','ISG15':'ENSG00000187608','RSAD2':'ENSG00000134321','SIGLEC1':'ENSG00000088827'}
S={}
def sh(c): return subprocess.run(['bash','-lc',c],capture_output=True,text=True).stdout
def fnum(x):
    try: return float(x)
    except (TypeError,ValueError): return 0.0
try:
    CSV=CSV_FULL if (os.path.exists(CSV_FULL) and os.path.getsize(CSV_FULL)>0) else CSV_STING
    g=pd.read_csv(CSV); g['research_id']=g.research_id.astype(str)
    for c in ['HAQ_d','AQ_d','R220H_d','R85H_d','copi_dmg']:
        if c in g.columns: g[c]=pd.to_numeric(g[c],errors='coerce').fillna(0).astype(int)
    has_copi='copi_dmg' in g.columns
    S['csv']=os.path.basename(CSV); S['n_loaded']=len(g); S['hap_freq']=dict(Counter(list(g.hap1)+list(g.hap2)))
    print(f"loaded {len(g)} from {S['csv']} | hap freq {S['hap_freq']} | phased-COPI col {has_copi}")
    print("== ISG (6-gene) + 5-gene (drop SIGLEC1) ==")
    raw=sh(f'zcat "{TPM}" 2>/dev/null | awk -F"\\t" \'NR==1 || $1 ~ /{"|".join(IFN.values())}/\'')
    m=pd.read_csv(io.StringIO(raw),sep="\t"); m['ensg']=m.gene_id.str.split('.').str[0]
    m=m[m.ensg.isin(IFN.values())].set_index('ensg').drop(columns=['gene_id','transcript_id(s)'])
    ex=m.T.astype(float); ex.index=ex.index.astype(str); lg=np.log2(ex+1); zz=(lg-lg.mean())/lg.std()
    isg6=zz.mean(axis=1); isg5=zz.drop(columns=[IFN['SIGLEC1']],errors='ignore').mean(axis=1)
    isgdf=pd.DataFrame({'research_id':isg6.index,'ifn':isg6.values,'ifn5':isg5.reindex(isg6.index).values})
    print("== R85H + damaging-COPI carriers (BigQuery) + ancestry ==")
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    def carriers(vids):
        out=set()
        for i in range(0,len(vids),900):
            inl=",".join("'"+v+"'" for v in vids[i:i+900])
            out|=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid IN ({inl})"))
        return out
    r85h_bq=carriers([R85H_VID])
    import pysam
    with gzip.open(VAT,'rt') as fh: COLS=fh.readline().rstrip("\n").split("\t")
    IX={c:COLS.index(c) for c in ['vid','gene_symbol','consequence','revel','LoF']}
    tbx=pysam.TabixFile(VAT); dmg=[]
    for gene,(ch,s,e) in COPI.items():
        for line in tbx.fetch(f'chr{ch}',s,e):
            f=line.split("\t")
            if f[IX['gene_symbol']]!=gene or f[IX['vid']]==R85H_VID: continue
            if f[IX['LoF']]=='HC' or ('missense' in f[IX['consequence']] and fnum(f[IX['revel']])>0.5): dmg.append(f[IX['vid']])
    copi_bq=carriers(sorted(set(dmg)))
    print(f"   R85H(BQ) {len(r85h_bq)} | damaging-COPI vids {len(set(dmg))} carriers {len(copi_bq)}")
    anc=pd.read_csv(ANC,sep="\t"); anc['research_id']=anc.research_id.astype(str)
    d=g.merge(isgdf,on='research_id',how='inner').merge(anc[['research_id','ancestry_pred']],on='research_id',how='left')
    d['R85H']=d.research_id.isin(r85h_bq).astype(int); d['HAQ']=(d.HAQ_d>0).astype(int)
    d['AQ']=(d.AQ_d>0).astype(int); d['HYPO']=((d.HAQ_d>0)|(d.AQ_d>0)).astype(int)   # AQ-strict = pre-reg PRIMARY; HYPO = either hypomorph (powered)
    d['COPImut']=(d.research_id.isin(copi_bq)|(d.R85H>0)).astype(int)
    # ---- covariates (each optional; degrade gracefully) ----
    COVcols=[]
    try:  # fine-scale ancestry PCs
        pcc=[c for c in anc.columns if c.lower().startswith('pc') and c[2:].isdigit()]
        if pcc: pcs=anc[['research_id']+pcc]; pcnames=pcc
        elif 'pca_features' in anc.columns:
            P=anc['pca_features'].apply(lambda x: ast.literal_eval(x) if isinstance(x,str) else (x if isinstance(x,list) else []))
            k=min(min(map(len,P[P.map(len)>0])) if (P.map(len)>0).any() else 0,16)
            pcs=pd.DataFrame([row[:k] for row in P],columns=[f'PC{i+1}' for i in range(k)]); pcs['research_id']=anc.research_id.values; pcnames=[f'PC{i+1}' for i in range(k)]
        else: pcs=None; pcnames=[]
        if pcs is not None and pcnames:
            d=d.merge(pcs,on='research_id',how='left'); COVcols+=pcnames; print(f"   +PCs {len(pcnames)}")
    except Exception as ex: print("   (PCs skipped:",str(ex)[:70],")")
    try:  # age + sex
        ppl=bq.query(f"SELECT person_id, DATE_DIFF(DATE '2025-01-01', DATE(birth_datetime), YEAR) age, CAST(sex_at_birth_concept_id AS STRING) sexc FROM `{PROJ}.{DS}.person`").to_dataframe()
        ppl['research_id']=ppl.person_id.astype(str)
        d=d.merge(ppl[['research_id','age','sexc']],on='research_id',how='left'); COVcols+=['age','C(sexc)']; print("   +age,sex")
    except Exception as ex: print("   (age/sex skipped:",str(ex)[:70],")")
    try:  # cell composition
        cf=pd.read_csv(CELLS); cf['research_id']=cf.research_id.astype(str); ctc=[c for c in cf.columns if c!='research_id']
        d=d.merge(cf,on='research_id',how='left'); COVcols+=ctc; print(f"   +cell {ctc}")
    except Exception as ex: print("   (cell comp skipped — run cell_deconv_v9.py:",str(ex)[:50],")")
    COV=(" + "+" + ".join(COVcols)) if COVcols else ""
    S['covariates']=COVcols; print(f"   COVARIATE SET: {COVcols}")
    a=d[d.ancestry_pred=='afr'].copy()
    S['n_ifn_join']=len(d); S['n_afr']=len(a); S['afr_R85H']=int(a.R85H.sum()); S['afr_HAQ']=int(a.HAQ.sum()); S['afr_AQ']=int(a.AQ.sum()); S['afr_HYPO']=int(a.HYPO.sum()); S['afr_COPImut']=int(a.COPImut.sum())
    if 'R85H_d' in g.columns: S['R85H_PAV_vs_BQ']={'pav':int((d.R85H_d>0).sum()),'bq':int(d.R85H.sum())}
    print(f"   ISG-joined {len(d)} | AFR {len(a)} | R85H {int(a.R85H.sum())} | HAQ {int(a.HAQ.sum())} | AQ {int(a.AQ.sum())} | HYPO {int(a.HYPO.sum())} | COPImut {int(a.COPImut.sum())}")
    import statsmodels.formula.api as smf
    def ols(f,dd,term):
        try: r=smf.ols(f,data=dd,missing='drop').fit(); return {'beta':round(float(r.params[term]),4),'p':round(float(r.pvalues[term]),4),'n':int(r.nobs)}
        except Exception as e: return {'error':str(e)[:70]}
    print("== LADDER: STING -> ISG, covariate-adjusted (AQ/HAQ hypomorph => NEGATIVE = validates) ==")
    for t in ['HAQ_d','AQ_d','R220H_d']:
        S[f'ladder_all_{t}']=ols(f'ifn ~ {t} + C(ancestry_pred){COV}',d,t); S[f'ladder_afr_{t}']=ols(f'ifn ~ {t}{COV}',a,t)
        print(f"   ifn~{t}: ALL {S[f'ladder_all_{t}']} | AFR {S[f'ladder_afr_{t}']}")
    print("== FLAGSHIP: R85H × modifier -> ISG (protective => interaction NEGATIVE). AQ-strict = PRE-REG PRIMARY ==")
    def cells(dd,mod,col='ifn'):
        out={}
        for r in (0,1):
            for h in (0,1):
                sub=dd[(dd.R85H==r)&(dd[mod]==h)][col].dropna()
                out[f'R85H{r}_{mod}{h}']={'mean':round(float(sub.mean()),3) if len(sub) else None,'n':int(len(sub))}
        return out
    for mod in ['HAQ','AQ','HYPO']:
        S[f'cells_afr_{mod}']=cells(a,mod); S[f'cells_all_{mod}']=cells(d,mod)
        S[f'inter_crude_afr_{mod}']=ols(f'ifn ~ R85H*{mod}',a,f'R85H:{mod}')
        S[f'inter_adj_afr_{mod}']=ols(f'ifn ~ R85H*{mod}{COV}',a,f'R85H:{mod}')
        S[f'inter_adj_all_{mod}']=ols(f'ifn ~ R85H*{mod} + C(ancestry_pred){COV}',d,f'R85H:{mod}')
        S[f'inter_adj5_afr_{mod}']=ols(f'ifn5 ~ R85H*{mod}{COV}',a,f'R85H:{mod}')
        print(f"   [{mod}] AFR 2x2 {S[f'cells_afr_{mod}']}")
        print(f"        AFR crude {S[f'inter_crude_afr_{mod}']} | AFR adj {S[f'inter_adj_afr_{mod}']} | ALL adj {S[f'inter_adj_all_{mod}']} | 5g {S[f'inter_adj5_afr_{mod}']}")
    print("== CORRECTED (common WT reference): mutually-exclusive STING classes + JOINT dosage model ==")
    S['n_compound_HAQ_AQ']=int(((d.HAQ_d>0)&(d.AQ_d>0)).sum()); S['n_compound_afr_R85Hp']=int(((a.R85H==1)&(a.HAQ_d>0)&(a.AQ_d>0)).sum())
    print(f"   compound HAQ+AQ persons (HAQ on one allele, AQ on the other): cohort {S['n_compound_HAQ_AQ']} | AFR R85H+ {S['n_compound_afr_R85Hp']}")
    def sting_class(r):
        h={r.hap1,r.hap2}
        if 'HAQ' in h and 'AQ' in h: return 'HAQ/AQ'
        if 'HAQ' in h: return 'HAQ'
        if 'AQ' in h: return 'AQ'
        return 'WT/other'
    rp=a[a.R85H==1].copy(); rp['sc']=rp.apply(sting_class,axis=1)
    S['R85Hp_by_sting_class']={c:{'n':int((rp.sc==c).sum()),'mean_ifn':round(float(rp[rp.sc==c].ifn.mean()),3) if (rp.sc==c).any() else None,'median':round(float(rp[rp.sc==c].ifn.median()),3) if (rp.sc==c).any() else None,'n_IFN_activated_gt1SD':int((rp[rp.sc==c].ifn>1).sum())} for c in ['WT/other','HAQ','AQ','HAQ/AQ']}
    print("   AFR R85H+ by CLEAN STING class (mean/median/penetrance=n>1SD, each vs WT):",S['R85Hp_by_sting_class'])
    for hap in ['AQ_d','HAQ_d']:   # 1 vs 2 copies (het vs hom dose)
        S[f'R85Hp_by_{hap}']={int(x):{'n':int((rp[hap]==x).sum()),'mean_ifn':round(float(rp[rp[hap]==x].ifn.mean()),3) if (rp[hap]==x).any() else None} for x in sorted(rp[hap].unique())}
        print(f"   R85H+ by {hap} (0/1/2 = WT/het/hom): {S[f'R85Hp_by_{hap}']}")
    for lab,form,dd in [('AFR crude','ifn ~ R85H*HAQ_d + R85H*AQ_d',a),('AFR adj',f'ifn ~ R85H*HAQ_d + R85H*AQ_d{COV}',a),('ALL adj',f'ifn ~ R85H*HAQ_d + R85H*AQ_d + C(ancestry_pred){COV}',d)]:
        S[f'joint_HAQ_{lab}']=ols(form,dd,'R85H:HAQ_d'); S[f'joint_AQ_{lab}']=ols(form,dd,'R85H:AQ_d')
        print(f"   [{lab}] R85H:HAQ_d {S[f'joint_HAQ_{lab}']} | R85H:AQ_d {S[f'joint_AQ_{lab}']}")
    print("== BROADER: any damaging COPI × HAQ -> ISG, adjusted ==")
    S['inter_COPImut_all_adj']=ols(f'ifn ~ COPImut*HAQ + C(ancestry_pred){COV}',d,'COPImut:HAQ')
    print("   COPImut:HAQ (ALL adj):",S['inter_COPImut_all_adj'])
    print("\n===== COPI×STING×ISG (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except FileNotFoundError:
    print("run failed: no STING CSV — run sting_only_extract_v9.py or copi_sting_extract_v9.py first")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1400:])
