#!/usr/bin/env python
"""AoU v9 -- THE decisive missing experiment (all 5 lit agents converge on it): in whole-blood RNA-seq, does R85H /
IRAK3-LoF / STING-AQ shift a TYPE-I IFN vs TYPE-II (IFN-gamma) vs NF-kB vs NEUTROPHIL transcriptional module? Beyond the
6-gene type-I ISG (runs 022-027), this adjudicates the literature's IRAK-M question: IRAK-M loss is documented to raise
NF-kB / IFN-gamma (type II), NOT necessarily type-I IFN -- so separating the modules tells us which axis R85H+2nd-hit
actually engages. Per-sample module score = mean z of log2(TPM+1) over detected genes; module ~ exposure + 16PC + age +
sex + cell-composition (CLR fractions). Reports per-module gene coverage (so any ENSG miss is visible). Standard app.
Ends 'run complete' / 'run failed'.
"""
import os, io, gzip, subprocess, json
import numpy as np, pandas as pd
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); RNADIR=f"{MNT}/v9/multiomics/rnaseq"; SNP=f"{MNT}/v9/wgs/short_read/snpindel"
TPM=f"{RNADIR}/rsem/aou_rnaseq_20260415.rsem_genes_tpm.txt.gz"; ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CELL=os.path.expanduser("~/cell_fractions_v9.csv")
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
R85H_VID='19-18911007-C-T'; STINGV={'R71H':'5-139481493-C-T','G230A':'5-139478340-C-G','R293Q':'5-139477397-C-T'}
MODULES={
 'typeI_ISG':{'IFI27':'ENSG00000165949','IFI44L':'ENSG00000137959','IFIT1':'ENSG00000185745','ISG15':'ENSG00000187608','RSAD2':'ENSG00000134321','SIGLEC1':'ENSG00000088827'},
 'typeII_IFNg':{'CXCL9':'ENSG00000138755','CXCL10':'ENSG00000169245','CXCL11':'ENSG00000169248','GBP1':'ENSG00000117228','GBP5':'ENSG00000154451','STAT1':'ENSG00000115415','IRF1':'ENSG00000125347','IDO1':'ENSG00000131203'},
 'NFkB':{'TNF':'ENSG00000232810','IL6':'ENSG00000136244','IL1B':'ENSG00000125538','NFKBIA':'ENSG00000100906','CXCL8':'ENSG00000169429','CCL2':'ENSG00000108691','NFKB1':'ENSG00000109320','TNFAIP3':'ENSG00000118503'},
 'neutrophil':{'ELANE':'ENSG00000197561','MPO':'ENSG00000005381','PRTN3':'ENSG00000196415','S100A8':'ENSG00000143546','S100A9':'ENSG00000163220','FCGR3B':'ENSG00000162747'},
}
ALL_ENSG={e for m in MODULES.values() for e in m.values()}
def sh(c): return subprocess.run(['bash','-lc',c],capture_output=True,text=True).stdout
def fnum(x):
    try: return float(x)
    except (TypeError,ValueError): return 0.0
S={}
try:
    print("== extract module genes from RSEM TPM ==")
    pat="|".join(ALL_ENSG)
    raw=sh(f'zcat "{TPM}" 2>/dev/null | awk -F"\\t" \'NR==1 || $1 ~ /{pat}/\'')
    m=pd.read_csv(io.StringIO(raw),sep="\t"); m['ensg']=m.gene_id.str.split('.').str[0]
    m=m[m.ensg.isin(ALL_ENSG)].drop_duplicates('ensg').set_index('ensg')
    dropcols=[c for c in ['gene_id','transcript_id(s)'] if c in m.columns]; ex=m.drop(columns=dropcols).T.astype(float)
    ex.index=ex.index.astype(str); lg=np.log2(ex+1); z=(lg-lg.mean())/lg.std()
    found=set(m.index); cov={}
    d=pd.DataFrame(index=z.index)
    for mod,genes in MODULES.items():
        got=[g for g in genes.values() if g in found]; cov[mod]=f"{len(got)}/{len(genes)}"
        d[mod]=z[got].mean(axis=1) if got else np.nan
    S['module_coverage']=cov; print("   module gene coverage:",cov)
    d['research_id']=d.index.astype(str); rna_ids=set(d.research_id)
    print(f"   RNA samples {len(d)}")
    # ---- exposures (carriers ∩ RNA) ----
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    def carr(vids):
        out=set()
        for i in range(0,len(vids),900):
            if not vids[i:i+900]: continue
            inl=",".join("'"+v+"'" for v in vids[i:i+900])
            out|=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid IN ({inl})"))
        return out
    import pysam
    with gzip.open(VAT,'rt') as fh: COLS=fh.readline().rstrip("\n").split("\t")
    IX={c:COLS.index(c) for c in ['vid','gene_symbol','gvs_afr_af','LoF']}
    tbx=pysam.TabixFile(VAT); irak3=[]
    for line in tbx.fetch('chr12',66183995,66259622):
        f=line.split("\t")
        if f[IX['gene_symbol']]=='IRAK3' and f[IX['LoF']]=='HC' and fnum(f[IX['gvs_afr_af']])<0.001: irak3.append(f[IX['vid']])
    C_r85h=carr([R85H_VID]); C_irak3=carr(irak3); st={k:carr([v]) for k,v in STINGV.items()}
    C_aq=(st['G230A']&st['R293Q'])-st['R71H']
    d['R85H']=d.research_id.isin(C_r85h).astype(int); d['IRAK3_LoF']=d.research_id.isin(C_irak3).astype(int)
    d['AQ']=d.research_id.isin(C_aq).astype(int); d['R85H_x_IRAK3']=(d.R85H&d.IRAK3_LoF).astype(int); d['R85H_x_AQ']=(d.R85H&d.AQ).astype(int)
    print(f"   RNA carriers: R85H {int(d.R85H.sum())} | IRAK3-LoF {int(d.IRAK3_LoF.sum())} | AQ {int(d.AQ.sum())} | R85H&IRAK3 {int(d.R85H_x_IRAK3.sum())} | R85H&AQ {int(d.R85H_x_AQ.sum())}")
    # ---- covariates: PCs, age, sex, cell fractions ----
    anc=pd.read_csv(ANC,sep="\t"); anc['research_id']=anc.research_id.astype(str)
    COV=[]
    try:
        import ast
        P=anc['pca_features'].apply(lambda x: ast.literal_eval(x) if isinstance(x,str) else (x if isinstance(x,list) else []))
        k=min(min(map(len,P[P.map(len)>0])) if (P.map(len)>0).any() else 0,16)
        pcf=pd.DataFrame([row[:k] for row in P],columns=[f'PC{i+1}' for i in range(k)]); pcf['research_id']=anc.research_id.values
        d=d.merge(pcf,on='research_id',how='left'); COV+=[f'PC{i+1}' for i in range(k)]
    except Exception: pass
    try:
        ppl=bq.query(f"SELECT person_id, DATE_DIFF(DATE '2025-01-01', DATE(birth_datetime), YEAR) age, CAST(sex_at_birth_concept_id AS STRING) sexc FROM `{PROJ}.{DS}.person`").to_dataframe()
        ppl['research_id']=ppl.person_id.astype(str); d=d.merge(ppl[['research_id','age','sexc']],on='research_id',how='left'); COV+=['age','C(sexc)']
    except Exception: pass
    try:
        cf=pd.read_csv(CELL); cf['research_id']=cf.research_id.astype(str); cc=[c for c in cf.columns if c!='research_id']
        d=d.merge(cf,on='research_id',how='left'); COV+=cc; print(f"   cell-composition covariates: {cc}")
    except Exception: print("   (cell_fractions_v9.csv not found — running without cell adjustment)")
    covs=(" + "+" + ".join(COV)) if COV else ""
    import statsmodels.formula.api as smf
    def ols(mod,ex):
        sub=d.dropna(subset=[mod])
        try:
            r=smf.ols(f'{mod} ~ {ex}{covs}',data=sub,missing='drop').fit()
            return round(float(r.params[ex]),3),round(float(r.pvalues[ex]),4),int(r.nobs),int(sub[ex].sum())
        except Exception: return None,None,None,int(sub[ex].sum())
    print("\n== module ~ exposure (beta, p, n, n_carrier) -- covariate + cell-composition adjusted ==")
    print(f"{'exposure':14s} {'module':12s} {'beta':>7} {'p':>8} {'n':>6} {'carr':>5}")
    res=[]
    for ex in ['R85H','IRAK3_LoF','AQ','R85H_x_IRAK3','R85H_x_AQ']:
        for mod in MODULES:
            b,p,n,nc=ols(mod,ex); res.append({'exposure':ex,'module':mod,'beta':b,'p':p,'n':n,'n_carrier':nc})
            print(f"{ex:14s} {mod:12s} {str(b):>7} {str(p):>8} {str(n):>6} {nc:>5}")
    S['results']=res
    print("\n== READ: IRAK3-LoF -> which module? (lit predicts NF-kB/IFN-gamma > type-I). R85H alone -> expect ~flat (benign). ==")
    print("\n===== IFN MODULES RNA (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
