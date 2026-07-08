#!/usr/bin/env python
"""AoU v9 -- STING-CLASS-RESOLVED IFN modules. Fixes the reference-contamination in runs 044/050: 'AQ' was tested vs
EVERYONE (a reference that CONTAINS HAQ [hypomorphic/protective] + WT), and HAQ was never tested. Proper design:
mutually-exclusive STING class per person -- HAQ (R71H+G230A+R293Q), AQ (G230A+R293Q, R71-WT), or WT/neither (R232 = the
clean reference) -- plus a FORMAL interaction model  module ~ R85H*AQ + R85H*HAQ  (WT the implicit reference).
Predictions of the conditional-modifier thesis: HAQ main effect <= 0 (dampened); R85H×AQ > 0 (type-I IFN unmasked on the
gain background); R85H×HAQ ~ 0 / protective (HAQ = the built-in negative control via R71H). Whole-blood RNA; reports the
6-way stratified type-I ISG means (WT/AQ/HAQ × R85H-/+). NB STING class is unphased (person-level cis proxy). Standard app.
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
}
ALL_ENSG={e for m in MODULES.values() for e in m.values()}
def sh(c): return subprocess.run(['bash','-lc',c],capture_output=True,text=True).stdout
def fnum(x):
    try: return float(x)
    except (TypeError,ValueError): return 0.0
S={}
try:
    pat="|".join(ALL_ENSG); raw=sh(f'zcat "{TPM}" 2>/dev/null | awk -F"\\t" \'NR==1 || $1 ~ /{pat}/\'')
    m=pd.read_csv(io.StringIO(raw),sep="\t"); m['ensg']=m.gene_id.str.split('.').str[0]
    m=m[m.ensg.isin(ALL_ENSG)].drop_duplicates('ensg').set_index('ensg')
    dropc=[c for c in ['gene_id','transcript_id(s)'] if c in m.columns]; ex=m.drop(columns=dropc).T.astype(float)
    ex.index=ex.index.astype(str); lg=np.log2(ex+1); z=(lg-lg.mean())/lg.std()
    found=set(m.index); d=pd.DataFrame(index=z.index)
    for mod,genes in MODULES.items():
        got=[g for g in genes.values() if g in found]; d[mod]=z[got].mean(axis=1) if got else np.nan
    d['research_id']=d.index.astype(str)
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    def carr(vids):
        out=set()
        for v in vids: out|=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid='{v}'"))
        return out
    import pysam
    with gzip.open(VAT,'rt') as fh: COLS=fh.readline().rstrip("\n").split("\t")
    IX={c:COLS.index(c) for c in ['vid','gene_symbol','gvs_afr_af','LoF']}
    tbx=pysam.TabixFile(VAT); irak3=[]
    for line in tbx.fetch('chr12',66183995,66259622):
        f=line.split("\t")
        if f[IX['gene_symbol']]=='IRAK3' and f[IX['LoF']]=='HC' and fnum(f[IX['gvs_afr_af']])<0.001: irak3.append(f[IX['vid']])
    C_r85h=carr([R85H_VID]); C_irak3=carr(irak3); st={k:carr([v]) for k,v in STINGV.items()}
    C_aq=(st['G230A']&st['R293Q'])-st['R71H']; C_haq=st['R71H']&st['G230A']&st['R293Q']
    d['R85H']=d.research_id.isin(C_r85h).astype(int); d['IRAK3_LoF']=d.research_id.isin(C_irak3).astype(int)
    d['AQ']=d.research_id.isin(C_aq).astype(int); d['HAQ']=d.research_id.isin(C_haq).astype(int)
    d['sting']=np.where(d.HAQ==1,'HAQ',np.where(d.AQ==1,'AQ','WT'))
    print(f"RNA {len(d)} | R85H {int(d.R85H.sum())} | AQ {int(d.AQ.sum())} | HAQ {int(d.HAQ.sum())} | WT {int((d.sting=='WT').sum())}")
    print("   STING class × R85H counts:")
    for cls in ['WT','AQ','HAQ']:
        for r in [0,1]: print(f"      {cls:3s} R85H={r}: n={int(((d.sting==cls)&(d.R85H==r)).sum())}")
    # covariates
    anc=pd.read_csv(ANC,sep="\t"); anc['research_id']=anc.research_id.astype(str); COV=[]
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
        d=d.merge(cf,on='research_id',how='left'); COV+=cc
    except Exception: print("   (cell_fractions_v9.csv not found)")
    covs=(" + "+" + ".join(COV)) if COV else ""
    import statsmodels.formula.api as smf
    def fit(fm,terms):
        try:
            r=smf.ols(fm,data=d,missing='drop').fit(); out={}
            for t in terms: out[t]={'beta':round(float(r.params[t]),3),'p':round(float(r.pvalues[t]),4)} if t in r.params else None
            return out
        except Exception as e: return {'err':str(e)[:60]}
    # ---- (1) MAIN EFFECTS vs clean WT reference: module ~ AQ + HAQ + R85H + IRAK3 ----
    print("\n== (1) MAIN EFFECTS (AQ, HAQ each vs WT reference; clean): module ~ AQ + HAQ + R85H + IRAK3_LoF ==")
    S['main']={}
    for mod in MODULES:
        o=fit(f'{mod} ~ AQ + HAQ + R85H + IRAK3_LoF{covs}',['AQ','HAQ','R85H','IRAK3_LoF']); S['main'][mod]=o
        fmt=lambda t: f"{o[t]['beta']}(p{o[t]['p']})" if o.get(t) else "na"
        print(f"   {mod:12s} AQ {fmt('AQ'):>15} | HAQ {fmt('HAQ'):>15} | R85H {fmt('R85H'):>15} | IRAK3 {fmt('IRAK3_LoF'):>15}")
    # ---- (2) FORMAL INTERACTION: module ~ R85H*AQ + R85H*HAQ (WT ref) -- the epistasis + the HAQ negative control ----
    print("\n== (2) INTERACTION: module ~ R85H*AQ + R85H*HAQ (predict R85H:AQ>0, R85H:HAQ~0/neg) ==")
    S['interaction']={}
    for mod in ['typeI_ISG','typeII_IFNg','NFkB']:
        o=fit(f'{mod} ~ R85H*AQ + R85H*HAQ{covs}',['R85H','AQ','HAQ','R85H:AQ','R85H:HAQ']); S['interaction'][mod]=o
        fmt=lambda t: f"{o[t]['beta']}(p{o[t]['p']})" if o.get(t) else "na"
        print(f"   {mod:12s} R85H {fmt('R85H'):>14} AQ {fmt('AQ'):>14} HAQ {fmt('HAQ'):>14} | R85H:AQ {fmt('R85H:AQ'):>14} R85H:HAQ {fmt('R85H:HAQ'):>14}")
    # ---- (3) 6-way stratified type-I ISG means ----
    print("\n== (3) type-I ISG mean (z) by STING class × R85H — the clean picture ==")
    S['strata']={}
    for cls in ['WT','AQ','HAQ']:
        for r in [0,1]:
            g=d[(d.sting==cls)&(d.R85H==r)]; mn=round(float(g.typeI_ISG.mean()),3) if len(g) else None
            S['strata'][f'{cls}_R85H{r}']={'n':int(len(g)),'mean_typeI_ISG':mn}
            print(f"   {cls:3s} R85H={r}: n={len(g):>5} mean typeI_ISG z = {mn}")
    print("\n== READ: HAQ main <=0? R85H:AQ>0 (unmasked)? R85H:HAQ~0 (HAQ negative control)? ==")
    print("\n===== STING-CLASS IFN (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
