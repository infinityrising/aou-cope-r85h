#!/usr/bin/env python
"""AoU v9 -- PHASED STING-class IFN (the definitive version). Fixes the cis-proof problem: unphased 'carries G230A+R293Q'
does NOT prove they are on the SAME allele -- G230A|R293Q in TRANS is not a real AQ allele. AQ/HAQ/non-HAQ conclusions
require PHASED sequence. This uses long-read CIS-phased AQ_d / HAQ_d (from sting_phenome_pav_v9.csv), restricted to RNA∩PAV.
The overlap is high (~8.3k of ~9k RNA samples are phased), so this is BOTH rigorous AND well-powered for the molecular
readout (unlike clinical endpoints, where PAV=13k is underpowered and unphased whole-cohort is the LD-justified proxy).
Model: module ~ R85H*cisAQ + R85H*cisHAQ (WT/neither = cis reference). Predictions: cisHAQ main <=0 (dampened);
R85H:cisAQ > 0 (unmasked on the gain allele); R85H:cisHAQ ~ 0/neg (HAQ = the PHASED negative control). Reports the 6-way
stratified type-I ISG means and the phased-vs-unphased class agreement. Standard app. Ends 'run complete'/'run failed'.
"""
import os, io, gzip, subprocess, json
import numpy as np, pandas as pd
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); RNADIR=f"{MNT}/v9/multiomics/rnaseq"; SNP=f"{MNT}/v9/wgs/short_read/snpindel"
TPM=f"{RNADIR}/rsem/aou_rnaseq_20260415.rsem_genes_tpm.txt.gz"; ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CELL=os.path.expanduser("~/cell_fractions_v9.csv"); PHPAV=os.path.expanduser("~/sting_phenome_pav_v9.csv")
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
R85H_VID='19-18911007-C-T'; STINGV={'R71H':'5-139481493-C-T','G230A':'5-139478340-C-G','R293Q':'5-139477397-C-T'}
MODULES={
 'typeI_specific':{'SIGLEC1':'ENSG00000088827','IFI27':'ENSG00000165949','USP18':'ENSG00000184979','IFI6':'ENSG00000126709','IFI44L':'ENSG00000137959'},  # type-I-PREFERENTIAL (Nombel 2025; Kim 2024)
 'typeII_IFNg':{'CXCL9':'ENSG00000138755','GBP5':'ENSG00000154451','GBP1':'ENSG00000117228','IDO1':'ENSG00000131203','ANKRD22':'ENSG00000152766'},  # IFN-γ-SPECIFIC (CXCL10/11 dropped -> shared)
 'sharedISG_tone':{'ISG15':'ENSG00000187608','RSAD2':'ENSG00000134321','IFIT1':'ENSG00000185745','MX1':'ENSG00000157601','OAS1':'ENSG00000089127'},  # pan type-I/II = IFN tone, not a type
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
    d['research_id']=d.index.astype(str); n_rna=len(d)
    # ---- PHASED STING class from PAV (cis) ----
    pav=pd.read_csv(PHPAV); pav['research_id']=pav.research_id.astype(str)
    d=d.merge(pav[['research_id','AQ_d','HAQ_d']],on='research_id',how='inner')   # RNA ∩ PAV (phased)
    d['cisAQ']=(d.AQ_d>=1).astype(int); d['cisHAQ']=(d.HAQ_d>=1).astype(int)
    d.loc[d.cisHAQ==1,'cisAQ']=0                                                  # HAQ takes precedence (mutually exclusive)
    d['sting']=np.where(d.cisHAQ==1,'HAQ',np.where(d.cisAQ==1,'AQ','WT'))
    print(f"RNA {n_rna} | RNA∩PAV(phased) {len(d)} ({100*len(d)/n_rna:.0f}% phased) | cisAQ {int(d.cisAQ.sum())} | cisHAQ {int(d.cisHAQ.sum())} | WT {int((d.sting=='WT').sum())}")
    # ---- exposures: R85H + IRAK3 (from BQ/VAT) ----
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
    C_r85h=carr([R85H_VID]); C_irak3=carr(irak3)
    d['R85H']=d.research_id.isin(C_r85h).astype(int); d['IRAK3_LoF']=d.research_id.isin(C_irak3).astype(int)
    print("   STING(cis) × R85H counts:")
    for cls in ['WT','AQ','HAQ']:
        for r in [0,1]: print(f"      {cls:3s} R85H={r}: n={int(((d.sting==cls)&(d.R85H==r)).sum())}")
    # ---- covariates ----
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
    try:  # smoking (EHR ever-smoker, ICD tobacco): whole-blood ISG/inflammatory expression confounder -- parity with EMR + RNA models
        smk=set(str(r[0]) for r in bq.query(f"SELECT DISTINCT co.person_id FROM `{PROJ}.{DS}.condition_occurrence` co JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id WHERE c.vocabulary_id LIKE 'ICD%' AND (c.concept_code LIKE 'Z72.0%' OR c.concept_code LIKE 'F17%' OR c.concept_code LIKE '305.1%' OR c.concept_code LIKE 'V15.82%')"))
        d['smoker']=d.research_id.isin(smk).astype(float); COV+=['smoker']; print(f"   adj: ever-smoker {int(d.smoker.sum())} added (smoking parity)")
    except Exception: pass
    covs=(" + "+" + ".join(COV)) if COV else ""
    import statsmodels.formula.api as smf
    def fit(fm,terms):
        try:
            r=smf.ols(fm,data=d,missing='drop').fit()
            return {t:({'beta':round(float(r.params[t]),3),'p':round(float(r.pvalues[t]),4)} if t in r.params else None) for t in terms}
        except Exception as e: return {'err':str(e)[:60]}
    print("\n== (1) PHASED MAIN EFFECTS (cisAQ, cisHAQ vs cis-WT reference): module ~ cisAQ + cisHAQ + R85H + IRAK3_LoF ==")
    S['main']={}
    for mod in MODULES:
        o=fit(f'{mod} ~ cisAQ + cisHAQ + R85H + IRAK3_LoF{covs}',['cisAQ','cisHAQ','R85H','IRAK3_LoF']); S['main'][mod]=o
        fmt=lambda t: f"{o[t]['beta']}(p{o[t]['p']})" if o.get(t) else "na"
        print(f"   {mod:12s} cisAQ {fmt('cisAQ'):>15} | cisHAQ {fmt('cisHAQ'):>15} | R85H {fmt('R85H'):>15} | IRAK3 {fmt('IRAK3_LoF'):>15}")
    print("\n== (2) PHASED INTERACTION: module ~ R85H*cisAQ + R85H*cisHAQ  (predict R85H:cisAQ>0, R85H:cisHAQ~0/neg) ==")
    S['interaction']={}
    for mod in ['typeI_specific','typeII_IFNg','sharedISG_tone','NFkB']:
        o=fit(f'{mod} ~ R85H*cisAQ + R85H*cisHAQ{covs}',['R85H','cisAQ','cisHAQ','R85H:cisAQ','R85H:cisHAQ']); S['interaction'][mod]=o
        fmt=lambda t: f"{o[t]['beta']}(p{o[t]['p']})" if o.get(t) else "na"
        print(f"   {mod:12s} R85H {fmt('R85H'):>13} cisAQ {fmt('cisAQ'):>13} cisHAQ {fmt('cisHAQ'):>13} | R85H:cisAQ {fmt('R85H:cisAQ'):>13} R85H:cisHAQ {fmt('R85H:cisHAQ'):>13}")
    print("\n== (3) TYPE-I-SPECIFIC score mean (z) by PHASED STING class × R85H ==")
    S['strata']={}
    for cls in ['WT','AQ','HAQ']:
        for r in [0,1]:
            g=d[(d.sting==cls)&(d.R85H==r)]; mn=round(float(g.typeI_specific.mean()),3) if len(g) else None
            S['strata'][f'{cls}_R85H{r}']={'n':int(len(g)),'mean_typeI_specific':mn}
            print(f"   {cls:3s} R85H={r}: n={len(g):>5} mean typeI_specific z = {mn}")
    print("\n== READ (PHASED = definitive): is R85H:cisAQ>0 on typeI_SPECIFIC (SIGLEC1/IFI27/USP18) but ~0 on typeII_IFNg -> TYPE-I-SELECTIVE unmasking? cisHAQ main<=0? R85H:cisHAQ~0 (phased neg control)? ==")
    print("\n===== PHASED STING-CLASS IFN (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
