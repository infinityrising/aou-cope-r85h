#!/usr/bin/env python
"""AoU v9 -- CELL-INTRINSIC STING arms vs IFN-recruitment (ZS's central mechanistic question). STING has a canonical arm
(TBK1->IRF3->type-I IFN-a/b -> ISG + chemokine immune RECRUITMENT) and cell-intrinsic arms (NF-kB, AUTOPHAGY, APOPTOSIS,
SENESCENCE/SASP). Our data: AQ (STING gain) -> NF-kB (not type-I) + fibrotic-ILD. Hypothesis: the fibrotic phenotype is
driven by the cell-intrinsic SENESCENCE/SASP + NF-kB arm, not IFN-a/recruitment. Test in whole-blood RNA (immune-cell
senescence/autophagy/apoptosis programs + circulating SASP) whether AQ / R85H / IRAK3-LoF / cis-R85H×AQ partition toward
cell-intrinsic (senescence/SASP/autophagy/apoptosis) vs recruitment (type-I ISG / IFN-g chemokines). PHASED (RNA∩PAV, cis).
Gene panels are a first pass (lit agent will refine); per-module coverage reported. Standard app. Ends 'run complete'/'run failed'.
"""
import os, io, gzip, subprocess, json
import numpy as np, pandas as pd
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); RNADIR=f"{MNT}/v9/multiomics/rnaseq"; SNP=f"{MNT}/v9/wgs/short_read/snpindel"
TPM=f"{RNADIR}/rsem/aou_rnaseq_20260415.rsem_genes_tpm.txt.gz"; ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CELL=os.path.expanduser("~/cell_fractions_v9.csv"); PHPAV=os.path.expanduser("~/sting_phenome_pav_v9.csv")
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
R85H_VID='19-18911007-C-T'; STINGV={'R71H':'5-139481493-C-T','G230A':'5-139478340-C-G','R293Q':'5-139477397-C-T'}
MODULES={  # panels per verified lit synthesis (SenMayo/Saul 2022; HALLMARK; Klionsky 2021). CXCL9/10/11 kept OUT of NF-kB/SASP (IFN-inducible).
 'typeI_ISG':{'IFI27':'ENSG00000165949','IFI44L':'ENSG00000137959','IFIT1':'ENSG00000185745','ISG15':'ENSG00000187608','RSAD2':'ENSG00000134321','SIGLEC1':'ENSG00000088827','MX1':'ENSG00000157601','OAS1':'ENSG00000089127'},
 'IFNg_recruit':{'CXCL9':'ENSG00000138755','CXCL10':'ENSG00000169245','CXCL11':'ENSG00000169248','GBP1':'ENSG00000117228','GBP5':'ENSG00000154451','STAT1':'ENSG00000115415','IDO1':'ENSG00000131203'},
 'NFkB_axis':{'NFKB1':'ENSG00000109320','RELA':'ENSG00000173039','NFKBIA':'ENSG00000100906','TNFAIP3':'ENSG00000118503','BIRC3':'ENSG00000023445','PTGS2':'ENSG00000073756','ICAM1':'ENSG00000090339','SOD2':'ENSG00000112096','TNF':'ENSG00000232810'},
 'senescence_SASP':{'CDKN1A':'ENSG00000124762','CDKN2A':'ENSG00000147889','SERPINE1':'ENSG00000106366','GDF15':'ENSG00000130513','IGFBP3':'ENSG00000146674','IGFBP7':'ENSG00000163453','TIMP1':'ENSG00000102265','IL6':'ENSG00000136244','IL1B':'ENSG00000125538','CXCL8':'ENSG00000169429','CCL2':'ENSG00000108691','MMP9':'ENSG00000100985','INHBA':'ENSG00000122641','GLB1':'ENSG00000170266'},
 'autophagy':{'MAP1LC3B':'ENSG00000140941','SQSTM1':'ENSG00000161011','GABARAPL1':'ENSG00000139112','ATG5':'ENSG00000057663','ATG7':'ENSG00000197548','ATG16L1':'ENSG00000085978','BECN1':'ENSG00000126581','ULK1':'ENSG00000177169','WIPI2':'ENSG00000157954','NBR1':'ENSG00000188554'},
 'apoptosis':{'BAX':'ENSG00000087088','BAK1':'ENSG00000030110','BBC3':'ENSG00000105327','PMAIP1':'ENSG00000141682','BCL2L11':'ENSG00000153094','BID':'ENSG00000015475','CASP3':'ENSG00000164305','CASP8':'ENSG00000064012','CASP9':'ENSG00000132906','GADD45A':'ENSG00000116717'},
 'fibrosis_TGFb':{'TGFB1':'ENSG00000105329','CCN2':'ENSG00000118523','FN1':'ENSG00000115414','COL1A1':'ENSG00000108821','ACTA2':'ENSG00000107796','TIMP1':'ENSG00000102265'},
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
    found=set(m.index); d=pd.DataFrame(index=z.index); cov={}
    for mod,genes in MODULES.items():
        got=[g for g in genes.values() if g in found]; cov[mod]=f"{len(got)}/{len(genes)}"; d[mod]=z[got].mean(axis=1) if got else np.nan
    print("module gene coverage:",cov); S['coverage']=cov
    d['research_id']=d.index.astype(str)
    pav=pd.read_csv(PHPAV); pav['research_id']=pav.research_id.astype(str)
    d=d.merge(pav[['research_id','AQ_d','HAQ_d']],on='research_id',how='inner')
    d['cisAQ']=(d.AQ_d>=1).astype(int); d['cisHAQ']=(d.HAQ_d>=1).astype(int); d.loc[d.cisHAQ==1,'cisAQ']=0
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
    print(f"RNA∩PAV {len(d)} | cisAQ {int(d.cisAQ.sum())} | cisHAQ {int(d.cisHAQ.sum())} | R85H {int(d.R85H.sum())} | IRAK3-LoF {int(d.IRAK3_LoF.sum())}")
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
    except Exception: print("   (cell_fractions not found)")
    try:  # smoking (EHR ever-smoker, ICD tobacco): whole-blood ISG/inflammatory expression confounder
        smk=set(str(r[0]) for r in bq.query(f"SELECT DISTINCT co.person_id FROM `{PROJ}.{DS}.condition_occurrence` co JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id WHERE c.vocabulary_id LIKE 'ICD%' AND (c.concept_code LIKE 'Z72.0%' OR c.concept_code LIKE 'F17%' OR c.concept_code LIKE '305.1%' OR c.concept_code LIKE 'V15.82%')"))
        d['smoker']=d.research_id.isin(smk).astype(float); COV+=['smoker']; print(f"   adj: ever-smoker {int(d.smoker.sum())} added to RNA covariates")
    except Exception: pass
    covs=(" + "+" + ".join(COV)) if COV else ""
    import statsmodels.formula.api as smf
    def fit(fm,terms):
        try:
            r=smf.ols(fm,data=d,missing='drop').fit()
            return {t:({'beta':round(float(r.params[t]),3),'p':round(float(r.pvalues[t]),4)} if t in r.params else None) for t in terms}
        except Exception as e: return {'err':str(e)[:50]}
    print("\n== MAIN EFFECTS by module (cisAQ, cisHAQ, R85H, IRAK3_LoF vs cis-WT ref) — CELL-INTRINSIC vs RECRUITMENT ==")
    print(f"{'module':16s} {'cisAQ':>16} {'cisHAQ':>14} {'R85H':>14} {'IRAK3_LoF':>16}")
    S['main']={}
    for mod in MODULES:
        o=fit(f'{mod} ~ cisAQ + cisHAQ + R85H + IRAK3_LoF{covs}',['cisAQ','cisHAQ','R85H','IRAK3_LoF']); S['main'][mod]=o
        g=lambda t: f"{o[t]['beta']}(p{o[t]['p']})" if o.get(t) else "na"
        print(f"{mod:16s} {g('cisAQ'):>16} {g('cisHAQ'):>14} {g('R85H'):>14} {g('IRAK3_LoF'):>16}")
    print("\n== R85H×cisAQ INTERACTION by module (does R85H unmask cell-intrinsic vs recruitment on the AQ allele?) ==")
    S['interaction']={}
    for mod in MODULES:
        o=fit(f'{mod} ~ R85H*cisAQ + R85H*cisHAQ{covs}',['R85H:cisAQ','R85H:cisHAQ']); S['interaction'][mod]=o
        g=lambda t: f"{o[t]['beta']}(p{o[t]['p']})" if o.get(t) else "na"
        print(f"   {mod:16s} R85H:cisAQ {g('R85H:cisAQ'):>15}  R85H:cisHAQ {g('R85H:cisHAQ'):>15}")
    print("\n== READ: cell-intrinsic (senescence_SASP/autophagy/apoptosis/fibrosis) vs recruitment (typeI_ISG/IFNg). Does AQ favor SENESCENCE+NF-kB over type-I? ==")
    print("\n===== CELL-INTRINSIC RNA (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
