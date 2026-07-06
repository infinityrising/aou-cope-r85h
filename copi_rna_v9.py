#!/usr/bin/env python
"""AoU v9 COPE-R85H — COPI coatomer subunit RNA expression (molecular pillar: COPI + STING -> IFN). STANDARD app.
Extract the 9 COPI subunit TPMs from RSEM; per-subunit z + a COPI-complex score; relate to R85H genotype and the
IFN score (does R85H alter COPE mRNA / COPI expression? does lower COPI expression track higher IFN?). Within-AFR.
Ends 'run complete' / 'run failed'.
"""
import os, subprocess, io, json
import numpy as np, pandas as pd
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
RNADIR=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9/v9/multiomics/rnaseq")
TPM=f"{RNADIR}/rsem/aou_rnaseq_20260415.rsem_genes_tpm.txt.gz"
ANC=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv")
COPI={'COPA':'ENSG00000122218','COPB1':'ENSG00000129083','COPB2':'ENSG00000184432','COPG1':'ENSG00000181789',
      'COPG2':'ENSG00000158623','COPZ1':'ENSG00000111481','COPZ2':'ENSG00000005243','COPE':'ENSG00000105669','ARCN1':'ENSG00000095139'}
IFNG={'IFI27':'ENSG00000165949','IFI44L':'ENSG00000137959','IFIT1':'ENSG00000185745','ISG15':'ENSG00000187608','RSAD2':'ENSG00000134321','SIGLEC1':'ENSG00000088827'}
R85H_VID='19-18911007-C-T'; S={}
def sh(c): return subprocess.run(['bash','-lc',c],capture_output=True,text=True).stdout
def extract(genes):
    raw=sh(f'zcat "{TPM}" 2>/dev/null | awk -F"\\t" \'NR==1 || $1 ~ /{"|".join(genes.values())}/\'')
    m=pd.read_csv(io.StringIO(raw),sep="\t"); m['ensg']=m.gene_id.str.split('.').str[0]
    m=m[m.ensg.isin(genes.values())].set_index('ensg').drop(columns=['gene_id','transcript_id(s)'])
    inv={v:k for k,v in genes.items()}; e=m.T.astype(float); e.index=e.index.astype('int64'); e.columns=[inv[c] for c in e.columns]
    return e
try:
    print("== A. COPI subunit TPM ==")
    copi=extract(COPI); print("   subunits found:", list(copi.columns))
    zc=(np.log2(copi+1)-np.log2(copi+1).mean())/np.log2(copi+1).std()
    ifn=extract(IFNG); zi=np.log2(ifn+1); ifn_score=((zi-zi.mean())/zi.std()).mean(axis=1)
    d=pd.DataFrame({'copi_score':zc.mean(axis=1),'ifn_score':ifn_score}); d['research_id']=d.index
    for k in copi.columns: d[k]=zc[k]
    print("== B. R85H + ancestry ==")
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    r85h=set(int(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid='{R85H_VID}'"))
    anc=pd.read_csv(ANC,sep="\t",usecols=['research_id','ancestry_pred']); anc['research_id']=anc.research_id.astype('int64')
    d=d.merge(anc,on='research_id',how='left'); d['R85H']=d.research_id.isin(r85h).astype(int)
    a=d[d.ancestry_pred=='afr']
    print("== C. tests (within-AFR) ==")
    import statsmodels.formula.api as smf
    def ols(f,dd,term):
        dd=dd.dropna(subset=[term.split()[0]] if ' ' not in term else None) if False else dd
        try: r=smf.ols(f,data=dd).fit(); return {'beta':round(float(r.params[term]),4),'p':round(float(r.pvalues[term]),4),'n':int(r.nobs)}
        except Exception as ex: return {'error':str(ex)[:60]}
    S['R85H_to_COPE']=ols('COPE ~ R85H',a,'R85H')
    S['R85H_to_COPIscore']=ols('copi_score ~ R85H',a,'R85H')
    S['R85H_per_subunit']={k:ols(f'Q("{k}") ~ R85H',a,'R85H') for k in copi.columns}
    S['COPIscore_to_IFN']=ols('ifn_score ~ copi_score',a,'copi_score')
    print("   R85H->COPE:",S['R85H_to_COPE']," | R85H->COPIscore:",S['R85H_to_COPIscore']," | COPIscore->IFN:",S['COPIscore_to_IFN'])
    print("   R85H per-subunit:",S['R85H_per_subunit'])
    print("\n===== COPI SUMMARY (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1000:])
