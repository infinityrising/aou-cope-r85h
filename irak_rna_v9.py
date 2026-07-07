#!/usr/bin/env python
"""AoU v9 — IRAK-axis mutations in the RNA-seq cohort + IRAK-mut -> ISG (immune-brake -> interferon).
Answers: are there IRAK1/3/4 (+MYD88) damaging variants among RNA-seq patients, how many, and do they move the ISG score?
Directional biology to check: IRAK3 (=IRAK-M) is a NEGATIVE regulator -> LoF should RAISE ISG; IRAK4 LoF impairs
MyD88 signaling -> may LOWER ISG. Also R85H × IRAK two-hit -> ISG. Damaging = VAT HC-pLoF OR missense REVEL>0.5, rare
(gvs_afr_af<0.1%); LoF-only tallied separately. Carriers via BigQuery ∩ RNA-seq cohort. Within-AFR + ancestry-adjusted.
STANDARD app. Ends 'run complete' / 'run failed'.
"""
import os, subprocess, io, json, gzip
from collections import Counter
import numpy as np, pandas as pd
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9")
RNADIR=f"{MNT}/v9/multiomics/rnaseq"; SNP=f"{MNT}/v9/wgs/short_read/snpindel"
TPM=f"{RNADIR}/rsem/aou_rnaseq_20260415.rsem_genes_tpm.txt.gz"; ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"
VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1); R85H_VID='19-18911007-C-T'
IRAK={'IRAK1':('X',154005501,154025650),'IRAK3':('12',66183995,66259622),'IRAK4':('12',43753938,43803307),'MYD88':('3',38133552,38148024)}
IFN={'IFI27':'ENSG00000165949','IFI44L':'ENSG00000137959','IFIT1':'ENSG00000185745','ISG15':'ENSG00000187608','RSAD2':'ENSG00000134321','SIGLEC1':'ENSG00000088827'}
S={}
def sh(c): return subprocess.run(['bash','-lc',c],capture_output=True,text=True).stdout
def fnum(x):
    try: return float(x)
    except (TypeError,ValueError): return 0.0
try:
    print("== ISG (6-gene) ==")
    raw=sh(f'zcat "{TPM}" 2>/dev/null | awk -F"\\t" \'NR==1 || $1 ~ /{"|".join(IFN.values())}/\'')
    m=pd.read_csv(io.StringIO(raw),sep="\t"); m['ensg']=m.gene_id.str.split('.').str[0]
    m=m[m.ensg.isin(IFN.values())].set_index('ensg').drop(columns=['gene_id','transcript_id(s)'])
    ex=m.T.astype(float); ex.index=ex.index.astype(str); lg=np.log2(ex+1)
    isg=((lg-lg.mean())/lg.std()).mean(axis=1); isgdf=pd.DataFrame({'research_id':isg.index,'ifn':isg.values})
    rna_ids=set(isgdf.research_id)
    print("== IRAK damaging variants (VAT: rare HC-pLoF or missense REVEL>0.5) ==")
    import pysam
    with gzip.open(VAT,'rt') as fh: COLS=fh.readline().rstrip("\n").split("\t")
    IX={c:COLS.index(c) for c in ['vid','gene_symbol','consequence','gvs_afr_af','revel','LoF']}
    tbx=pysam.TabixFile(VAT); per={g:{'lof':[],'mis':[]} for g in IRAK}
    for g,(ch,s,e) in IRAK.items():
        for line in tbx.fetch(f'chr{ch}',s,e):
            f=line.split("\t")
            if f[IX['gene_symbol']]!=g: continue
            if fnum(f[IX['gvs_afr_af']])>=0.001: continue
            if f[IX['LoF']]=='HC': per[g]['lof'].append(f[IX['vid']])
            elif 'missense' in f[IX['consequence']] and fnum(f[IX['revel']])>0.5: per[g]['mis'].append(f[IX['vid']])
    S['damaging_vids_per_gene']={g:{'LoF':len(v['lof']),'mis_REVEL':len(v['mis'])} for g,v in per.items()}
    print("   damaging vids:",S['damaging_vids_per_gene'])
    print("== carriers via BigQuery ∩ RNA-seq ==")
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    def carriers(vids):
        out=set()
        for i in range(0,len(vids),900):
            if not vids[i:i+900]: continue
            inl=",".join("'"+v+"'" for v in vids[i:i+900])
            out|=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid IN ({inl})"))
        return out
    anc=pd.read_csv(ANC,sep="\t",usecols=['research_id','ancestry_pred']); anc['research_id']=anc.research_id.astype(str)
    amap=dict(zip(anc.research_id,anc.ancestry_pred))
    r85h=carriers([R85H_VID]); r85h_rna=r85h & rna_ids
    gene_car={}; gene_car_rna={}
    for g in IRAK:
        c=carriers(per[g]['lof']+per[g]['mis']); gene_car[g]=c; gene_car_rna[g]=c & rna_ids
        by_anc=Counter(amap.get(x,'NA') for x in gene_car_rna[g])
        print(f"   {g}: cohort carriers {len(c)} | ∩RNA {len(gene_car_rna[g])} | RNA by anc {dict(by_anc)}")
    S['carriers_rna_per_gene']={g:len(gene_car_rna[g]) for g in IRAK}
    irak_any_rna=set().union(*(gene_car_rna[g] for g in ['IRAK1','IRAK3','IRAK4']))
    S['IRAK_any_RNA']=len(irak_any_rna); S['R85H_RNA']=len(r85h_rna)
    S['R85H_x_IRAK_RNA_doublecarriers']=len(r85h_rna & irak_any_rna)
    print(f"   ANY IRAK1/3/4 ∩RNA = {len(irak_any_rna)} | R85H∩RNA = {len(r85h_rna)} | R85H×IRAK double ∩RNA = {len(r85h_rna & irak_any_rna)}")
    print("== IRAK-mut -> ISG (IRAK3-LoF negative regulator => expect POSITIVE beta) ==")
    d=isgdf.merge(anc,on='research_id',how='left'); d['R85H']=d.research_id.isin(r85h).astype(int)
    d['IRAK_any']=d.research_id.isin(irak_any_rna).astype(int)
    for g in IRAK: d[g]=d.research_id.isin(gene_car_rna[g]).astype(int)
    a=d[d.ancestry_pred=='afr']
    import statsmodels.formula.api as smf
    def ols(f,dd,term):
        try: r=smf.ols(f,data=dd.dropna(subset=['ifn'])).fit(); return {'beta':round(float(r.params[term]),4),'p':round(float(r.pvalues[term]),4),'n':int(r.nobs),'n_carr':int(dd[term].sum()) if term in dd else None}
        except Exception as e: return {'error':str(e)[:60]}
    for term in ['IRAK_any','IRAK3','IRAK4','IRAK1','MYD88']:
        S[f'isg_{term}_all']=ols(f'ifn ~ {term} + C(ancestry_pred)',d,term); S[f'isg_{term}_afr']=ols(f'ifn ~ {term}',a,term)
        print(f"   ISG~{term}: ALL {S[f'isg_{term}_all']} | AFR {S[f'isg_{term}_afr']}")
    print("== R85H × IRAK two-hit -> ISG ==")
    d['inter']=d.R85H*d.IRAK_any
    S['R85H_x_IRAK']=ols('ifn ~ R85H + IRAK_any + inter + C(ancestry_pred)',d,'inter')
    print("   R85H×IRAK_any:",S['R85H_x_IRAK'])
    print("\n===== IRAK-in-RNA (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1000:])
