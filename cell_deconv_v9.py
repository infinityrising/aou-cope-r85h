#!/usr/bin/env python
"""AoU v9 — cell-composition estimation for whole-blood RNA-seq (covariate for the ISG regression; ZS: full deconvolution).
Bulk whole-blood ISG is confounded by leukocyte composition (monocytes + pDCs make type-I IFN; neutrophils dilute).
Computes per-person reference-MARKER cell-type scores = mean z of log2(TPM+1) over canonical lineage markers, for the
major blood lineages -> ~/cell_fractions_v9.csv, to adjust ISG ~ R85H*HAQ. (True LM22/CIBERSORT fractions need an
external signature not downloadable on the locked AoU VM; marker scores span the same composition axes.) Self-QC prints
markers-found per type (catches a wrong ENSG -> that marker silently drops, score uses the rest). NOTE SIGLEC1 is BOTH
an ISG-panel gene AND a monocyte marker -> adjusting for monocyte score partly removes it; the analyzer also reports a
5-gene ISG (SIGLEC1 dropped) as sensitivity. STANDARD app. Ends 'run complete' / 'run failed'.
"""
import os, subprocess, io
import numpy as np, pandas as pd
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); RNADIR=f"{MNT}/v9/multiomics/rnaseq"
TPM=f"{RNADIR}/rsem/aou_rnaseq_20260415.rsem_genes_tpm.txt.gz"
OUT=os.path.expanduser("~/cell_fractions_v9.csv")
MARKERS={
 'neut':{'FCGR3B':'ENSG00000162747','CSF3R':'ENSG00000119535','S100A8':'ENSG00000143546','S100A9':'ENSG00000163220','CXCR2':'ENSG00000180871','FUT4':'ENSG00000196371'},
 'mono':{'CD14':'ENSG00000170458','LYZ':'ENSG00000090382','FCN1':'ENSG00000085265','VCAN':'ENSG00000038427','S100A12':'ENSG00000163221','CSF1R':'ENSG00000182578'},
 'cd4t':{'CD3D':'ENSG00000167286','CD3E':'ENSG00000198851','IL7R':'ENSG00000168685','CD4':'ENSG00000010610'},
 'cd8t':{'CD8A':'ENSG00000153563','CD8B':'ENSG00000172116','GZMK':'ENSG00000113088'},
 'bcell':{'MS4A1':'ENSG00000156738','CD79A':'ENSG00000105369','CD79B':'ENSG00000007312','CD19':'ENSG00000177455','BANK1':'ENSG00000153064'},
 'nk':{'NKG7':'ENSG00000105374','KLRD1':'ENSG00000134539','GNLY':'ENSG00000115523','NCAM1':'ENSG00000149294'},
 'pdc':{'LILRA4':'ENSG00000239961','CLEC4C':'ENSG00000198178','IL3RA':'ENSG00000185291','GZMB':'ENSG00000100453'},
}
def sh(c): return subprocess.run(['bash','-lc',c],capture_output=True,text=True).stdout
try:
    allids={e for d in MARKERS.values() for e in d.values()}; pat="|".join(sorted(allids))
    raw=sh(f'zcat "{TPM}" 2>/dev/null | awk -F"\\t" \'NR==1 || $1 ~ /{pat}/\'')
    m=pd.read_csv(io.StringIO(raw),sep="\t"); m['ensg']=m.gene_id.str.split('.').str[0]
    m=m[m.ensg.isin(allids)].set_index('ensg').drop(columns=['gene_id','transcript_id(s)'])
    expr=m.T.astype(float); expr.index=expr.index.astype(str); lg=np.log2(expr+1); z=(lg-lg.mean())/lg.std()
    raw=pd.DataFrame(index=expr.index); qc={}
    for ct,dd in MARKERS.items():
        present=[e for e in dd.values() if e in z.columns]; qc[ct]=f"{len(present)}/{len(dd)}"
        raw[ct]=z[present].mean(axis=1) if present else np.nan
    print("markers found per type:",qc)
    print("RAW (cross-sample z) corr — a global technical axis inflates ALL pairs (myeloid-lymphoid NOT negative):")
    print(raw.corr().round(2).to_string())
    # --- compositional fix (CLR-like): remove the per-sample global axis so scores become RELATIVE composition ---
    libaxis=raw.mean(axis=1)                 # global technical/RNA-quality/detection axis (kept as a covariate)
    cent=raw.sub(libaxis,axis=0)             # per-sample centered => relative composition (myeloid vs lymphoid now trade off)
    print("CENTERED (compositional) corr — myeloid (neut/mono) vs lymphoid (cd4t/cd8t/bcell/nk) should now be NEGATIVE:")
    print(cent.corr().round(2).to_string())
    REF='neut'                               # drop one (reference) to avoid sum-to-0 collinearity; keep mono+pdc (ISG-relevant)
    out=cent.drop(columns=[REF]).copy(); out['libaxis']=libaxis; out['research_id']=out.index
    out.to_csv(OUT,index=False)
    print(f"cell composition (6 centered + libaxis; ref={REF}) -> {OUT} ({len(out)} persons)")
    print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1000:])
