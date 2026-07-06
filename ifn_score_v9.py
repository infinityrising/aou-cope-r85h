#!/usr/bin/env python
"""AoU v9 COPE-R85H — IFN-score arm (the plan's SENSITIVE molecular primary). STANDARD app (pandas + BigQuery).
6-gene type-I interferon score (Rice/Crow) on whole-blood RNA-seq (RSEM gene TPM) -> test R85H -> IFN score, within-AFR.
Extracts only the 6 IFN genes (memory-light awk, never loads the whole matrix). Ends 'run complete' / 'run failed'.
CAVEATS (first pass): raw z-score IFN (no inverse-normal / cell-composition / batch adjust — refine if signal). n is modest (39 AFR R85H w/ RNA).
"""
import os, subprocess, io, json
import numpy as np, pandas as pd
CDR = "wb-silky-artichoke-2408.C2025Q4R6"; PROJ, DS = CDR.split(".", 1)
RNADIR = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9/v9/multiomics/rnaseq")
TPM = f"{RNADIR}/rsem/aou_rnaseq_20260415.rsem_genes_tpm.txt.gz"
MANIFEST = f"{RNADIR}/manifest.tsv"
ANC = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv")
IFN = {'IFI27':'ENSG00000165949','IFI44L':'ENSG00000137959','IFIT1':'ENSG00000185745',
       'ISG15':'ENSG00000187608','RSAD2':'ENSG00000134321','SIGLEC1':'ENSG00000088827'}
R85H_VID = '19-18911007-C-T'
S = {}
def sh(c): return subprocess.run(['bash','-lc',c], capture_output=True, text=True).stdout

try:
    print("== A. extract the 6 IFN genes from RSEM TPM (memory-light) ==")
    pat = "|".join(IFN.values())
    raw = sh(f'zcat "{TPM}" 2>/dev/null | awk -F"\\t" \'NR==1 || $1 ~ /{pat}/\'')
    m = pd.read_csv(io.StringIO(raw), sep="\t")
    m['ensg'] = m['gene_id'].str.split('.').str[0]
    m = m[m.ensg.isin(IFN.values())].set_index('ensg').drop(columns=['gene_id', 'transcript_id(s)'])
    print("   genes found:", sorted({k for k,v in IFN.items() if v in m.index}), "| samples:", m.shape[1])
    expr = m.T.astype(float)                    # samples x 6 genes
    expr.index = expr.index.astype('int64')

    print("== B. IFN score = mean z of log2(TPM+1) ==")
    lg = np.log2(expr + 1.0)
    z = (lg - lg.mean()) / lg.std()
    ifn = z.mean(axis=1).rename('ifn_score')
    print(f"   IFN score for {len(ifn):,} samples | mean {ifn.mean():.3f} sd {ifn.std():.3f}")

    print("== C. join genotype + ancestry ==")
    from google.cloud import bigquery
    bq = bigquery.Client(); T = f"`{PROJ}.{DS}.cb_variant_to_person`"
    r85h = set(int(r.pid) for r in bq.query(
        f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid='{R85H_VID}'"))
    anc = pd.read_csv(ANC, sep="\t", usecols=['research_id','ancestry_pred']); anc['research_id'] = anc.research_id.astype('int64')
    d = pd.DataFrame({'ifn_score': ifn}); d['research_id'] = d.index
    ov = d.research_id.isin(set(anc.research_id)).mean()
    print(f"   sample-id ∩ research_id overlap = {ov:.1%}  (if low, columns are sampleid -> need manifest map)")
    if ov < 0.5:                                # columns are sampleid, map via manifest
        man = pd.read_csv(MANIFEST, sep="\t", usecols=['sampleid','research_id'])
        s2r = dict(zip(man.sampleid.astype('int64'), man.research_id.astype('int64')))
        d['research_id'] = d.research_id.map(s2r)
    d = d.dropna(subset=['research_id']).merge(anc, on='research_id', how='left')
    d['R85H'] = d.research_id.isin(r85h).astype(int)
    print(f"   RNA-seq w/ ancestry: {d.ancestry_pred.notna().sum():,} | R85H w/ RNA: {int(d.R85H.sum())} | AFR R85H w/ RNA: {int(((d.R85H==1)&(d.ancestry_pred=='afr')).sum())}")

    print("== D. test R85H -> IFN score ==")
    import statsmodels.formula.api as smf
    def test(dd, label):
        dd = dd.dropna(subset=['ifn_score','R85H'])
        if dd.R85H.sum() < 5: print(f"   {label}: too few R85H carriers ({int(dd.R85H.sum())})"); return None
        r = smf.ols('ifn_score ~ R85H', data=dd).fit()
        res = {'beta': round(float(r.params['R85H']),4), 'p': float(r.pvalues['R85H']), 'n': int(len(dd)), 'n_R85H': int(dd.R85H.sum()),
               'ifn_R85H_mean': round(float(dd[dd.R85H==1].ifn_score.mean()),3), 'ifn_ref_mean': round(float(dd[dd.R85H==0].ifn_score.mean()),3)}
        print(f"   {label}: beta={res['beta']} p={res['p']:.3g} | IFN mean R85H={res['ifn_R85H_mean']} vs ref={res['ifn_ref_mean']} (n={res['n']}, R85H={res['n_R85H']})")
        return res
    S['IFN_all'] = test(d, 'R85H->IFN (all)')
    S['IFN_afr'] = test(d[d.ancestry_pred == 'afr'], 'R85H->IFN (within-AFR)')

    print("\n===== IFN SUMMARY (paste back) =====")
    print(json.dumps(S, indent=1, default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:", type(e).__name__, str(e)[:300]); print(traceback.format_exc()[-1000:])
