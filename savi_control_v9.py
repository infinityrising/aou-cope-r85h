#!/usr/bin/env python
"""AoU v9 — SAVI POSITIVE CONTROL for the ISG score (PRE-REG self-validating ladder = SAVI↑ / AQ↓).
I had run only the hypomorph (down) end (AQ/HAQ->ISG); this adds the GAIN-OF-FUNCTION anchor the pre-reg requires.
SAVI = constitutively-active STING1 (recurrent GoF hotspot residues 147/154/155/158/166/206/281/284, R232 numbering)
=> known HIGH type-I IFN. Carriers (VAT->BQ) ∩ RNA -> ISG (expect HIGH => validates the score reads STING ACTIVATION,
not just cell state). SAVI is ULTRA-RARE (~17 cohort-wide) so ∩RNA is likely n=0-2 = QUALITATIVE control (even 1-2
GoF carriers with a coordinated high ISG validates the score). Also rare STING1-missense -> ISG as a broader anchor.
Auto-detects the VAT protein-change column + prints examples to confirm numbering. STANDARD app. Ends 'run complete'/'run failed'.
"""
import os, subprocess, io, gzip, json, re
import numpy as np, pandas as pd
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9")
RNADIR=f"{MNT}/v9/multiomics/rnaseq"; SNP=f"{MNT}/v9/wgs/short_read/snpindel"
TPM=f"{RNADIR}/rsem/aou_rnaseq_20260415.rsem_genes_tpm.txt.gz"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
IFN={'IFI27':'ENSG00000165949','IFI44L':'ENSG00000137959','IFIT1':'ENSG00000185745','ISG15':'ENSG00000187608','RSAD2':'ENSG00000134321','SIGLEC1':'ENSG00000088827'}
SAVI_RES={147,154,155,158,166,206,281,284}          # R232-numbering recurrent SAVI GoF hotspots
STING_IV=('chr5',139475528,139482935); S={}
def sh(c): return subprocess.run(['bash','-lc',c],capture_output=True,text=True).stdout
def fnum(x):
    try: return float(x)
    except (TypeError,ValueError): return 0.0
try:
    print("== ISG (6-gene + per-gene z) ==")
    raw=sh(f'zcat "{TPM}" 2>/dev/null | awk -F"\\t" \'NR==1 || $1 ~ /{"|".join(IFN.values())}/\'')
    m=pd.read_csv(io.StringIO(raw),sep="\t"); m['ensg']=m.gene_id.str.split('.').str[0]
    m=m[m.ensg.isin(IFN.values())].set_index('ensg').drop(columns=['gene_id','transcript_id(s)'])
    ex=m.T.astype(float); ex.index=ex.index.astype(str); lg=np.log2(ex+1); zz=(lg-lg.mean())/lg.std()
    isg=zz.mean(axis=1); inv={v:k for k,v in IFN.items()}; zz.columns=[inv[c] for c in zz.columns]
    zz['ifn']=isg; zz['research_id']=zz.index; rna_ids=set(zz.research_id)
    print("== VAT: STING1 missense + SAVI hotspots ==")
    import pysam
    with gzip.open(VAT,'rt') as fh: COLS=fh.readline().rstrip("\n").split("\t")
    pcol=next((c for c in COLS if re.search(r'protein|aa_change|hgvsp|amino',c,re.I)),None)
    print("   protein-change column:",pcol)
    IX={c:COLS.index(c) for c in (['vid','gene_symbol','consequence','gvs_afr_af']+([pcol] if pcol else [])) if c in COLS}
    tbx=pysam.TabixFile(VAT); miss=[]
    for line in tbx.fetch(*STING_IV):
        f=line.split("\t")
        if f[IX['gene_symbol']] not in ('STING1','TMEM173'): continue
        if 'missense' not in f[IX['consequence']]: continue
        aa=f[IX[pcol]] if (pcol and pcol in IX) else ''
        mm=re.search(r'p\.\(?[A-Za-z]{3}(\d+)',aa or ''); res=int(mm.group(1)) if mm else None  # residue #, not the ENSP id
        miss.append((f[IX['vid']],aa,res,fnum(f[IX['gvs_afr_af']])))
    print("   example STING1 missense (vid, aa):", [(v,aa) for (v,aa,r,af) in miss[:5]])
    savi=[v for (v,aa,r,af) in miss if r in SAVI_RES]
    S['n_sting1_missense']=len(miss); S['savi_hits']=[(v,aa) for (v,aa,r,af) in miss if r in SAVI_RES]
    print(f"   STING1 missense vids {len(miss)} | SAVI-hotspot vids {len(savi)}: {S['savi_hits']}")
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    def carriers(vids):
        out=set()
        for i in range(0,len(vids),900):
            if not vids[i:i+900]: continue
            inl=",".join("'"+v+"'" for v in vids[i:i+900])
            out|=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid IN ({inl})"))
        return out
    savi_car=carriers(savi); savi_rna=savi_car & rna_ids
    S['savi_carriers_cohort']=len(savi_car); S['savi_carriers_RNA']=len(savi_rna)
    print(f"   SAVI carriers: cohort {len(savi_car)} | ∩RNA {len(savi_rna)}")
    genes=list(IFN.keys())
    print("== ★ POSITIVE CONTROL: SAVI carrier ISG (expect HIGH + coordinated across all 6) ==")
    sav=zz[zz.research_id.isin(savi_rna)].sort_values('ifn',ascending=False)
    if len(sav):
        print(sav[['research_id','ifn']+genes].round(2).to_string(index=False))
        S['savi_isg']={'n':len(sav),'mean':round(float(sav.ifn.mean()),3),'median':round(float(sav.ifn.median()),3),'max':round(float(sav.ifn.max()),3),'min':round(float(sav.ifn.min()),3)}
    else: S['savi_isg']='(no SAVI carriers in the RNA cohort — control not directly powered)'
    print("   SAVI ISG summary:",S['savi_isg'])
    print("== broader anchor: RARE STING1 missense (afr_af<0.1%) -> ISG ==")
    rare_miss=[v for (v,aa,r,af) in miss if af<0.001]
    rmc=carriers(rare_miss) & rna_ids
    b=zz.copy(); b['rare_sting_miss']=b.research_id.isin(rmc).astype(int); b['savi']=b.research_id.isin(savi_rna).astype(int)
    import statsmodels.formula.api as smf
    def ols(f,dd,t):
        try: r=smf.ols(f,data=dd).fit(); return {'beta':round(float(r.params[t]),4),'p':round(float(r.pvalues[t]),4),'n':int(r.nobs),'n_carr':int(dd[t].sum())}
        except Exception as e: return {'error':str(e)[:60]}
    S['rare_STING1_missense_to_ISG']=ols('ifn ~ rare_sting_miss',b,'rare_sting_miss'); S['savi_to_ISG']=ols('ifn ~ savi',b,'savi')
    print("   rare STING1 missense -> ISG:",S['rare_STING1_missense_to_ISG']); print("   SAVI -> ISG:",S['savi_to_ISG'])
    print("\n===== SAVI POSITIVE CONTROL (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1100:])
