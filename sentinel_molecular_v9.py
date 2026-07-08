#!/usr/bin/env python
"""AoU v9 — SENTINEL (1348881) full molecular characterization = the R85H+AQ hub-branch worked example.
 (1) Is she in the long-read (PAV) / RNA cohorts? (2) PHASE her AQ — hap1/hap2: are G230A+R293Q in CIS (true AQ-strict)?
 (3) Her ISG (6-gene z + cohort percentile) IF in RNA — N-of-1: does the AQ-route COPA patient have ELEVATED interferon?
 (4) R85H dose, (5) COPA-spectrum EHR status, (6) confirm her documented COPI-γ modifiers (COPG1 I55T, COPZ2 G22R).
INTERNAL n=1 — suppress on export. STANDARD app. Ends 'run complete' / 'run failed'.
"""
import os, io, gzip, subprocess, json, re
import numpy as np, pandas as pd
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); SNP=f"{MNT}/v9/wgs/short_read/snpindel"; RNADIR=f"{MNT}/v9/multiomics/rnaseq"
TPM=f"{RNADIR}/rsem/aou_rnaseq_20260415.rsem_genes_tpm.txt.gz"; ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1); SENT='1348881'; R85H_VID='19-18911007-C-T'
RNAPAV=os.path.expanduser("~/copi_sting_pav_v9.csv"); PHPAV=os.path.expanduser("~/sting_phenome_pav_v9.csv")
IFN={'IFI27':'ENSG00000165949','IFI44L':'ENSG00000137959','IFIT1':'ENSG00000185745','ISG15':'ENSG00000187608','RSAD2':'ENSG00000134321','SIGLEC1':'ENSG00000088827'}
STINGV={'R71H':'5-139481493-C-T','G230A':'5-139478340-C-G','R293Q':'5-139477397-C-T','R220H':'5-139478370-C-T'}
COPG1=('3',129249575,129278068); COPZ2=('17',48026142,48038030)
PHENO={'ILD':['J84','515','516','J98.4'],'infl_arthritis':['M05','M06','M07','M08','714'],'vasculitis':['M31.3','M31.7','M30','M31','446','I77.6']}
S={'sentinel':SENT}
def sh(c): return subprocess.run(['bash','-lc',c],capture_output=True,text=True).stdout
try:
    print("== (1)+(2) cohort membership + PHASED STING (AQ cis?) ==")
    for path,lab in [(RNAPAV,'RNA_PAV'),(PHPAV,'phenome_PAV')]:
        try:
            g=pd.read_csv(path); g['research_id']=g.research_id.astype(str); row=g[g.research_id==SENT]
            if len(row):
                r=row.iloc[0]; S[lab]={'present':True,'hap1':r.hap1,'hap2':r.hap2,'AQ_cis':bool(r.hap1=='AQ' or r.hap2=='AQ'),'HAQ_cis':bool(r.hap1=='HAQ' or r.hap2=='HAQ'),'R85H_dose':int(r.R85H_d) if 'R85H_d' in g.columns else 'n/a'}
            else: S[lab]={'present':False,'cohort_n':len(g)}
        except FileNotFoundError: S[lab]='(csv not found)'
        print(f"   {lab}: {S[lab]}")
    print("== (3) her ISG (6-gene z + percentile), if in RNA ==")
    raw=sh(f'zcat "{TPM}" 2>/dev/null | awk -F"\\t" \'NR==1 || $1 ~ /{"|".join(IFN.values())}/\'')
    m=pd.read_csv(io.StringIO(raw),sep="\t"); m['ensg']=m.gene_id.str.split('.').str[0]
    m=m[m.ensg.isin(IFN.values())].set_index('ensg').drop(columns=['gene_id','transcript_id(s)'])
    ex=m.T.astype(float); ex.index=ex.index.astype(str); lg=np.log2(ex+1); z=(lg-lg.mean())/lg.std(); isg=z.mean(axis=1)
    if SENT in isg.index:
        val=float(isg.loc[SENT]); pct=float((isg<val).mean()*100)
        inv={v:k for k,v in IFN.items()}; pg={inv[c]:round(float(z.loc[SENT,c]),2) for c in z.columns}
        S['sentinel_ISG']={'z':round(val,3),'percentile':round(pct,1),'per_gene_z':pg}
        print(f"   ★ sentinel ISG z={val:.3f} (percentile {pct:.1f}) per-gene {pg}")
    else: S['sentinel_ISG']='(sentinel not in RNA-seq cohort)'; print("   "+S['sentinel_ISG'])
    print("== (4)+(6) BQ carrier status: R85H, STING, COPI-γ modifiers ==")
    import pysam
    with gzip.open(VAT,'rt') as fh: COLS=fh.readline().rstrip("\n").split("\t")
    hasaa='aa_change' in COLS; IX={c:COLS.index(c) for c in ['vid','gene_symbol','aa_change'] if c in COLS}
    tbx=pysam.TabixFile(VAT); mods={}
    for gene,(ch,s,e),pat in [('COPG1',COPG1,r'Ile55Thr'),('COPZ2',COPZ2,r'Gly22Arg')]:
        for line in tbx.fetch('chr'+ch,s,e):
            f=line.split("\t")
            if f[IX['gene_symbol']]!=gene: continue
            if hasaa and re.search(pat,f[IX['aa_change']]): mods[f"{gene}_{pat}"]=f[IX['vid']]
    print(f"   found modifier vids: {mods}")
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    allv=[R85H_VID]+list(STINGV.values())+list(mods.values())
    inl=",".join("'"+v+"'" for v in allv)
    hits=set(r.vid for r in bq.query(f"SELECT DISTINCT vid FROM {T}, UNNEST(person_ids.list) e WHERE vid IN ({inl}) AND e.element={SENT}"))
    S['R85H']=R85H_VID in hits; S['STING']={k:(v in hits) for k,v in STINGV.items()}; S['COPI_modifiers']={k:(v in hits) for k,v in mods.items()}
    print(f"   R85H {S['R85H']} | STING {S['STING']} | COPI-γ {S['COPI_modifiers']}")
    print("== (5) COPA-spectrum EHR ==")
    codes=[c for v in PHENO.values() for c in v]; lk=" OR ".join([f"c.concept_code LIKE '{p}%'" for p in codes])
    dom={}
    for k,cc in PHENO.items():
        lkk=" OR ".join([f"c.concept_code LIKE '{p}%'" for p in cc])
        n=len([r for r in bq.query(f"SELECT 1 FROM `{PROJ}.{DS}.condition_occurrence` co JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id WHERE co.person_id={SENT} AND c.vocabulary_id LIKE 'ICD%' AND ({lkk}) LIMIT 1")])
        dom[k]=n>0
    S['COPA_domains']=dom; print(f"   COPA-spectrum: {dom}")
    print("\n== SUMMARY: sentinel = R85H + AQ-STING hub branch ==")
    print(f"   R85H+ / AQ(G230A+R293Q,no R71H) / COPI-γ(I55T,G22R) — the STING/COPI route; NO IRAK3-LoF (run 035).")
    print("\n===== SENTINEL MOLECULAR (paste back; INTERNAL n=1) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1100:])
