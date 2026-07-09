#!/usr/bin/env python
"""AoU v9 -- IRAK3 SIGLEC1-PROTEIN neutral band = the LAST IRAK3 molecular test (ZS). The RNA type-I main effect fell inside
the neutral floor (run 076); the only remaining IRAK3 type-I signal is plasma SIGLEC1 (the type-I-SPECIFIC Olink analyte):
IRAK3-LoF -> SIGLEC1 beta 0.555 p=0.0002 (run 067, n41). Question: does that exceed the FOUNDER/TECHNICAL floor of rare
NEUTRAL variants (synonymous/intron, non-immune, matched to IRAK3's Olink carrier count) through the IDENTICAL main-effect
model (SIGLEC1_z ~ variant + 16PC+age+sex+smoking)? VERDICT: IRAK3 beyond the neutral band's 95th pct -> REAL/above-floor ->
IRAK3->type-I rests on protein; INSIDE -> drop IRAK3->type-I molecular entirely (IRAK3 stays EMR-protective only).
p0.0002 is ~10x more significant than the RNA p0.02, so it may well survive. Standard app. Ends 'run complete'/'run failed'.
"""
import os, io, gzip, subprocess, json
import numpy as np, pandas as pd
import warnings; warnings.filterwarnings("ignore"); np.seterr(all="ignore")
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9")
import glob as _glob
_ON="Olink_10k_aou_v9_normalized.tsv.gz"; _cand=[f"{MNT}/v9/multiomics/proteomics/npx/batch_normalized/{_ON}", f"{MNT}/v9/multiomics/proteomics/normalized/{_ON}"]
OLINK=next((p for p in _cand if os.path.exists(p)), None) or (_glob.glob(f"{MNT}/v9/multiomics/proteomics/**/{_ON}",recursive=True) or [_cand[0]])[0]
MANIFEST=f"{MNT}/v9/multiomics/proteomics/manifest.tsv"
SNP=f"{MNT}/v9/wgs/short_read/snpindel"; ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
IRAK3=('12',66183995,66259622); AFMAX=0.001; NEUT_LO,NEUT_HI=0.002,0.008
NEUTRAL=('synonymous','intron','intergenic','non_coding','upstream','downstream','5_prime','3_prime'); NONNEUTRAL=('missense','stop','frameshift','splice','start_lost','stop_lost','inframe')
WINDOWS=[('1',150000000,155000000),('2',85000000,90000000),('3',50000000,55000000),('4',80000000,85000000),('5',40000000,45000000),('6',40000000,45000000),('7',80000000,85000000),('8',70000000,75000000),('9',90000000,95000000),('10',70000000,75000000),('11',65000000,70000000),('12',50000000,55000000),('14',60000000,65000000),('15',65000000,70000000),('16',60000000,65000000),('17',45000000,50000000),('19',45000000,50000000),('20',35000000,40000000)]
EXCLUDE={'RORC','RORA','HAX1','CAMP','LCN2','LTF','ZBP1','CGAS','MB21D1','TREX1','SAMHD1','ADAR','SP100','SP110','SP140','NMI','C1QA','C1QB','C1QC','C1R','C1S','C2','C3','C4A','C4B','C5','C6','C7','C8A','C8B','C9','CFB','CFH','CFI','CFD','SERPING1','IRAK1','IRAK2','IRAK3','IRAK4','MYD88','STING1','TMEM173'}
IMMUNE_PREFIX=('IFI','IFIT','OAS','MX','GBP','RSAD','USP18','HERC','DDX58','IFIH','ISG','PGLYRP','S100A','S100B','DEFA','DEFB','CXCL','CXCR','CCL','CCR','XCL','CX3C','TNFAIP','TNFSF','TNFRSF','TLR','NLRP','NLRC','NAIP','NOD1','NOD2','CARD','AIM2','IRF','STAT','SOCS','JAK','TRIM','SIGLEC','LY6','BST2','FCGR','FCER','FCRL','HLA','TAP','PSMB8','PSMB9','MASP','FCN','CLEC4','CLEC7','TREM','SLAMF','KLR','KIR','NCR','CTLA','PDCD1','LAG3','HAVCR','TIGIT','LTB','LTA')
_IMM={}
def is_immune(g):
    r=_IMM.get(g)
    if r is None:
        r=(g in EXCLUDE) or (g[:2]=='IL' and len(g)>2 and g[2].isdigit()) or any(g.startswith(p) for p in IMMUNE_PREFIX); _IMM[g]=r
    return r
NTARGET=30; OLINK_CARR=(15,170)
def sh(c): return subprocess.run(['bash','-lc',c],capture_output=True,text=True).stdout
def fnum(x):
    try: return float(x)
    except (TypeError,ValueError): return 0.0
S={}
try:
    # ---- Olink SIGLEC1 ----
    hdr=sh(f'zcat "{OLINK}" 2>/dev/null | head -1').rstrip("\n"); cols=hdr.split("\t"); low=[c.lower() for c in cols]
    def findcol(cands):
        for i,c in enumerate(low):
            if any(k==c for k in cands): return i,cols[i]
        for i,c in enumerate(low):
            if any(k in c for k in cands): return i,cols[i]
        return None,None
    si,sc=findcol(['sampleid','sample_id','research_id','person_id']); ai,ac=findcol(['assay','gene','protein']); vi,vc=findcol(['npx','value','normalized'])
    print(f"Olink cols: sample={sc} assay={ac} value={vc}")
    raw=sh(f'zcat "{OLINK}" 2>/dev/null | awk -F"\\t" \'NR==1 || ${ai+1} ~ /^SIGLEC1$/\'')
    m=pd.read_csv(io.StringIO(raw),sep="\t",dtype=str); m['__npx']=pd.to_numeric(m[vc],errors='coerce'); m['__s']=m[sc].astype(str); m['__a']=m[ac].astype(str)
    m=m[m.__a=='SIGLEC1'].dropna(subset=['__npx']); w=m.pivot_table(index='__s',columns='__a',values='__npx',aggfunc='mean'); w.index=w.index.astype(str)
    man=pd.read_csv(MANIFEST,sep="\t",dtype=str); pidcol=next((c for c in man.columns if c.lower() in ('research_id','person_id','pid')),None)
    ols=set(w.index); sidcol=max(man.columns,key=lambda c:len(set(man[c].dropna().astype(str))&ols))
    bridge=dict(zip(man[sidcol].astype(str),man[pidcol].astype(str))); w=w.rename(index=lambda s: bridge.get(str(s))); w=w[w.index.notna()]; w=w.groupby(level=0).mean(numeric_only=True)
    d=pd.DataFrame(index=w.index); d['SIGLEC1_z']=(w.SIGLEC1-w.SIGLEC1.mean())/w.SIGLEC1.std(); d['research_id']=d.index.astype(str)
    OLK=set(d.research_id); print(f"Olink persons with SIGLEC1: {len(d)}")
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    def carr_multi(vids):
        dd={}
        for i in range(0,len(vids),700):
            chunk=[v for v in vids[i:i+700] if v]
            if not chunk: continue
            inl=",".join("'"+v+"'" for v in chunk)
            for r in bq.query(f"SELECT vid, e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid IN ({inl})"): dd.setdefault(str(r.vid),set()).add(str(r.pid))
        return dd
    import pysam
    with gzip.open(VAT,'rt') as fh: COLS=fh.readline().rstrip("\n").split("\t")
    IX={c:COLS.index(c) for c in ['vid','gene_symbol','gvs_afr_af','LoF','consequence']}
    tbx=pysam.TabixFile(VAT)
    irak3=[l.split("\t")[IX['vid']] for l in tbx.fetch('chr'+IRAK3[0],IRAK3[1],IRAK3[2]) for f in [l.split("\t")] if f[IX['gene_symbol']]=='IRAK3' and f[IX['LoF']]=='HC' and fnum(f[IX['gvs_afr_af']])<AFMAX]
    neut=[]; nseen=set()
    for (ch,s,e) in WINDOWS:
        try: it=tbx.fetch('chr'+ch,int(s),int(e))
        except Exception: continue
        nc=0; nl=0
        for line in it:
            nl+=1
            if nl>25000: break
            f=line.split("\t"); g=f[IX['gene_symbol']]
            if g=='' or is_immune(g): continue
            af=fnum(f[IX['gvs_afr_af']]); vid=f[IX['vid']]; cons=f[IX['consequence']]
            if vid not in nseen and nc<8 and NEUT_LO<=af<=NEUT_HI and any(x in cons for x in NEUTRAL) and not any(x in cons for x in NONNEUTRAL):
                nseen.add(vid); neut.append(vid); nc+=1
            if nc>=8: break
        if len(neut)>=NTARGET+15: break
    CARR=carr_multi(list(set(irak3)|set(neut)))
    def burden(vids):
        out=set()
        for v in vids: out|=CARR.get(v,set())
        return out&OLK
    d['IRAK3_LoF']=d.research_id.isin(burden(irak3)).astype(int)
    print(f"IRAK3-LoF vids {len(irak3)} -> Olink carriers {int(d.IRAK3_LoF.sum())} | neutral candidates {len(neut)}")
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
        smk=set(str(r[0]) for r in bq.query(f"SELECT DISTINCT co.person_id FROM `{PROJ}.{DS}.condition_occurrence` co JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id WHERE c.vocabulary_id LIKE 'ICD%' AND (c.concept_code LIKE 'Z72.0%' OR c.concept_code LIKE 'F17%' OR c.concept_code LIKE '305.1%' OR c.concept_code LIKE 'V15.82%')"))
        d['smoker']=d.research_id.isin(smk).astype(float); COV+=['smoker']
    except Exception: pass
    covs=(" + "+" + ".join(COV)) if COV else ""
    import statsmodels.formula.api as smf
    def main_beta(cset,name='B'):
        dd=d.copy(); dd[name]=dd.research_id.isin(cset).astype(int)
        if int(dd[name].sum())<5: return None
        try: r=smf.ols(f'SIGLEC1_z ~ {name}{covs}',data=dd,missing='drop').fit(); return round(float(r.params[name]),3)
        except Exception: return None
    try:
        rr=smf.ols(f'SIGLEC1_z ~ IRAK3_LoF{covs}',data=d,missing='drop').fit(); real_b=round(float(rr.params['IRAK3_LoF']),3); real_p=round(float(rr.pvalues['IRAK3_LoF']),4)
    except Exception: real_b=None; real_p=None
    print(f"\n★ REAL: IRAK3-LoF -> SIGLEC1 (plasma, type-I-specific) beta={real_b} p={real_p} (Olink carriers={int(d.IRAK3_LoF.sum())})")
    S['real']={'beta':real_b,'p':real_p,'n_carr':int(d.IRAK3_LoF.sum())}
    nband=[]
    for v in neut:
        cs=CARR.get(v,set())&OLK; n=len(cs)
        if not (OLINK_CARR[0]<=n<=OLINK_CARR[1]): continue
        b=main_beta(cs)
        if b is not None: nband.append({'vid':v,'beta':b,'n':n})
        if len(nband)>=NTARGET: break
    nbetas=[x['beta'] for x in nband]; S['neutral_band']=nband
    print(f"\n== ★ NEUTRAL-VARIANT null ({len(nband)} rare non-immune synonymous/intron variants, Olink carriers {OLINK_CARR}): does a NEUTRAL variant move SIGLEC1? = the founder/technical FLOOR ==")
    if nbetas and real_b is not None:
        narr=np.array(nbetas); npct=round(100*float((narr<real_b).mean()),1); napct=round(100*float((np.abs(narr)<abs(real_b)).mean()),1)
        S['neutral_summary']={'n':len(nbetas),'median':round(float(np.median(narr)),3),'p95':round(float(np.percentile(narr,95)),3),'p975_abs':round(float(np.percentile(np.abs(narr),97.5)),3),'max_abs':round(float(np.max(np.abs(narr))),3),'IRAK3_signed_pct':npct,'IRAK3_abs_pct':napct}
        print(f"   neutral floor: median={round(float(np.median(narr)),3)} | 95th(signed)={round(float(np.percentile(narr,95)),3)} | 97.5th(|beta|)={round(float(np.percentile(np.abs(narr),97.5)),3)} | max|beta|={round(float(np.max(np.abs(narr))),3)}")
        print(f"   ★★ REAL IRAK3-LoF SIGLEC1 beta={real_b} = {npct}th percentile (signed) / {napct}th (|beta|) of the NEUTRAL floor")
        v=("REAL / ABOVE FOUNDER-FLOOR (IRAK3-LoF->SIGLEC1 beyond 95th pct of neutral variants -> IRAK3->type-I rests on protein; KEEP)" if npct>=95 else "INSIDE the neutral floor -> SIGLEC1 not distinguishable from founder/technical noise -> DROP IRAK3->type-I molecular")
        S['verdict']=v; print(f"\n== VERDICT: {v} ==")
    else: print("   neutral band empty / real missing"); S['verdict']='empty'
    print("\n== This is the LAST IRAK3 molecular test. Above floor -> IRAK3 keeps a (protein-level) type-I signature; inside -> IRAK3 is EMR-protective ONLY (asthma/ILD, Lipsitch+cross-ancestry). Either way the EMR arm stands. ==")
    print("\n===== IRAK3 SIGLEC1-PROTEIN NEUTRAL BAND (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1400:])
