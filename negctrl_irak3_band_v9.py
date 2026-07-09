#!/usr/bin/env python
"""AoU v9 -- IRAK3 co-route SPECIFICITY band (the IRAK3 analogue of the R85H×AQ pair band; ZS: 'what about IRAK3?').
IRAK3's RNA finding is a MAIN effect (IRAK3-LoF rare-truncating burden -> TYPE-I-SPECIFIC ISG, beta~0.32 p0.02) -- not an
interaction, so it can't go in the rare×common pair band, and R85H×IRAK3 is untestable (~2-4 doubles). Question here: is
IRAK3-LoF's type-I elevation SPECIFIC, or would ANY random rare-LoF GENE burden raise the type-I score similarly (a generic
'rare-LoF-burden' or founder artifact)? Build a NULL BAND of ~25 random NON-IMMUNE genes' rare-LoF(HC) burdens, matched to
IRAK3-LoF's carrier frequency, each pushed through the IDENTICAL main-effect model (typeI_specific ~ burden + 16PC+age+sex+
cell-fractions) on the SAME RNA∩PAV cohort. VERDICT: IRAK3-LoF beta beyond the band's 95th percentile -> SPECIFIC; inside ->
not distinguishable from a random rare-LoF gene burden. Complements the run-061 Lipsitch (EMR) specificity. Ends 'run complete'/'run failed'.
"""
import os, io, gzip, subprocess, json
import numpy as np, pandas as pd
import warnings; warnings.filterwarnings("ignore"); np.seterr(all="ignore")
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); RNADIR=f"{MNT}/v9/multiomics/rnaseq"; SNP=f"{MNT}/v9/wgs/short_read/snpindel"
TPM=f"{RNADIR}/rsem/aou_rnaseq_20260415.rsem_genes_tpm.txt.gz"; ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CELL=os.path.expanduser("~/cell_fractions_v9.csv"); PHPAV=os.path.expanduser("~/sting_phenome_pav_v9.csv")
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
TYPEI={'SIGLEC1':'ENSG00000088827','IFI27':'ENSG00000165949','USP18':'ENSG00000184979','IFI6':'ENSG00000126709','IFI44L':'ENSG00000137959'}
IRAK3=('12',66183995,66259622)
AFMAX=0.001                        # rare (matches IRAK3-LoF selection)
# broad windows to harvest genes carrying rare LoF-HC variants (avoid the IFN/immune/COPI loci via EXCLUDE)
WINDOWS=[('1',150000000,155000000),('2',85000000,90000000),('3',50000000,55000000),('4',80000000,85000000),('5',50000000,55000000),('6',40000000,45000000),('7',80000000,85000000),('8',70000000,75000000),('9',90000000,95000000),('10',70000000,75000000),('11',65000000,70000000),('12',50000000,55000000),('14',60000000,65000000),('15',65000000,70000000),('16',60000000,65000000),('17',45000000,50000000),('19',45000000,50000000),('20',35000000,40000000)]  # ~5Mb each + outer break -> bounded scan
EXCLUDE={'COPA','COPB1','COPB2','COPG1','COPG2','COPZ1','COPZ2','COPE','ARCN1','STING1','TMEM173','IRAK1','IRAK2','IRAK3','IRAK4','MYD88','IFIH1','DDX58','TBK1','IRF3','IRF7','TLR7','TLR9','TICAM1','TRAF3','TRAF6','RORC','RORA','HAX1','CAMP','LCN2','LTF','ZBP1','CGAS','MB21D1','TREX1','SAMHD1','ADAR','SP100','SP110','SP140','NMI','C1QA','C1QB','C1QC','C1R','C1S','C2','C3','C4A','C4B','C5','C6','C7','C8A','C8B','C9','CFB','CFH','CFI','CFD','SERPING1'}
# CLEANED null: exclude ANY immune/inflammation/antimicrobial/ISG gene (run-074 band was contaminated by PGLYRP/RORC/HAX1/S100)
IMMUNE_PREFIX=('IFI','IFIT','OAS','MX','GBP','RSAD','USP18','HERC','DDX58','IFIH','ISG','PGLYRP','S100A','S100B','DEFA','DEFB','CXCL','CXCR','CCL','CCR','XCL','CX3C','TNFAIP','TNFSF','TNFRSF','TLR','NLRP','NLRC','NAIP','NOD1','NOD2','CARD','AIM2','IRF','STAT','SOCS','JAK','TRIM','SIGLEC','LY6','BST2','FCGR','FCER','FCRL','HLA','TAP','PSMB8','PSMB9','MASP','FCN','CLEC4','CLEC7','TREM','SLAMF','KLR','KIR','NCR','CTLA','PDCD1','LAG3','HAVCR','TIGIT','LTB','LTA')
def is_immune(g):
    if g in EXCLUDE: return True
    if g[:2]=='IL' and len(g)>2 and g[2].isdigit(): return True   # interleukins (not ILK/ILDR/ILVBL)
    return any(g.startswith(p) for p in IMMUNE_PREFIX)
NTARGET=40; RNA_CARR=(20,150); MIN_LOF_VIDS=2   # tighter carrier match to IRAK3(35) -> less-noisy burden betas; more genes for a stable percentile
def sh(c): return subprocess.run(['bash','-lc',c],capture_output=True,text=True).stdout
def fnum(x):
    try: return float(x)
    except (TypeError,ValueError): return 0.0
S={}
try:
    pat="|".join(TYPEI.values()); raw=sh(f'zcat "{TPM}" 2>/dev/null | awk -F"\\t" \'NR==1 || $1 ~ /{pat}/\'')
    m=pd.read_csv(io.StringIO(raw),sep="\t"); m['ensg']=m.gene_id.str.split('.').str[0]
    m=m[m.ensg.isin(TYPEI.values())].drop_duplicates('ensg').set_index('ensg')
    dropc=[c for c in ['gene_id','transcript_id(s)'] if c in m.columns]; ex=m.drop(columns=dropc).T.astype(float)
    ex.index=ex.index.astype(str); lg=np.log2(ex+1); z=(lg-lg.mean())/lg.std(); got=[g for g in TYPEI.values() if g in set(m.index)]
    d=pd.DataFrame(index=z.index); d['typeI_specific']=z[got].mean(axis=1); d['research_id']=d.index.astype(str)
    pav=pd.read_csv(PHPAV); pav['research_id']=pav.research_id.astype(str); d=d.merge(pav[['research_id']],on='research_id',how='inner')
    RNA=set(d.research_id); print(f"type-I genes {len(got)}/{len(TYPEI)} | RNA∩PAV {len(d)}")
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
    NEUTRAL=('synonymous','intron','intergenic','non_coding','upstream','downstream','5_prime','3_prime'); NONNEUTRAL=('missense','stop','frameshift','splice','start_lost','stop_lost','inframe')
    tbx=pysam.TabixFile(VAT)
    # ---- IRAK3-LoF (the real burden) ----
    irak3=[l.split("\t")[IX['vid']] for l in tbx.fetch('chr'+IRAK3[0],IRAK3[1],IRAK3[2]) for f in [l.split("\t")] if f[IX['gene_symbol']]=='IRAK3' and f[IX['LoF']]=='HC' and fnum(f[IX['gvs_afr_af']])<AFMAX]
    # ---- harvest random NON-IMMUNE genes with rare LoF-HC burdens ----
    genevids={}; neut=[]; nseen=set()
    for (ch,s,e) in WINDOWS:
        try: it=tbx.fetch('chr'+ch,int(s),int(e))
        except Exception: continue
        nc=0
        for line in it:
            f=line.split("\t"); g=f[IX['gene_symbol']]
            if g=='' or is_immune(g): continue
            af=fnum(f[IX['gvs_afr_af']]); vid=f[IX['vid']]; cons=f[IX['consequence']]
            if f[IX['LoF']]=='HC' and af<AFMAX: genevids.setdefault(g,set()).add(vid)                                     # LoF-gene-burden null (gene-specificity)
            elif vid not in nseen and nc<8 and 0.002<=af<=0.008 and any(x in cons for x in NEUTRAL) and not any(x in cons for x in NONNEUTRAL):
                nseen.add(vid); neut.append(vid); nc+=1                                                                   # NEUTRAL-variant null (founder/technical floor)
        if sum(1 for vs in genevids.values() if len(vs)>=MIN_LOF_VIDS)>=NTARGET*2 and len(neut)>=NTARGET*3: break
    genes=[g for g,vs in genevids.items() if len(vs)>=MIN_LOF_VIDS]
    print(f"IRAK3-LoF vids {len(irak3)} | candidate non-immune LoF-burden genes (>= {MIN_LOF_VIDS} rare-LoF vids): {len(genes)}")
    # batched carriers for IRAK3 + all candidate gene burdens
    allv=list(set(irak3)|{v for g in genes for v in genevids[g]}|set(neut)); CARR=carr_multi(allv)
    def burden(vids):
        out=set()
        for v in vids: out|=CARR.get(v,set())
        return out&RNA
    d['IRAK3_LoF']=d.research_id.isin(burden(irak3)).astype(int)
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
    except Exception: print("   (cell_fractions not found)")
    covs=(" + "+" + ".join(COV)) if COV else ""
    import statsmodels.formula.api as smf
    def main_beta(cset,name='B'):
        dd=d.copy(); dd[name]=dd.research_id.isin(cset).astype(int)
        if int(dd[name].sum())<5: return None
        try:
            r=smf.ols(f'typeI_specific ~ {name}{covs}',data=dd,missing='drop').fit(); return round(float(r.params[name]),3)
        except Exception: return None
    # real
    try:
        rr=smf.ols(f'typeI_specific ~ IRAK3_LoF{covs}',data=d,missing='drop').fit(); real_b=round(float(rr.params['IRAK3_LoF']),3); real_p=round(float(rr.pvalues['IRAK3_LoF']),4)
    except Exception: real_b=None; real_p=None
    print(f"\n★ REAL: IRAK3-LoF -> typeI_specific beta={real_b} p={real_p} (RNA carriers={int(d.IRAK3_LoF.sum())})")
    S['real']={'beta':real_b,'p':real_p,'n_carr':int(d.IRAK3_LoF.sum())}
    # ---- NULL BAND: random gene LoF burdens matched to IRAK3 RNA-carrier count ----
    band=[]
    for g in genes:
        cs=burden(genevids[g]); n=len(cs)
        if not (RNA_CARR[0]<=n<=RNA_CARR[1]): continue
        b=main_beta(cs)
        if b is not None: band.append({'gene':g,'beta':b,'n':n,'n_lof_vids':len(genevids[g])})
        if len(band)>=NTARGET: break
    S['band']=band; betas=[x['beta'] for x in band]
    print(f"\n== NULL BAND: {len(band)} random non-immune rare-LoF GENE burdens (RNA carriers {RNA_CARR}, matched to IRAK3-LoF {int(d.IRAK3_LoF.sum())}) through the identical typeI_specific main-effect model ==")
    if betas and real_b is not None:
        arr=np.array(betas); pct=round(100*float((arr<real_b).mean()),1); apct=round(100*float((np.abs(arr)<abs(real_b)).mean()),1)
        S['band_summary']={'n':len(betas),'median':round(float(np.median(arr)),3),'p95':round(float(np.percentile(arr,95)),3),'p975_abs':round(float(np.percentile(np.abs(arr),97.5)),3),'max_abs':round(float(np.max(np.abs(arr))),3),'IRAK3_signed_pct':pct,'IRAK3_abs_pct':apct}
        print(f"   band betas: median={round(float(np.median(arr)),3)} | 95th(signed)={round(float(np.percentile(arr,95)),3)} | 97.5th(|beta|)={round(float(np.percentile(np.abs(arr),97.5)),3)} | max|beta|={round(float(np.max(np.abs(arr))),3)}")
        print(f"   ★ REAL IRAK3-LoF beta={real_b} sits at the {pct}th percentile (signed) / {apct}th (|beta|, two-sided) of the random-gene-LoF-burden null band")
        verdict=("SPECIFIC (IRAK3-LoF -> type-I beyond 95th pct of random rare-LoF gene burdens -> not a generic LoF-burden/founder effect)" if pct>=95 else "INSIDE THE BAND -> IRAK3-LoF's type-I elevation is NOT distinguishable from a random rare-LoF gene burden -> non-specific")
        S['verdict']=verdict; print(f"\n== VERDICT: {verdict} ==")
    else: print("   band empty (no gene burdens matched) — widen windows / carrier range"); S['verdict']='band empty'
    # ==== NEUTRAL-VARIANT null (PRIMARY founder/technical de-confounder; ZS: synonymous/intron like the pair band) ====
    nband=[]
    for v in neut:
        cs=CARR.get(v,set())&RNA; n=len(cs)
        if not (RNA_CARR[0]<=n<=RNA_CARR[1]): continue
        b=main_beta(cs)
        if b is not None: nband.append({'vid':v,'beta':b,'n':n})
        if len(nband)>=NTARGET: break
    nbetas=[x['beta'] for x in nband]; S['neutral_band']=nband
    print(f"\n== ★ NEUTRAL-VARIANT null ({len(nband)} rare non-immune synonymous/intron/intergenic variants, RNA carriers {RNA_CARR}): does a FUNCTIONALLY-NEUTRAL variant show a type-I association? = the founder/technical FLOOR ==")
    if nbetas and real_b is not None:
        narr=np.array(nbetas); npct=round(100*float((narr<real_b).mean()),1); napct=round(100*float((np.abs(narr)<abs(real_b)).mean()),1)
        S['neutral_summary']={'n':len(nbetas),'median':round(float(np.median(narr)),3),'p95':round(float(np.percentile(narr,95)),3),'p975_abs':round(float(np.percentile(np.abs(narr),97.5)),3),'IRAK3_signed_pct':npct,'IRAK3_abs_pct':napct}
        print(f"   neutral floor: median={round(float(np.median(narr)),3)} | 95th(signed)={round(float(np.percentile(narr,95)),3)} | 97.5th(|beta|)={round(float(np.percentile(np.abs(narr),97.5)),3)}")
        print(f"   ★★ REAL IRAK3-LoF beta={real_b} = {npct}th percentile (signed) / {napct}th (|beta|) of the NEUTRAL floor")
        nverdict=("REAL / ABOVE FOUNDER-FLOOR (IRAK3-LoF→type-I beyond 95th pct of functionally-NEUTRAL variants → NOT a founder/technical artifact)" if npct>=95 else "within the neutral floor → not distinguishable from founder/technical noise")
        S['neutral_verdict']=nverdict; print(f"   == NEUTRAL-NULL VERDICT: {nverdict} ==")
    else: print("   neutral band empty")
    print("\n== TWO NULLS: (1) NEUTRAL-variant floor = is the effect REAL vs founder/technical noise (PRIMARY; parallels the pair band + run-061 Lipsitch synonymous-null). (2) LoF-gene-burden = is IRAK3 SPECIAL among knockouts (harder; generic-LoF-sensitive). Real+above-neutral is the load-bearing claim; gene-specificity is the bonus. ==")
    print("\n===== IRAK3 BANDS (neutral + LoF-gene) (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1400:])
