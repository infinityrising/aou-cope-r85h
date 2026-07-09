#!/usr/bin/env python
"""AoU v9 -- NEGATIVE-CONTROL-variant PAIR BAND (prereg's decisive specificity de-confounder; task#4, never built until now).
Question: is the CENTERPIECE interaction -- R85H(rare) x cis-AQ(common) -> TYPE-I-SPECIFIC ISG (beta +0.77, p0.01) --
SPECIFIC, or just what any rare x common variant PAIR produces on the type-I score (a founder/structure/random-pair artifact)?
Build a NULL BAND: ~25 pairs of (rare AFR-AF-matched neutral variant) x (common neutral variant), NEITHER in an IFN/immune
gene, each pushed through the IDENTICAL phased-RNA interaction model (typeI_specific ~ V_rare * V_common + 16PC+age+sex+cell
fractions) on the SAME RNA∩PAV cohort. VERDICT: R85H:cisAQ beta beyond the band's 95th percentile -> SPECIFIC; inside -> the
interaction is NOT distinguishable from random neutral pairs -> artifact. This is the molecular analogue of the single-variant
band and the de-confounder the prereg promised. Standard app. Ends 'run complete' / 'run failed'.
"""
import os, io, gzip, subprocess, json
import numpy as np, pandas as pd
import warnings; warnings.filterwarnings("ignore"); np.seterr(all="ignore")
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); RNADIR=f"{MNT}/v9/multiomics/rnaseq"; SNP=f"{MNT}/v9/wgs/short_read/snpindel"
TPM=f"{RNADIR}/rsem/aou_rnaseq_20260415.rsem_genes_tpm.txt.gz"; ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CELL=os.path.expanduser("~/cell_fractions_v9.csv"); PHPAV=os.path.expanduser("~/sting_phenome_pav_v9.csv")
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
R85H_VID='19-18911007-C-T'
TYPEI={'SIGLEC1':'ENSG00000088827','IFI27':'ENSG00000165949','USP18':'ENSG00000184979','IFI6':'ENSG00000126709','IFI44L':'ENSG00000137959'}
# control-variant selection (matched to R85H rare / cis-AQ common; neutral; non-immune)
RARE_LO,RARE_HI=0.008,0.02        # match R85H AFR AF ~1.3%
COMMON_LO,COMMON_HI=0.03,0.12     # common (cis-AQ-like); enough carriers for an interaction
EUR_MAX=0.05
NEUTRAL=('synonymous','intron','intergenic','non_coding','upstream','downstream','5_prime_utr','3_prime_utr','5_prime_UTR','3_prime_UTR')
NONNEUTRAL=('missense','stop','frameshift','splice','start_lost','stop_lost','inframe')
WINDOWS=[('1',60000000,66000000),('2',60000000,66000000),('3',70000000,76000000),('4',60000000,66000000),('5',40000000,46000000),('6',60000000,66000000),('7',62000000,68000000),('8',60000000,66000000),('9',80000000,86000000),('10',60000000,66000000),('11',70000000,76000000),('12',40000000,46000000),('13',50000000,56000000),('14',50000000,56000000),('15',60000000,66000000),('16',50000000,56000000),('17',40000000,46000000),('18',30000000,36000000),('19',40000000,46000000),('20',30000000,36000000),('21',30000000,36000000),('22',30000000,36000000)]
EXCLUDE={'COPA','COPB1','COPB2','COPG1','COPG2','COPZ1','COPZ2','COPE','ARCN1','STING1','TMEM173','IRAK1','IRAK3','IRAK4','MYD88','IFIH1','DDX58','TBK1','IRF3','IRF7'}
NTARGET=45; RNA_RARE=(8,250); RNA_COMMON=(200,3000); MIN_DOUBLE=3   # widened: more windows + higher target + higher common-carrier floor -> ~30+ pairs above the >=3-double threshold for a tighter empirical p
def sh(c): return subprocess.run(['bash','-lc',c],capture_output=True,text=True).stdout
def fnum(x):
    try: return float(x)
    except (TypeError,ValueError): return 0.0
S={}
try:
    # ---- type-I-specific RNA score ----
    pat="|".join(TYPEI.values()); raw=sh(f'zcat "{TPM}" 2>/dev/null | awk -F"\\t" \'NR==1 || $1 ~ /{pat}/\'')
    m=pd.read_csv(io.StringIO(raw),sep="\t"); m['ensg']=m.gene_id.str.split('.').str[0]
    m=m[m.ensg.isin(TYPEI.values())].drop_duplicates('ensg').set_index('ensg')
    dropc=[c for c in ['gene_id','transcript_id(s)'] if c in m.columns]; ex=m.drop(columns=dropc).T.astype(float)
    ex.index=ex.index.astype(str); lg=np.log2(ex+1); z=(lg-lg.mean())/lg.std()
    got=[g for g in TYPEI.values() if g in set(m.index)]
    d=pd.DataFrame(index=z.index); d['typeI_specific']=z[got].mean(axis=1); d['research_id']=d.index.astype(str)
    print(f"type-I-specific genes found: {len(got)}/{len(TYPEI)} | RNA samples {len(d)}")
    pav=pd.read_csv(PHPAV); pav['research_id']=pav.research_id.astype(str)
    d=d.merge(pav[['research_id','AQ_d','HAQ_d']],on='research_id',how='inner'); d['cisAQ']=(d.AQ_d>=1).astype(int); d['cisHAQ']=(d.HAQ_d>=1).astype(int)
    RNA=set(d.research_id); print(f"RNA∩PAV {len(d)} | cisAQ {int(d.cisAQ.sum())} | cisHAQ {int(d.cisHAQ.sum())}")
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    def carr(vid): return set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid='{vid}'"))
    d['R85H']=d.research_id.isin(carr(R85H_VID)).astype(int)
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
    def inter_beta(rare_set,common_set,rname='Vr',cname='Vc'):
        dd=d.copy(); dd[rname]=dd.research_id.isin(rare_set).astype(int); dd[cname]=dd.research_id.isin(common_set).astype(int)
        nd=int(((dd[rname]==1)&(dd[cname]==1)).sum())
        if int(dd[rname].sum())<5 or int(dd[cname].sum())<20 or nd<MIN_DOUBLE: return None,nd
        try:
            r=smf.ols(f'typeI_specific ~ {rname}*{cname}{covs}',data=dd,missing='drop').fit()
            t=f'{rname}:{cname}'; return (round(float(r.params[t]),3) if t in r.params else None),nd
        except Exception: return None,nd
    # ---- REAL interaction ----
    real_b,real_nd=inter_beta(carr(R85H_VID),set(d[d.cisAQ==1].research_id),'R85H2','cisAQ2')
    # use the canonical fit for the reference (R85H*cisAQ directly)
    try:
        rr=smf.ols(f'typeI_specific ~ R85H*cisAQ{covs}',data=d,missing='drop').fit(); real_b=round(float(rr.params['R85H:cisAQ']),3); real_p=round(float(rr.pvalues['R85H:cisAQ']),4)
    except Exception: real_p=None
    try:
        rh=smf.ols(f'typeI_specific ~ R85H*cisHAQ{covs}',data=d,missing='drop').fit(); realh_b=round(float(rh.params['R85H:cisHAQ']),3); realh_p=round(float(rh.pvalues['R85H:cisHAQ']),4)
    except Exception: realh_b=None; realh_p=None
    print(f"\n★ REAL: R85H:cisAQ -> typeI_specific beta={real_b} p={real_p} (doubles={int(((d.R85H==1)&(d.cisAQ==1)).sum())}) | NEG-CTRL R85H:cisHAQ beta={realh_b} p={realh_p} (doubles={int(((d.R85H==1)&(d.cisHAQ==1)).sum())})")
    S['real']={'beta':real_b,'p':real_p}; S['real_HAQ']={'beta':realh_b,'p':realh_p}
    # ---- sample matched neutral control variants ----
    with gzip.open(VAT,'rt') as fh: COLS=fh.readline().rstrip("\n").split("\t")
    has_eur='gvs_eur_af' in COLS
    IX={c:COLS.index(c) for c in (['vid','gene_symbol','consequence','gvs_afr_af']+(['gvs_eur_af'] if has_eur else []))}
    import pysam; tbx=pysam.TabixFile(VAT); rare_cand=[]; common_cand=[]; seen=set()
    for (ch,s,e) in WINDOWS:
        try: it=tbx.fetch('chr'+ch,s,e)
        except Exception: continue
        rn=cn=0
        for line in it:
            f=line.split("\t"); cons=f[IX['consequence']]; vid=f[IX['vid']]
            if vid in seen or any(x in cons for x in NONNEUTRAL) or not any(x in cons for x in NEUTRAL): continue
            if f[IX['gene_symbol']] in EXCLUDE: continue
            af=fnum(f[IX['gvs_afr_af']])
            if has_eur and fnum(f[IX['gvs_eur_af']])>EUR_MAX: continue
            if RARE_LO<=af<=RARE_HI and rn<6: seen.add(vid); rare_cand.append(vid); rn+=1
            elif COMMON_LO<=af<=COMMON_HI and cn<6: seen.add(vid); common_cand.append(vid); cn+=1
            if rn>=6 and cn>=6: break     # stop scanning this window once full (was scanning the whole 6Mb -> the hang)
    print(f"candidate neutral controls: rare {len(rare_cand)} | common {len(common_cand)} (BATCHED carrier lookup + RNA-carrier match)")
    # BATCHED carrier lookup: all candidates in ~1-2 queries (was 1 BQ/candidate -> the slowdown)
    def carr_multi(vids):
        dd={}
        for i in range(0,len(vids),700):
            chunk=vids[i:i+700]
            if not chunk: continue
            inl=",".join("'"+v+"'" for v in chunk)
            for r in bq.query(f"SELECT vid, e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid IN ({inl})"): dd.setdefault(str(r.vid),set()).add(str(r.pid))
        return dd
    CARR=carr_multi(rare_cand+common_cand)
    def screen(cands,lo,hi,target):
        out=[]
        for v in cands:
            cs=CARR.get(v,set())&RNA; n=len(cs)
            if lo<=n<=hi: out.append((v,cs))
            if len(out)>=target: break
        return out
    rare=screen(rare_cand,*RNA_RARE,NTARGET); common=screen(common_cand,*RNA_COMMON,NTARGET)
    print(f"matched controls in RNA∩PAV: rare {len(rare)} (carriers {RNA_RARE}) | common {len(common)} (carriers {RNA_COMMON})")
    S['n_rare']=len(rare); S['n_common']=len(common)
    # ---- NULL BAND of rare x common interaction betas ----
    npair=min(len(rare),len(common)); band=[]
    for i in range(npair):
        rv,rs=rare[i]; cv,cs=common[i]; b,nd=inter_beta(rs,cs)
        if b is not None: band.append({'rare':rv,'common':cv,'beta':b,'double':nd})
    S['band']=band; betas=[x['beta'] for x in band]
    print(f"\n== NULL BAND: {len(band)} rare x common neutral pairs through the identical typeI_specific interaction model ==")
    if betas:
        arr=np.array(betas); pct=round(100*float((arr<real_b).mean()),1); apct=round(100*float((np.abs(arr)<abs(real_b)).mean()),1)
        S['band_summary']={'n':len(betas),'median':round(float(np.median(arr)),3),'p95':round(float(np.percentile(arr,95)),3),'p975_abs':round(float(np.percentile(np.abs(arr),97.5)),3),'max_abs':round(float(np.max(np.abs(arr))),3),'real_signed_pct':pct,'real_abs_pct':apct}
        print(f"   band betas: median={round(float(np.median(arr)),3)} | 95th(signed)={round(float(np.percentile(arr,95)),3)} | 97.5th(|beta|)={round(float(np.percentile(np.abs(arr),97.5)),3)} | max|beta|={round(float(np.max(np.abs(arr))),3)}")
        print(f"   ★ REAL R85H:cisAQ beta={real_b} sits at the {pct}th percentile (signed) / {apct}th (|beta|, two-sided) of the null pair band")
        if realh_b is not None:
            hpct=round(100*float((arr<realh_b).mean()),1); hapct=round(100*float((np.abs(arr)<abs(realh_b)).mean()),1)
            S['band_summary']['HAQ_signed_pct']=hpct; S['band_summary']['HAQ_abs_pct']=hapct
            print(f"   ☆ NEG-CTRL R85H:cisHAQ beta={realh_b} sits at the {hpct}th percentile (signed) / {hapct}th (|beta|) -> EXPECT INSIDE the band (null), the built-in HAQ specificity contrast to cisAQ")
        verdict=("SPECIFIC (interaction beyond 95th pct of the neutral-pair null -> not a random rare×common / founder artifact)" if (real_b is not None and pct>=95) else "INSIDE THE BAND -> the interaction is NOT distinguishable from random neutral variant pairs -> weakens the centerpiece / declare non-specific")
        S['verdict']=verdict; print(f"\n== VERDICT: {verdict} ==")
    else: print("   band empty (no pairs converged) — widen AF windows / carrier ranges"); S['verdict']='band empty'
    print("\n== NOTE: this de-confounds SPECIFICITY (is the pair special?); it does NOT fix the n=11 POWER fragility (LOO/rank) -- both caveats stand together. ==")
    print("\n===== NEG-CONTROL PAIR BAND (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1400:])
