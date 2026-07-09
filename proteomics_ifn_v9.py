#!/usr/bin/env python
"""AoU v9 -- PROTEOMICS squeeze: Olink 10k plasma protein-level corroboration of the RNA IFN signal (run 044). Tests
whether R85H / AQ / IRAK3-LoF / R85H×AQ elevate IFN-inducible + inflammatory PLASMA proteins -- CXCL9/CXCL10/CXCL11 =
the plasma readout of type-I/II ISG; IL6/TNF/IL18/etc = NF-kB/inflammatory. A SECOND modality for the conditional-modifier
claim; the key test is whether R85H×AQ elevates the IFN-protein composite as it did the RNA type-I ISG module. Auto-detects
the Olink normalized-TSV structure (long format: SampleID / Assay / NPX). Reports coverage + carrier counts (so an
ID-space mismatch is visible). Standard app. Ends 'run complete' / 'run failed'.
"""
import os, io, gzip, subprocess, json
import numpy as np, pandas as pd
import warnings; warnings.filterwarnings("ignore"); np.seterr(all="ignore")   # small-cell fit warnings are expected/benign; keep output readable
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9")
import glob as _glob
_ON="Olink_10k_aou_v9_normalized.tsv.gz"; _cand=[f"{MNT}/v9/multiomics/proteomics/npx/batch_normalized/{_ON}", f"{MNT}/v9/multiomics/proteomics/normalized/{_ON}"]
OLINK=next((p for p in _cand if os.path.exists(p)), None) or (_glob.glob(f"{MNT}/v9/multiomics/proteomics/**/{_ON}",recursive=True) or [_cand[0]])[0]  # self-locating: subdir was renamed normalized/ -> npx/batch_normalized/
SNP=f"{MNT}/v9/wgs/short_read/snpindel"; ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
R85H_VID='19-18911007-C-T'; STINGV={'R71H':'5-139481493-C-T','G230A':'5-139478340-C-G','R293Q':'5-139477397-C-T'}
# Three interferon pathways resolved at PROTEIN level (ZS: measure the LIGANDS + Type I/II/III distinctly, not just downstream type-II-biased chemokines).
TYPE_I  =['IFNA1','IFNA2','IFNB1','SIGLEC1','ISG15','LGALS3BP','TNFSF10','BST2','LY6E','MX1']  # ligands(<LOD/off-Explore) + SIGLEC1(only type-I-SPECIFIC Olink analyte; York 2007) + secreted type-I-IFN-induced to DISCOVER on panel: ISG15, LGALS3BP(Mac-2BP/90K), TNFSF10(TRAIL; MeMed-BV), BST2(tetherin), LY6E
TYPE_II =['IFNG','IL18BP']                                  # IFNG (MAS/flare sensor, often <LOD) + IL18BP (IFN-gamma-inducible; Moller 2003)
TYPE_III=['IFNL1','IFNL2','IFNL3','IL29']                   # type-III ligands (lambda / IL28-29) -- epithelial; expect near-LOD in plasma
CHEMO_II=['CXCL9']                                          # IFN-gamma-SPECIFIC chemokine (type-II readout; Bracaglia 2017)
CHEMO_SHARED=['CXCL10','CXCL11','LGALS9','B2M']             # dual type-I/II-inducible -- NOT type-I-specific (the run-055/062 confound)
INFLAM_PROT=['IL6','TNF','IL18','IL1B','CXCL8','CCL2']      # NF-kB / inflammatory
LIGANDS=['IFNA1','IFNA2','IFNB1','IFNG','IFNL1','IFNL2','IFNL3','IL29']   # the actual IFN ligands (discovery/coverage)
TARGETS=list(dict.fromkeys(TYPE_I+TYPE_II+TYPE_III+CHEMO_II+CHEMO_SHARED+INFLAM_PROT))
def sh(c): return subprocess.run(['bash','-lc',c],capture_output=True,text=True).stdout
def fnum(x):
    try: return float(x)
    except (TypeError,ValueError): return 0.0
S={}
try:
    hdr=sh(f'zcat "{OLINK}" 2>/dev/null | head -1').rstrip("\n"); cols=hdr.split("\t")
    print("OLINK path:",OLINK); print("Olink header cols:",cols[:20])
    low=[c.lower() for c in cols]
    def findcol(cands):
        for i,c in enumerate(low):
            if any(k==c for k in cands): return i,cols[i]
        for i,c in enumerate(low):
            if any(k in c for k in cands): return i,cols[i]
        return None,None
    si,sc=findcol(['research_id','person_id','sampleid','sample_id','pid','sample'])
    ai,ac=findcol(['assay','gene_symbol','genesymbol','gene','protein','target'])
    vi,vc=findcol(['npx','value','normalized','level'])
    print(f"detected: sample={sc}(col#{si}) assay={ac}(col#{ai}) value={vc}(col#{vi})")
    if si is None or ai is None or vi is None:
        print("run failed: could not auto-detect Olink columns (see header above)"); print("run complete")
    else:
        pat="|".join(TARGETS)
        # substring /IFN/ discovers every interferon-named assay actually on the panel (ligands + receptors) in one pass; exact-match the chemokine/inflammatory set
        raw=sh(f'zcat "{OLINK}" 2>/dev/null | awk -F"\\t" \'NR==1 || ${ai+1} ~ /IFN/ || ${ai+1} ~ /^({pat})$/\'')
        m=pd.read_csv(io.StringIO(raw),sep="\t",dtype=str)
        m['__npx']=pd.to_numeric(m[vc],errors='coerce'); m['__s']=m[sc].astype(str); m['__a']=m[ac].astype(str)
        ifn_on_panel=sorted(a for a in m.__a.unique() if 'IFN' in a.upper())
        print("INTERFERON assays present on this Olink panel:",ifn_on_panel)   # <- tells us whether IFN-alpha/beta/lambda are even measurable in plasma
        m=m[m.__a.isin(TARGETS)].dropna(subset=['__npx'])
        print("target proteins used:",sorted(m.__a.unique())); S['ifn_on_panel']=ifn_on_panel
        w=m.pivot_table(index='__s',columns='__a',values='__npx',aggfunc='mean'); w.index=w.index.astype(str)
        # ---- BRIDGE Olink SampleID -> person_id via the proteomics manifest (SampleID is a plate barcode, not person_id) ----
        MANIFEST=f"{MNT}/v9/multiomics/proteomics/manifest.tsv"
        try:
            man=pd.read_csv(MANIFEST,sep="\t",dtype=str); print("manifest cols:",list(man.columns))
            pidcol=next((c for c in man.columns if c.lower() in ('person_id','research_id','participant_id','pid')),None)
            ols=set(w.index); sidcol=None; best=0
            for c in man.columns:
                ov=len(set(man[c].dropna().astype(str))&ols)
                if ov>best: best=ov; sidcol=c
            print(f"   manifest bridge: SampleID col='{sidcol}' overlap {best}/{len(ols)} -> person col='{pidcol}'")
            if sidcol and pidcol and best>0:
                bridge=dict(zip(man[sidcol].astype(str),man[pidcol].astype(str)))
                w=w.rename(index=lambda s: bridge.get(str(s))); w=w[w.index.notna()]
                w=w.groupby(level=0).mean(numeric_only=True); print(f"   mapped to {len(w)} persons")
            else: print("   ** no manifest bridge column matched Olink SampleIDs -> inspect manifest cols above **")
        except Exception as ee: print("   manifest bridge failed:",str(ee)[:120])
        def zc(df,cc):
            cc=[c for c in cc if c in df.columns]
            if not cc: return None
            z=(df[cc]-df[cc].mean())/df[cc].std(); return z.mean(axis=1)
        d=pd.DataFrame(index=w.index)
        d['typeI_IFN']=zc(w,TYPE_I); d['typeII_IFN']=zc(w,TYPE_II+CHEMO_II); d['typeIII_IFN']=zc(w,TYPE_III)
        d['shared_chemokine']=zc(w,CHEMO_SHARED); d['inflam_prot']=zc(w,INFLAM_PROT)
        pcov={k:[c for c in v if c in w.columns] for k,v in [('typeI_ligands',TYPE_I),('typeII(IFNG+CXCL9)',TYPE_II+CHEMO_II),('typeIII_ligands',TYPE_III),('shared_chemokine',CHEMO_SHARED),('inflam',INFLAM_PROT)]}
        print("pathway coverage (assays used per composite):",pcov); S['pathway_coverage']=pcov
        for c in TARGETS:
            if c in w.columns: d[c]=(w[c]-w[c].mean())/w[c].std()
        d['research_id']=d.index.astype(str); print(f"Olink persons with target proteins: {len(d)}")
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
        C_r85h=carr([R85H_VID]); C_irak3=carr(irak3); st={k:carr([v]) for k,v in STINGV.items()}; C_aq=(st['G230A']&st['R293Q'])-st['R71H']; C_haq=st['R71H']&st['G230A']&st['R293Q']
        d['R85H']=d.research_id.isin(C_r85h).astype(int); d['IRAK3_LoF']=d.research_id.isin(C_irak3).astype(int); d['AQ']=d.research_id.isin(C_aq).astype(int); d['HAQ']=d.research_id.isin(C_haq).astype(int); d['R85H_x_AQ']=(d.R85H&d.AQ).astype(int); d['R85H_x_HAQ']=(d.R85H&d.HAQ).astype(int)
        print(f"Olink carriers: R85H {int(d.R85H.sum())} | IRAK3-LoF {int(d.IRAK3_LoF.sum())} | AQ {int(d.AQ.sum())} | HAQ {int(d.HAQ.sum())} | R85H&AQ {int(d.R85H_x_AQ.sum())} | R85H&HAQ {int(d.R85H_x_HAQ.sum())}")
        if int(d.R85H.sum())==0: print("   ** 0 R85H carriers in Olink -> SampleID likely != person_id; check proteomics/manifest.tsv for the ID bridge **")
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
        try:  # smoking (EHR ever-smoker, ICD tobacco): plasma inflammatory-protein confounder -- parity with EMR + RNA models
            smk=set(str(r[0]) for r in bq.query(f"SELECT DISTINCT co.person_id FROM `{PROJ}.{DS}.condition_occurrence` co JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id WHERE c.vocabulary_id LIKE 'ICD%' AND (c.concept_code LIKE 'Z72.0%' OR c.concept_code LIKE 'F17%' OR c.concept_code LIKE '305.1%' OR c.concept_code LIKE 'V15.82%')"))
            d['smoker']=d.research_id.isin(smk).astype(float); COV+=['smoker']; print(f"   adj: ever-smoker {int(d.smoker.sum())} added (smoking parity)")
        except Exception: pass
        covs=(" + "+" + ".join(COV)) if COV else ""
        import statsmodels.formula.api as smf
        def ols(y,ex):
            sub=d.dropna(subset=[y])
            try: r=smf.ols(f'{y} ~ {ex}{covs}',data=sub,missing='drop').fit(); return round(float(r.params[ex]),3),round(float(r.pvalues[ex]),4),int(sub[ex].sum())
            except Exception: return None,None,int(sub[ex].sum())
        outc=['typeI_IFN','typeII_IFN','typeIII_IFN','shared_chemokine','inflam_prot']+[c for c in ['SIGLEC1','ISG15','LGALS3BP','TNFSF10','BST2','LY6E','IFNG','CXCL9','IL18BP','CXCL10','CXCL11','LGALS9','B2M'] if c in d.columns]
        print("\n== plasma protein ~ exposure (beta(p)) -- Type I/II/III IFN RESOLVED; 16PC + age + sex + smoking adjusted ==")
        print(f"{'exposure':12s} "+" ".join(f"{o:>13s}" for o in outc))
        S['results']={}
        for ex in ['R85H','IRAK3_LoF','AQ','HAQ','R85H_x_AQ','R85H_x_HAQ']:
            row=[]; S['results'][ex]={}
            for o in outc:
                b,p,nc=ols(o,ex); S['results'][ex][o]={'beta':b,'p':p,'n_carr':nc}; row.append(f"{b}({p})")
            print(f"{ex:12s} "+" ".join(f"{x:>13s}" for x in row))
        print("\n== READ: is TYPE-I (IFN-a/b ligands) even measurable in plasma (see coverage)? does R85H×AQ elevate TYPE-I specifically, vs type-II(IFNG+CXCL9)/type-III(IFNL)? Prior CXCL9/10/11 composite was type-II-weighted -> the run-055/062 R85H×AQ null may be a pathway mismatch. ==")
        print("\n===== PROTEOMICS IFN (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1400:])
