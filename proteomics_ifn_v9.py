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
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9")
OLINK=f"{MNT}/v9/multiomics/proteomics/normalized/Olink_10k_aou_v9_normalized.tsv.gz"
SNP=f"{MNT}/v9/wgs/short_read/snpindel"; ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1)
R85H_VID='19-18911007-C-T'; STINGV={'R71H':'5-139481493-C-T','G230A':'5-139478340-C-G','R293Q':'5-139477397-C-T'}
IFN_PROT=['CXCL9','CXCL10','CXCL11']                       # IFN-inducible chemokines = plasma ISG readout
INFLAM_PROT=['IL6','TNF','IL18','IL1B','CXCL8','CCL2']
TARGETS=IFN_PROT+INFLAM_PROT
def sh(c): return subprocess.run(['bash','-lc',c],capture_output=True,text=True).stdout
def fnum(x):
    try: return float(x)
    except (TypeError,ValueError): return 0.0
S={}
try:
    hdr=sh(f'zcat "{OLINK}" 2>/dev/null | head -1').rstrip("\n"); cols=hdr.split("\t")
    print("Olink header cols:",cols[:20])
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
        raw=sh(f'zcat "{OLINK}" 2>/dev/null | awk -F"\\t" \'NR==1 || ${ai+1} ~ /^({pat})$/\'')
        m=pd.read_csv(io.StringIO(raw),sep="\t",dtype=str)
        m['__npx']=pd.to_numeric(m[vc],errors='coerce'); m['__s']=m[sc].astype(str); m['__a']=m[ac].astype(str)
        m=m[m.__a.isin(TARGETS)].dropna(subset=['__npx'])
        print("target proteins found:",sorted(m.__a.unique()))
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
        d=pd.DataFrame(index=w.index); d['IFN_prot']=zc(w,IFN_PROT); d['inflam_prot']=zc(w,INFLAM_PROT)
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
        C_r85h=carr([R85H_VID]); C_irak3=carr(irak3); st={k:carr([v]) for k,v in STINGV.items()}; C_aq=(st['G230A']&st['R293Q'])-st['R71H']
        d['R85H']=d.research_id.isin(C_r85h).astype(int); d['IRAK3_LoF']=d.research_id.isin(C_irak3).astype(int); d['AQ']=d.research_id.isin(C_aq).astype(int); d['R85H_x_AQ']=(d.R85H&d.AQ).astype(int)
        print(f"Olink carriers: R85H {int(d.R85H.sum())} | IRAK3-LoF {int(d.IRAK3_LoF.sum())} | AQ {int(d.AQ.sum())} | R85H&AQ {int(d.R85H_x_AQ.sum())}")
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
        covs=(" + "+" + ".join(COV)) if COV else ""
        import statsmodels.formula.api as smf
        def ols(y,ex):
            sub=d.dropna(subset=[y])
            try: r=smf.ols(f'{y} ~ {ex}{covs}',data=sub,missing='drop').fit(); return round(float(r.params[ex]),3),round(float(r.pvalues[ex]),4),int(sub[ex].sum())
            except Exception: return None,None,int(sub[ex].sum())
        outc=['IFN_prot','inflam_prot']+[c for c in TARGETS if c in d.columns]
        print("\n== plasma protein ~ exposure (beta(p)) -- 16PC + age + sex adjusted ==")
        print(f"{'exposure':12s} "+" ".join(f"{o:>13s}" for o in outc))
        S['results']={}
        for ex in ['R85H','IRAK3_LoF','AQ','R85H_x_AQ']:
            row=[]; S['results'][ex]={}
            for o in outc:
                b,p,nc=ols(o,ex); S['results'][ex][o]={'beta':b,'p':p,'n_carr':nc}; row.append(f"{b}({p})")
            print(f"{ex:12s} "+" ".join(f"{x:>13s}" for x in row))
        print("\n== READ: does R85H×AQ elevate IFN_prot (CXCL9/10/11) at PROTEIN level, mirroring the RNA type-I ISG signal? ==")
        print("\n===== PROTEOMICS IFN (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1400:])
