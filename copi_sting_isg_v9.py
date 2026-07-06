#!/usr/bin/env python
"""AoU v9 — COPI-subunit + STING PHASED variants (long-read PAV) -> ISG panel, stratified by COPI-mutation × STING.
THE experiment: obtain haplotype-resolved variants across the 9 COPI subunits + STING from PAV (assembly-phased),
intersect with RNA-seq, then test whether a hypomorphic STING haplotype (HAQ) BLUNTS the ISG elevation of a COPI
mutation. Flagship: R85H carriers WITH vs WITHOUT HAQ -> HAQ protective (lower ISG). Broader: any damaging COPI × HAQ.
Source: manifest grch38_pav_vcf (phased, PS=None global assembly phase). Damaging = VAT HC-pLoF or missense REVEL>0.5.
STANDARD app (pysam + BigQuery + statsmodels). Threaded PAV reads. Ends 'run complete' / 'run failed'.
"""
import os, sys, subprocess, io, json, time, gzip
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import Counter
import numpy as np, pandas as pd
subprocess.run([sys.executable,'-m','pip','install','-q','pysam'],capture_output=True)
import pysam
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); GS="gs://vwb-aou-datasets-controlled/"
LR=f"{MNT}/v9/wgs/long_read"; RNADIR=f"{MNT}/v9/multiomics/rnaseq"; SNP=f"{MNT}/v9/wgs/short_read/snpindel"
TPM=f"{RNADIR}/rsem/aou_rnaseq_20260415.rsem_genes_tpm.txt.gz"; ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"
VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1); R85H_VID='19-18911007-C-T'
# GRCh38 gene intervals (Ensembl, verified). chrom as plain number; 'chr' added per source.
COPI={'COPA':('1',160288580,160343566),'COPB1':('11',14443357,14500027),'COPB2':('3',139353942,139389736),
      'COPG1':('3',129249575,129278068),'COPG2':('7',130505553,130668755),'COPZ1':('12',54301202,54351848),
      'COPZ2':('17',48026142,48038030),'COPE':('19',18899511,18919407),'ARCN1':('11',118572384,118603033)}
STING_IV=('5',139475528,139482935)
STING_POS={139481493:('R71H','C','T'),139478340:('G230A','C','G'),139477397:('R293Q','C','T'),139478370:('R220H','C','T')}
STNAMES=['R71H','G230A','R293Q','R220H']
IFN={'IFI27':'ENSG00000165949','IFI44L':'ENSG00000137959','IFIT1':'ENSG00000185745','ISG15':'ENSG00000187608','RSAD2':'ENSG00000134321','SIGLEC1':'ENSG00000088827'}
NW=12; S={}
def sh(c): return subprocess.run(['bash','-lc',c],capture_output=True,text=True).stdout
def gs2mnt(p): return (MNT+"/"+p[len(GS):]) if isinstance(p,str) and p.startswith(GS) else None
def fnum(x):
    try: return float(x)
    except (TypeError,ValueError): return 0.0
def fetchiv(v,ch,s,e):
    for c in (f'chr{ch}',ch):
        try: return list(v.fetch(c,s,e))
        except (ValueError,KeyError): continue
    return []
def classify(h):
    R,G,Q=h['R71H'],h['G230A'],h['R293Q']
    if R and G and Q: return 'HAQ'
    if G and Q and not R: return 'AQ'
    if R and not G and not Q: return 'H'
    if not R and not G and not Q: return 'R232'
    return 'partial'
# ---- built once, read-only in threads ----
ANNOT={}; COPI_DMG=set()
def parse_patient(row):
    try:
        mp=gs2mnt(row['grch38_pav_vcf'])
        if not (mp and os.path.exists(mp)): return None
        v=pysam.VariantFile(mp); s0=list(v.header.samples)[0]
        h1={n:0 for n in STNAMES}; h2={n:0 for n in STNAMES}
        for rec in fetchiv(v,*STING_IV):
            k=STING_POS.get(rec.pos)
            if not k: continue
            nm,ref,alt=k; alts=[str(x) for x in (rec.alts or [])]
            if alt not in alts: continue
            ai=alts.index(alt)+1; gt=rec.samples[s0]['GT']
            if gt and len(gt)==2:
                if gt[0]==ai: h1[nm]=1
                if gt[1]==ai: h2[nm]=1
        c1,c2=classify(h1),classify(h2)
        r85h_d=0; dmg=[]; novel=0
        for gene,(ch,s,e) in COPI.items():
            for rec in fetchiv(v,ch,s,e):
                gt=rec.samples[s0]['GT']
                if not gt: continue
                for i,alt in enumerate([str(x) for x in (rec.alts or [])]):
                    ai=i+1
                    if ai not in gt: continue
                    vid=f"{ch}-{rec.pos}-{rec.ref}-{alt}"; dose=sum(1 for g in gt if g==ai)
                    if vid==R85H_VID: r85h_d=dose; continue
                    if vid in COPI_DMG: dmg.append((gene,vid))
                    elif vid not in ANNOT: novel+=1
        v.close()
        return {'research_id':row['research_id'],'hap1':c1,'hap2':c2,
                'HAQ_d':int(c1=='HAQ')+int(c2=='HAQ'),'AQ_d':int(c1=='AQ')+int(c2=='AQ'),
                'R220H_d':h1['R220H']+h2['R220H'],'R85H_d':r85h_d,'copi_dmg':len(dmg),
                'dmg_genes':';'.join(sorted({g for g,_ in dmg})),'dmg_vids':';'.join(v for _,v in dmg),'pav_novel':novel}
    except Exception: return None
try:
    print("== 0. VAT annotation for COPI+STING intervals (damaging set) ==")
    with gzip.open(VAT,'rt') as fh: COLS=fh.readline().rstrip("\n").split("\t")
    IX={c:COLS.index(c) for c in ['vid','gene_symbol','consequence','gvs_afr_af','revel','LoF']}
    tbx=pysam.TabixFile(VAT); pergene=Counter()
    for gene,(ch,s,e) in {**COPI,'STING1':STING_IV}.items():
        for line in tbx.fetch(f'chr{ch}',s,e):
            f=line.split("\t")
            if f[IX['gene_symbol']]!=gene: continue
            vid=f[IX['vid']]; cons=f[IX['consequence']]; lof=f[IX['LoF']]; rev=fnum(f[IX['revel']])
            ANNOT[vid]=(gene,cons,lof,rev,fnum(f[IX['gvs_afr_af']]))
            if lof=='HC' or ('missense' in cons and rev>0.5): COPI_DMG.add(vid); pergene[gene]+=1
    S['copi_damaging_per_gene']=dict(pergene); S['R85H_vat']=ANNOT.get(R85H_VID)
    print(f"   damaging COPI vids {len(COPI_DMG)} per gene {dict(pergene)} | R85H VAT {ANNOT.get(R85H_VID)}")

    man=pd.read_csv(f"{LR}/manifest.tsv",sep="\t"); man['research_id']=man.research_id.astype(str)
    man=man[man.grch38_pav_vcf.notna()].drop_duplicates('research_id')
    rna_ids=set(pd.read_csv(f"{RNADIR}/manifest.tsv",sep="\t",usecols=['research_id']).research_id.astype(str))
    work=man[man.research_id.isin(rna_ids)]
    print(f"== 1. phase PAV∩RNA: {len(work)} (threads={NW}) ==")
    rows=[]; t0=time.time(); done=0
    with ThreadPoolExecutor(max_workers=NW) as ex:
        futs=[ex.submit(parse_patient,r) for _,r in work.iterrows()]
        for fu in as_completed(futs):
            r=fu.result(); done+=1
            if r: rows.append(r)
            if done==200 or done%2000==0: print(f"   {done}/{len(work)} ok={len(rows)} {(time.time()-t0)/done:.3f}s/samp -> ~{(time.time()-t0)/done*len(work)/60:.1f}min")
    g=pd.DataFrame(rows); print(f"   parsed {len(g)} in {time.time()-t0:.0f}s")
    g.to_csv(os.path.expanduser("~/copi_sting_pav_v9.csv"),index=False)
    S['n_parsed']=len(g); S['hap_freq']=dict(Counter(list(g.hap1)+list(g.hap2)))
    S['R85H_from_PAV']=int((g.R85H_d>0).sum()); S['copi_dmg_carriers']=int((g.copi_dmg>0).sum())
    print(f"   hap freq {S['hap_freq']} | R85H(PAV) {S['R85H_from_PAV']} | COPI-dmg carriers {S['copi_dmg_carriers']}")

    print("== 2. ISG score + ancestry + R85H(BQ cross-check) ==")
    raw=sh(f'zcat "{TPM}" 2>/dev/null | awk -F"\\t" \'NR==1 || $1 ~ /{"|".join(IFN.values())}/\'')
    m=pd.read_csv(io.StringIO(raw),sep="\t"); m['ensg']=m.gene_id.str.split('.').str[0]
    m=m[m.ensg.isin(IFN.values())].set_index('ensg').drop(columns=['gene_id','transcript_id(s)'])
    ex2=m.T.astype(float); ex2.index=ex2.index.astype(str); lg=np.log2(ex2+1)
    isg=((lg-lg.mean())/lg.std()).mean(axis=1); isgdf=pd.DataFrame({'research_id':isg.index,'ifn':isg.values})
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    r85h_bq=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid='{R85H_VID}'"))
    anc=pd.read_csv(ANC,sep="\t",usecols=['research_id','ancestry_pred']); anc['research_id']=anc.research_id.astype(str)
    d=g.merge(isgdf,on='research_id',how='inner').merge(anc,on='research_id',how='left')
    d['R85H_bq']=d.research_id.isin(r85h_bq).astype(int)
    S['R85H_PAV_vs_BQ']={'pav':int((d.R85H_d>0).sum()),'bq':int(d.R85H_bq.sum()),'concord':int(((d.R85H_d>0)==(d.R85H_bq>0)).sum())}
    print(f"   ISG-joined {len(d)} | R85H PAV vs BQ {S['R85H_PAV_vs_BQ']}")
    d['R85H']=((d.R85H_d>0)|(d.R85H_bq>0)).astype(int); d['HAQ']=(d.HAQ_d>0).astype(int)
    d['COPImut']=((d.copi_dmg>0)|(d.R85H>0)).astype(int)
    a=d[d.ancestry_pred=='afr'].copy()
    S['n_afr']=len(a); S['afr_R85H']=int(a.R85H.sum()); S['afr_HAQ']=int(a.HAQ.sum())

    import statsmodels.formula.api as smf
    def ols(f,dd,term):
        try: r=smf.ols(f,data=dd.dropna(subset=['ifn'])).fit(); return {'beta':round(float(r.params[term]),4),'p':round(float(r.pvalues[term]),4),'n':int(r.nobs)}
        except Exception as e: return {'error':str(e)[:70]}
    print("== 3. LADDER (validation): STING -> ISG. AQ/HAQ hypomorph => NEGATIVE ==")
    for t in ['HAQ_d','AQ_d','R220H_d']:
        S[f'ladder_all_{t}']=ols(f'ifn ~ {t} + C(ancestry_pred)',d,t); S[f'ladder_afr_{t}']=ols(f'ifn ~ {t}',a,t)
        print(f"   ifn~{t}: ALL(+anc) {S[f'ladder_all_{t}']} | AFR {S[f'ladder_afr_{t}']}")
    print("== 4. FLAGSHIP: R85H × HAQ -> ISG (HAQ protective => R85H+HAQ+ lowest; interaction NEGATIVE) ==")
    def cells(dd):
        out={}
        for r in (0,1):
            for h in (0,1):
                sub=dd[(dd.R85H==r)&(dd.HAQ==h)].ifn.dropna()
                out[f'R85H{r}_HAQ{h}']={'mean':round(float(sub.mean()),3) if len(sub) else None,'n':int(len(sub))}
        return out
    S['cells_afr']=cells(a); S['cells_all']=cells(d)
    print("   AFR 2x2:",S['cells_afr']); print("   ALL 2x2:",S['cells_all'])
    xa=a[a.R85H==1]
    S['haq_protection_afr']={'R85H+HAQ+_mean':round(float(xa[xa.HAQ==1].ifn.mean()),3) if (xa.HAQ==1).any() else None,'n+':int((xa.HAQ==1).sum()),
                             'R85H+HAQ-_mean':round(float(xa[xa.HAQ==0].ifn.mean()),3) if (xa.HAQ==0).any() else None,'n-':int((xa.HAQ==0).sum())}
    S['inter_afr']=ols('ifn ~ R85H*HAQ',a,'R85H:HAQ')
    S['inter_all_adj']=ols('ifn ~ R85H*HAQ + C(ancestry_pred)',d,'R85H:HAQ')
    print("   HAQ-protection within R85H+ (AFR):",S['haq_protection_afr'])
    print("   interaction R85H:HAQ  AFR",S['inter_afr']," | ALL(+anc)",S['inter_all_adj'])
    print("== 5. BROADER: any damaging COPI × HAQ -> ISG ==")
    S['inter_COPImut_all_adj']=ols('ifn ~ COPImut*HAQ + C(ancestry_pred)',d,'COPImut:HAQ')
    xc=d[d.COPImut==1]
    S['COPImut_HAQ_means']={'mut+HAQ+':round(float(xc[xc.HAQ==1].ifn.mean()),3) if (xc.HAQ==1).any() else None,'n+':int((xc.HAQ==1).sum()),
                            'mut+HAQ-':round(float(xc[xc.HAQ==0].ifn.mean()),3) if (xc.HAQ==0).any() else None,'n-':int((xc.HAQ==0).sum())}
    print("   COPImut×HAQ (ALL+anc):",S['inter_COPImut_all_adj']," means",S['COPImut_HAQ_means'])
    print("\n===== COPI×STING×ISG (paste back) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1400:])
