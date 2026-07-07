#!/usr/bin/env python
"""AoU v9 — PHASE 1 (RESUMABLE): extract phased COPI-subunit + STING variants from long-read PAV -> checkpoint CSV.
RUN ON THE STANDARD JupyterLab app (the Dataproc/Spark driver OOM-killed the 12-thread version; standard app has the
full node RAM). Crash-proof + resumable:
  - writes each patient row to ~/copi_sting_pav_v9.csv immediately (append),
  - write-ahead log ~/copi_sting_pav_v9.attempted marks a sample BEFORE the risky pysam read, so a re-run SKIPS
    both finished and any deterministic-crasher sample (resume, never re-hit the crash),
  - R85H carriers processed FIRST (flagship-critical), flushed stderr progress, NW=4.
To restart from scratch: delete ~/copi_sting_pav_v9.csv and ~/copi_sting_pav_v9.attempted. Ends 'run complete'/'run failed'.
"""
import os, sys, subprocess, gzip, csv, time, threading
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
subprocess.run([sys.executable,'-m','pip','install','-q','pysam'],capture_output=True)
import pysam
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); GS="gs://vwb-aou-datasets-controlled/"
LR=f"{MNT}/v9/wgs/long_read"; RNADIR=f"{MNT}/v9/multiomics/rnaseq"; SNP=f"{MNT}/v9/wgs/short_read/snpindel"
VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1); R85H_VID='19-18911007-C-T'
OUT=os.path.expanduser("~/copi_sting_pav_v9.csv"); ATT=os.path.expanduser("~/copi_sting_pav_v9.attempted")
COPI={'COPA':('1',160288580,160343566),'COPB1':('11',14443357,14500027),'COPB2':('3',139353942,139389736),
      'COPG1':('3',129249575,129278068),'COPG2':('7',130505553,130668755),'COPZ1':('12',54301202,54351848),
      'COPZ2':('17',48026142,48038030),'COPE':('19',18899511,18919407),'ARCN1':('11',118572384,118603033)}
STING_IV=('5',139475528,139482935)
STING_POS={139481493:('R71H','C','T'),139478340:('G230A','C','G'),139477397:('R293Q','C','T'),139478370:('R220H','C','T')}
STNAMES=['R71H','G230A','R293Q','R220H']
FIELDS=['research_id','hap1','hap2','HAQ_d','AQ_d','R220H_d','R85H_d','copi_dmg','dmg_genes','dmg_vids','pav_novel']
NW=4
def log(m): print(m,file=sys.stderr,flush=True)
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
ANNOT=set(); COPI_DMG=set(); _lock=threading.Lock(); _att=None
def mark(rid):
    with _lock: _att.write(rid+"\n"); _att.flush()
def parse_patient(row):
    rid=row['research_id']; mark(rid)                          # write-ahead BEFORE the risky read
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
        c1,c2=classify(h1),classify(h2); r85h_d=0; dmg=[]; novel=0
        for gene,(ch,s,e) in COPI.items():
            for rec in fetchiv(v,ch,s,e):
                gt=rec.samples[s0]['GT']
                if not gt: continue
                for i,alt in enumerate([str(x) for x in (rec.alts or [])]):
                    ai=i+1
                    if ai not in gt: continue
                    vid=f"{ch}-{rec.pos}-{rec.ref}-{alt}"; dose=sum(1 for gg in gt if gg==ai)
                    if vid==R85H_VID: r85h_d=dose; continue
                    if vid in COPI_DMG: dmg.append((gene,vid))
                    elif vid not in ANNOT: novel+=1
        v.close()
        return {'research_id':rid,'hap1':c1,'hap2':c2,'HAQ_d':int(c1=='HAQ')+int(c2=='HAQ'),'AQ_d':int(c1=='AQ')+int(c2=='AQ'),
                'R220H_d':h1['R220H']+h2['R220H'],'R85H_d':r85h_d,'copi_dmg':len(dmg),
                'dmg_genes':';'.join(sorted({g for g,_ in dmg})),'dmg_vids':';'.join(vv for _,vv in dmg),'pav_novel':novel}
    except Exception: return None
try:
    log("build VAT damaging-COPI set...")
    with gzip.open(VAT,'rt') as fh: COLS=fh.readline().rstrip("\n").split("\t")
    IX={c:COLS.index(c) for c in ['vid','gene_symbol','consequence','gvs_afr_af','revel','LoF']}
    tbx=pysam.TabixFile(VAT)
    for gene,(ch,s,e) in COPI.items():
        for line in tbx.fetch(f'chr{ch}',s,e):
            f=line.split("\t")
            if f[IX['gene_symbol']]!=gene: continue
            ANNOT.add(f[IX['vid']])
            if f[IX['LoF']]=='HC' or ('missense' in f[IX['consequence']] and fnum(f[IX['revel']])>0.5): COPI_DMG.add(f[IX['vid']])
    log(f"damaging COPI vids {len(COPI_DMG)} | annotated {len(ANNOT)}")
    man=pd.read_csv(f"{LR}/manifest.tsv",sep="\t"); man['research_id']=man.research_id.astype(str)
    man=man[man.grch38_pav_vcf.notna()].drop_duplicates('research_id')
    rna=set(pd.read_csv(f"{RNADIR}/manifest.tsv",sep="\t",usecols=['research_id']).research_id.astype(str))
    man=man[man.research_id.isin(rna)]
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    r85h=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid='{R85H_VID}'"))
    man['pri']=man.research_id.isin(r85h).astype(int); man=man.sort_values('pri',ascending=False)
    done=set()
    if os.path.exists(ATT): done|=set(x for x in open(ATT).read().split() if x)
    if os.path.exists(OUT) and os.path.getsize(OUT)>0: done|=set(pd.read_csv(OUT,usecols=['research_id']).research_id.astype(str))
    todo=[r for _,r in man.iterrows() if r['research_id'] not in done]
    log(f"PAV∩RNA {len(man)} | already attempted/done {len(done)} | todo {len(todo)} | R85H∩PAV∩RNA {int(man.pri.sum())}")
    newf=not(os.path.exists(OUT) and os.path.getsize(OUT)>0)
    fout=open(OUT,'a',newline=''); w=csv.DictWriter(fout,fieldnames=FIELDS)
    if newf: w.writeheader(); fout.flush()
    _att=open(ATT,'a'); t0=time.time(); n=0; ok=0
    with ThreadPoolExecutor(max_workers=NW) as ex:
        futs=[ex.submit(parse_patient,r) for r in todo]
        for fu in as_completed(futs):
            res=fu.result(); n+=1
            if res: w.writerow(res); ok+=1
            if n%100==0:
                fout.flush(); el=time.time()-t0
                log(f"  {n}/{len(todo)} ok={ok} {el/n:.3f}s/samp | ~{el/n*(len(todo)-n)/60:.1f}min left")
    fout.flush(); fout.close(); _att.close()
    tot=len(pd.read_csv(OUT)) if os.path.exists(OUT) else 0
    log(f"EXTRACT COMPLETE this pass: +{ok} | total rows in CSV {tot}")
    print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
