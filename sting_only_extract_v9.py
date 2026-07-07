#!/usr/bin/env python
"""AoU v9 — FAST STING-only phased-haplotype extraction from long-read PAV (1 region fetch/sample vs 10).
The full COPI+STING extractor is FUSE-mount throughput-bound (~130 min; threads don't help — the requester-pays mount
serializes reads). The FLAGSHIP answer (R85H × HAQ -> ISG) needs only the STING cis-haplotype from PAV; R85H and
damaging-COPI status come from BigQuery. So read ONLY the STING interval -> ~10x less mount I/O -> answer in minutes.
Checkpointed/resumable (~/sting_only_pav_v9.csv), R85H carriers first, MAXN cap for a quick subset. STANDARD app.
NW=COPI_NW, MAXN=MAXN (0=all) via env. Ends 'run complete' / 'run failed'.
"""
import os, sys, subprocess, csv, time, threading
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
subprocess.run([sys.executable,'-m','pip','install','-q','pysam'],capture_output=True)
import pysam
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); GS="gs://vwb-aou-datasets-controlled/"
LR=f"{MNT}/v9/wgs/long_read"; RNADIR=f"{MNT}/v9/multiomics/rnaseq"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1); R85H_VID='19-18911007-C-T'
PHENOME=int(os.environ.get("PHENOME","0"))   # 1 = ALL PAV (phenome arm — EHR, not RNA-gated); 0 = RNA∩PAV
OUT=os.path.expanduser("~/sting_phenome_pav_v9.csv" if PHENOME else "~/sting_only_pav_v9.csv")
ATT=os.path.expanduser("~/sting_phenome_pav_v9.attempted" if PHENOME else "~/sting_only_pav_v9.attempted")
STING_IV=('5',139475528,139482935)
STING_POS={139481493:('R71H','C','T'),139478340:('G230A','C','G'),139477397:('R293Q','C','T'),139478370:('R220H','C','T')}
STNAMES=['R71H','G230A','R293Q','R220H']
FIELDS=['research_id','hap1','hap2','HAQ_d','AQ_d','R220H_d']
NW=int(os.environ.get("COPI_NW","8")); MAXN=int(os.environ.get("MAXN","0"))
def log(m): print(m,file=sys.stderr,flush=True)
def gs2mnt(p): return (MNT+"/"+p[len(GS):]) if isinstance(p,str) and p.startswith(GS) else None
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
_lock=threading.Lock(); _att=None
def mark(rid):
    with _lock: _att.write(rid+"\n"); _att.flush()
def parse_patient(rid,pav):
    mark(rid)
    try:
        mp=gs2mnt(pav)
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
        v.close(); c1,c2=classify(h1),classify(h2)
        return {'research_id':rid,'hap1':c1,'hap2':c2,'HAQ_d':int(c1=='HAQ')+int(c2=='HAQ'),
                'AQ_d':int(c1=='AQ')+int(c2=='AQ'),'R220H_d':h1['R220H']+h2['R220H']}
    except Exception: return None
try:
    man=pd.read_csv(f"{LR}/manifest.tsv",sep="\t"); man['research_id']=man.research_id.astype(str)
    man=man[man.grch38_pav_vcf.notna()].drop_duplicates('research_id')
    if not PHENOME:
        rna=set(pd.read_csv(f"{RNADIR}/manifest.tsv",sep="\t",usecols=['research_id']).research_id.astype(str))
        man=man[man.research_id.isin(rna)]
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    r85h=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid='{R85H_VID}'"))
    man['pri']=man.research_id.isin(r85h).astype(int); man=man.sort_values('pri',ascending=False)
    if MAXN>0: man=man.head(MAXN)
    done=set()
    if os.path.exists(ATT): done|=set(x for x in open(ATT).read().split() if x)
    if os.path.exists(OUT) and os.path.getsize(OUT)>0: done|=set(pd.read_csv(OUT,usecols=['research_id']).research_id.astype(str))
    todo=[(r.research_id,r.grch38_pav_vcf) for r in man.itertuples() if r.research_id not in done]
    log(f"STING-only PAV∩RNA {len(man)}{' (MAXN cap)' if MAXN else ''} | done {len(done)} | todo {len(todo)} | R85H∩ {int(man.pri.sum())} | NW {NW}")
    newf=not(os.path.exists(OUT) and os.path.getsize(OUT)>0)
    fout=open(OUT,'a',newline=''); w=csv.DictWriter(fout,fieldnames=FIELDS)
    if newf: w.writeheader(); fout.flush()
    _att=open(ATT,'a'); t0=time.time(); n=0; ok=0
    with ThreadPoolExecutor(max_workers=NW) as ex:
        futs=[ex.submit(parse_patient,rid,pav) for rid,pav in todo]
        for fu in as_completed(futs):
            res=fu.result(); n+=1
            if res: w.writerow(res); ok+=1
            if n%200==0:
                fout.flush(); el=time.time()-t0
                log(f"  {n}/{len(todo)} ok={ok} {el/n:.3f}s/samp | ~{el/n*(len(todo)-n)/60:.1f}min left")
    fout.flush(); fout.close(); _att.close()
    tot=len(pd.read_csv(OUT)) if os.path.exists(OUT) else 0
    log(f"STING EXTRACT COMPLETE: +{ok} this pass | total {tot} rows in {OUT}")
    print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
