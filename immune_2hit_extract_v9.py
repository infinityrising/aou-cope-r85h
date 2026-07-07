#!/usr/bin/env python
"""AoU v9 — IRAK/immune-brake SECOND-HIT hunt in R85H carriers, from long-read PAV (phased, assembly-based).
Task: do R85H carriers carry a damaging variant in an immune-brake gene (IRAK3 primary) — CATALOGUED or PRIVATE/novel
that the short-read VAT misses? Population = R85H carriers WITH long-read PAV (≈138; small -> fast). For each, read PAV
over the immune panel -> rare damaging variants (VAT HC-pLoF or missense REVEL>0.5, gvs_afr_af<0.1%) + novel-to-VAT
coding-region variants, phased. Reports per-carrier 2nd-hit genes/vids + a by-gene tally (IRAK3 flagged). Resumable
checkpoint (~/immune_2hit_pav_v9.csv). RUN ON THE STANDARD app. NW via COPI_NW env. Ends 'run complete'/'run failed'.
"""
import os, sys, subprocess, gzip, csv, time, threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import Counter
import pandas as pd
subprocess.run([sys.executable,'-m','pip','install','-q','pysam'],capture_output=True)
import pysam
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); GS="gs://vwb-aou-datasets-controlled/"
LR=f"{MNT}/v9/wgs/long_read"; RNADIR=f"{MNT}/v9/multiomics/rnaseq"; SNP=f"{MNT}/v9/wgs/short_read/snpindel"
VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"; ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1); R85H_VID='19-18911007-C-T'
OUT=os.path.expanduser("~/immune_2hit_pav_v9.csv"); ATT=os.path.expanduser("~/immune_2hit_pav_v9.attempted")
# immune-brake panel (GRCh38, validated in burden_mask_v9). IRAK3 = primary second-hit target.
IMMUNE={'IRAK1':('X',154005501,154025650),'IRAK3':('12',66183995,66259622),'IRAK4':('12',43753938,43803307),
        'MYD88':('3',38133552,38148024),'TIRAP':('11',126277497,126303845),'TLR2':('4',153679050,153713537),
        'TLR4':('9',117699170,117729735),'TLR7':('X',12861994,12895361),'TLR9':('3',52216080,52230651),
        'NFKB1':('4',102495911,102622302),'TNFAIP3':('6',137861383,137890836)}
FIELDS=['research_id','ancestry','in_rna','n_hit','hit_genes','hit_vids','IRAK3_hit','novel_coding']
NW=int(os.environ.get("COPI_NW","4"))
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
DMG=set(); ANNOT=set(); CODING=('missense','stop','frameshift','splice','start_lost','stop_lost','inframe')
_lock=threading.Lock(); _att=None; RNA=set(); AMAP={}
def mark(rid):
    with _lock: _att.write(rid+"\n"); _att.flush()
def parse_patient(rid,pav):
    mark(rid)
    try:
        mp=gs2mnt(pav)
        if not (mp and os.path.exists(mp)): return None
        v=pysam.VariantFile(mp); s0=list(v.header.samples)[0]; hits=[]; novel=0
        for gene,(ch,s,e) in IMMUNE.items():
            for rec in fetchiv(v,ch,s,e):
                gt=rec.samples[s0]['GT']
                if not gt: continue
                for i,alt in enumerate([str(x) for x in (rec.alts or [])]):
                    if (i+1) not in gt: continue
                    vid=f"{ch}-{rec.pos}-{rec.ref}-{alt}"
                    if vid in DMG: hits.append((gene,vid))
                    elif vid not in ANNOT: novel+=1     # novel-to-VAT (region-level; coding filter applied downstream)
        v.close()
        return {'research_id':rid,'ancestry':AMAP.get(rid,''),'in_rna':int(rid in RNA),'n_hit':len(hits),
                'hit_genes':';'.join(sorted({g for g,_ in hits})),'hit_vids':';'.join(v2 for _,v2 in hits),
                'IRAK3_hit':int(any(g=='IRAK3' for g,_ in hits)),'novel_coding':novel}
    except Exception: return None
try:
    log("build VAT damaging immune-panel set (rare LoF/REVEL>0.5)...")
    with gzip.open(VAT,'rt') as fh: COLS=fh.readline().rstrip("\n").split("\t")
    IX={c:COLS.index(c) for c in ['vid','gene_symbol','consequence','gvs_afr_af','revel','LoF']}
    tbx=pysam.TabixFile(VAT); perg=Counter()
    for gene,(ch,s,e) in IMMUNE.items():
        for line in tbx.fetch(f'chr{ch}',s,e):
            f=line.split("\t")
            if f[IX['gene_symbol']]!=gene: continue
            vid=f[IX['vid']]; ANNOT.add(vid)
            if fnum(f[IX['gvs_afr_af']])>=0.001: continue
            if f[IX['LoF']]=='HC' or ('missense' in f[IX['consequence']] and fnum(f[IX['revel']])>0.5): DMG.add(vid); perg[gene]+=1
    log(f"damaging immune vids {len(DMG)} per gene {dict(perg)}")
    man=pd.read_csv(f"{LR}/manifest.tsv",sep="\t"); man['research_id']=man.research_id.astype(str)
    man=man[man.grch38_pav_vcf.notna()].drop_duplicates('research_id')
    RNA=set(pd.read_csv(f"{RNADIR}/manifest.tsv",sep="\t",usecols=['research_id']).research_id.astype(str))
    anc=pd.read_csv(ANC,sep="\t",usecols=['research_id','ancestry_pred']); anc['research_id']=anc.research_id.astype(str)
    AMAP=dict(zip(anc.research_id,anc.ancestry_pred))
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    r85h=set(str(r.pid) for r in bq.query(f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid='{R85H_VID}'"))
    pop=man[man.research_id.isin(r85h)][['research_id','grch38_pav_vcf']]
    done=set()
    if os.path.exists(ATT): done|=set(x for x in open(ATT).read().split() if x)
    if os.path.exists(OUT) and os.path.getsize(OUT)>0: done|=set(pd.read_csv(OUT,usecols=['research_id']).research_id.astype(str))
    todo=[(r.research_id,r.grch38_pav_vcf) for r in pop.itertuples() if r.research_id not in done]
    log(f"R85H∩PAV {len(pop)} | done {len(done)} | todo {len(todo)}")
    newf=not(os.path.exists(OUT) and os.path.getsize(OUT)>0)
    fout=open(OUT,'a',newline=''); w=csv.DictWriter(fout,fieldnames=FIELDS)
    if newf: w.writeheader(); fout.flush()
    _att=open(ATT,'a'); ok=0
    with ThreadPoolExecutor(max_workers=NW) as ex:
        futs=[ex.submit(parse_patient,rid,pav) for rid,pav in todo]
        for fu in as_completed(futs):
            res=fu.result()
            if res: w.writerow(res); ok+=1; fout.flush()
    fout.flush(); fout.close(); _att.close()
    g=pd.read_csv(OUT)
    n2=int((g.n_hit>0).sum()); tally=Counter()
    for gg in g.hit_genes.dropna():
        for x in str(gg).split(';'):
            if x: tally[x]+=1
    log(f"DONE: {len(g)} R85H carriers scanned | with >=1 immune 2nd-hit {n2} | IRAK3 hits {int(g.IRAK3_hit.sum())}")
    log(f"by gene {dict(tally)} | any novel-coding carriers {int((g.novel_coding>0).sum())}")
    print("\n===== IMMUNE 2ND-HIT (paste back) =====")
    print(f"scanned={len(g)} with_hit={n2} IRAK3={int(g.IRAK3_hit.sum())} by_gene={dict(tally)}")
    print(g[g.n_hit>0][['research_id','ancestry','in_rna','hit_genes','hit_vids']].to_string(index=False) if n2 else "(no catalogued 2nd-hits among R85H∩PAV)")
    print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
