#!/usr/bin/env python
"""AoU v9 — does the SENTINEL (person_id 1348881, the COPA-phenocopy R85H carrier) carry an IRAK3 LoF like OA37 (L210*)?
Comprehensive head-on scan: ALL her CODING variants in the IRAK-brake panel (IRAK1/3/4, MYD88) — any consequence, not
just the rare-damaging set — plus her STING haplotype (AQ?) and R85H. Two informative outcomes:
  • HAS IRAK3-LoF -> CONVERGENT: two independent COPA-phenocopy patients both via R85H + IRAK3-LoF -> strengthens IRAK3.
  • LACKS it -> DIVERGENT routes: her documented 2nd hits = COPG1 I55T + STING-AQ; OA37 = IRAK3 -> supports the
    multi-route conditional-penetrance model (same phenotype, different second hits).
NB short-read only (cb_variant_to_person); an indel LoF short-read misses would need her long-read (flagged if absent).
INTERNAL (n=1). STANDARD app. Ends 'run complete' / 'run failed'.
"""
import os, gzip, json
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1); SENT='1348881'; R85H_VID='19-18911007-C-T'
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); VAT=f"{MNT}/v9/wgs/short_read/snpindel/aux/vat/vat_complete.bgz.tsv.gz"
PANEL={'IRAK3':('12',66183995,66259622),'IRAK1':('X',154005501,154025650),'IRAK4':('12',43753938,43803307),'MYD88':('3',38133552,38148024)}
STING={'R71H':'5-139481493-C-T','G230A':'5-139478340-C-G','R293Q':'5-139477397-C-T','R220H':'5-139478370-C-T'}
CODING=('missense','stop','frameshift','splice','start_lost','stop_lost','inframe')
S={'sentinel':SENT}
try:
    import pysam
    with gzip.open(VAT,'rt') as fh: COLS=fh.readline().rstrip("\n").split("\t")
    hasaa='aa_change' in COLS; IX={c:COLS.index(c) for c in ['vid','gene_symbol','consequence','gvs_afr_af','revel','LoF']+(['aa_change'] if hasaa else [])}
    tbx=pysam.TabixFile(VAT); ann={}
    for g,(ch,s,e) in PANEL.items():
        for line in tbx.fetch('chr'+ch,s,e):
            f=line.split("\t")
            if f[IX['gene_symbol']]!=g or not any(c in f[IX['consequence']] for c in CODING): continue
            ann[f[IX['vid']]]={'gene':g,'cons':f[IX['consequence']],'LoF':f[IX['LoF']],'revel':f[IX['revel']],'aa':(f[IX['aa_change']] if hasaa else ''),'af':f[IX['gvs_afr_af']]}
    print(f"immune-brake CODING vids in VAT (IRAK1/3/4+MYD88): {len(ann)}")
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    vids=list(ann.keys()); carried=[]
    for i in range(0,len(vids),900):
        inl=",".join("'"+v+"'" for v in vids[i:i+900])
        for r in bq.query(f"SELECT DISTINCT vid FROM {T}, UNNEST(person_ids.list) e WHERE vid IN ({inl}) AND e.element={SENT}"): carried.append(r.vid)
    print("== sentinel's immune-brake CODING variants ==")
    if carried:
        for v in carried: a=ann[v]; print(f"   {a['gene']:6s} {v}  {a['aa']}  LoF={a['LoF']} REVEL={a['revel']} AF={a['af']}  ({a['cons']})")
    else: print("   NONE — no coding IRAK1/3/4/MYD88 variant")
    S['immune_brake_variants']=[{'vid':v,**ann[v]} for v in carried]
    S['has_IRAK3_LoF']=any(ann[v]['gene']=='IRAK3' and ann[v]['LoF']=='HC' for v in carried)
    S['has_any_IRAK3_coding']=any(ann[v]['gene']=='IRAK3' for v in carried)
    stv=list(STING.values())+[R85H_VID]; inl=",".join("'"+v+"'" for v in stv)
    hits=set(r.vid for r in bq.query(f"SELECT DISTINCT vid FROM {T}, UNNEST(person_ids.list) e WHERE vid IN ({inl}) AND e.element={SENT}"))
    S['R85H']=R85H_VID in hits; S['STING']={k:(v in hits) for k,v in STING.items()}
    aq=S['STING']['G230A'] and S['STING']['R293Q'] and not S['STING']['R71H']
    print(f"   R85H carrier: {S['R85H']} | STING: {S['STING']}  (AQ pattern={aq})")
    print("\n== ANSWER ==")
    print(f"   Sentinel 1348881 has an IRAK3 truncating LoF (like OA37 L210*)?  {S['has_IRAK3_LoF']}")
    print(f"   Sentinel has ANY coding IRAK3 variant?  {S['has_any_IRAK3_coding']}")
    print("   => CONVERGENT (both via IRAK3-LoF)" if S['has_IRAK3_LoF'] else "   => DIVERGENT: sentinel's route is NOT IRAK3-LoF (her documented 2nd hits = COPG1 I55T + STING-AQ); OA37's is IRAK3 -> multi-route model.")
    print("   (short-read only; an indel LoF short-read misses would need her long-read PAV.)")
    print("\n===== SENTINEL IRAK3 SCAN (paste back; INTERNAL n=1) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1000:])
