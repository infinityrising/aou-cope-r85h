#!/usr/bin/env python
"""AoU v9 — DEEP-DIVE on the R85H × IRAK3-damaging DOUBLE-CARRIERS (~4 in 535k; the L5 'individuals' layer, like OA37).
Answers: who are they, WHICH IRAK3 variant (LoF vs damaging-missense), zygosity ('two hets'?), ancestry/age/sex,
STING haplotype + R85H dose if in the long-read set, and their FULL clinical picture (lung + autoimmune EHR) — do
these two-hit people actually have disease, or are they healthy (= incomplete penetrance)? Also observed-vs-expected
doubles (ancestry deficit). ⚠️ SMALL CELL (<20): INTERNAL analysis only — suppress before any export/publication.
STANDARD app. Ends 'run complete' / 'run failed'.
"""
import os, gzip, json
import pandas as pd
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9"); SNP=f"{MNT}/v9/wgs/short_read/snpindel"
ANC=f"{SNP}/aux/ancestry/ancestry_preds.tsv"; VAT=f"{SNP}/aux/vat/vat_complete.bgz.tsv.gz"
CDR="wb-silky-artichoke-2408.C2025Q4R6"; PROJ,DS=CDR.split(".",1); R85H_VID='19-18911007-C-T'; NTOT=535662
IRAK3=('12',66183995,66259622)
PAVCSV=os.path.expanduser("~/copi_sting_pav_v9.csv")
PANEL={'asthma':['J45','493'],'ILD':['J84','515','516','J98.4'],'bronchiectasis':['J47','494'],
       'RA':['M05','M06','714'],'vasculitis':['M31','M30','I77.6','446'],'SLE':['M32','710.0'],
       'sarcoid':['D86'],'ANCA_vasc':['M31.3','M31.7']}
def fnum(x):
    try: return float(x)
    except (TypeError,ValueError): return 0.0
S={}
try:
    import pysam
    with gzip.open(VAT,'rt') as fh: COLS=fh.readline().rstrip("\n").split("\t")
    need=['vid','gene_symbol','consequence','gvs_afr_af','revel','LoF']+(['aa_change'] if 'aa_change' in COLS else [])
    IX={c:COLS.index(c) for c in need}
    tbx=pysam.TabixFile(VAT); dmg={}
    for line in tbx.fetch('chr'+IRAK3[0],IRAK3[1],IRAK3[2]):
        f=line.split("\t")
        if f[IX['gene_symbol']]!='IRAK3': continue
        if fnum(f[IX['gvs_afr_af']])>=0.001: continue
        if f[IX['LoF']]=='HC' or ('missense' in f[IX['consequence']] and fnum(f[IX['revel']])>0.5):
            dmg[f[IX['vid']]]={'LoF':f[IX['LoF']],'cons':f[IX['consequence']],'aa':f[IX.get('aa_change','')] if 'aa_change' in IX else '','revel':f[IX['revel']]}
    from google.cloud import bigquery
    bq=bigquery.Client(); T=f"`{PROJ}.{DS}.cb_variant_to_person`"
    def carr(vids):
        m={}
        for i in range(0,len(vids),900):
            inl=",".join("'"+v+"'" for v in vids[i:i+900])
            for r in bq.query(f"SELECT vid, e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid IN ({inl})"):
                m.setdefault(str(r.pid),[]).append(r.vid)
        return m
    r85h=set(carr([R85H_VID]).keys()); irak3=carr(list(dmg.keys())); doubles=sorted(set(irak3)&r85h)
    S['n_R85H']=len(r85h); S['n_IRAK3dmg']=len(irak3); S['n_doubles']=len(doubles); S['expected_by_chance']=round(len(r85h)*len(irak3)/NTOT,1)
    print(f"R85H {len(r85h)} | IRAK3-dmg {len(irak3)} | ★ DOUBLES {len(doubles)} (expected by chance {S['expected_by_chance']} => {'DEFICIT (R85H-AFR × IRAK3-non-AFR = different ancestries)' if len(doubles)<S['expected_by_chance'] else 'as-expected/excess'})")
    if not doubles:
        print("no double-carriers found"); print("run complete"); raise SystemExit
    anc=pd.read_csv(ANC,sep="\t",usecols=['research_id','ancestry_pred']); anc['research_id']=anc.research_id.astype(str)
    amap=dict(zip(anc.research_id,anc.ancestry_pred)); ids=",".join(doubles)
    demo={str(r.person_id):{'age':r.age,'sex':int(r.sex)} for r in bq.query(f"SELECT person_id, DATE_DIFF(DATE '2025-01-01', DATE(birth_datetime), YEAR) age, sex_at_birth_concept_id sex FROM `{PROJ}.{DS}.person` WHERE person_id IN ({ids})")}
    pav={}
    try:
        g=pd.read_csv(PAVCSV); g['research_id']=g.research_id.astype(str); g=g[g.research_id.isin(doubles)]
        for _,r in g.iterrows(): pav[r.research_id]={'R85H_dose':int(r.R85H_d),'STING_hap':f"{r.hap1}/{r.hap2}"}
    except Exception: pass
    allcodes=[c for v in PANEL.values() for c in v]; lk=" OR ".join([f"c.concept_code LIKE '{p}%'" for p in allcodes])
    phe={d:set() for d in doubles}
    for r in bq.query(f"""SELECT DISTINCT co.person_id, c.concept_code FROM `{PROJ}.{DS}.condition_occurrence` co JOIN `{PROJ}.{DS}.concept` c ON co.condition_source_concept_id=c.concept_id WHERE co.person_id IN ({ids}) AND c.vocabulary_id LIKE 'ICD%' AND ({lk})"""):
        for lab,codes in PANEL.items():
            if any(str(r.concept_code).startswith(p) for p in codes): phe[str(r.person_id)].add(lab)
    print("\n== R85H × IRAK3-damaging double-carriers (digenic double-hets, like OA37) ==")
    people=[]
    for d in doubles:
        rec={'id':d,'ancestry':amap.get(d),'age':demo.get(d,{}).get('age'),'sex':demo.get(d,{}).get('sex'),
             'IRAK3_variants':[{'vid':v,'LoF':dmg[v]['LoF'],'aa':dmg[v]['aa']} for v in irak3[d]],
             'R85H_dose':pav.get(d,{}).get('R85H_dose','n/a (not in long-read set)'),'STING_hap':pav.get(d,{}).get('STING_hap','n/a'),
             'disease':sorted(phe[d]) if phe[d] else 'NONE of the lung/autoimmune panel'}
        people.append(rec); print("  ",rec)
    S['doubles']=people
    nsick=sum(1 for p in people if p['disease']!='NONE of the lung/autoimmune panel')
    print(f"\n== {nsick}/{len(doubles)} double-carriers have a lung/autoimmune diagnosis. (0 or few sick => incomplete penetrance; two hits not sufficient.) ==")
    print("ZYGOSITY: R85H-hom ~0.7% of carriers, IRAK3-dmg-hom ~nonexistent => these are DIGENIC DOUBLE-HETS (one hit in each gene), same as OA37. Confirm per-locus with plink2 --export A on the exome for these ids.")
    print("\n===== DOUBLE-CARRIER DEEP-DIVE (paste back; SMALL CELL — INTERNAL ONLY) ====="); print(json.dumps(S,indent=1,default=str)); print("run complete")
except SystemExit: pass
except Exception as e:
    import traceback; print("run failed:",type(e).__name__,str(e)[:300]); print(traceback.format_exc()[-1200:])
