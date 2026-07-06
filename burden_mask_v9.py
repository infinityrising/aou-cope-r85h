#!/usr/bin/env python
"""AoU v9 COPE-R85H / IRAK3 — innate-pathway rare-LoF BURDEN MASK (Q2 second-hit).
RUN ON THE STANDARD JupyterLab app (BigQuery + pysam + pandas) — NO Hail/Spark. cb_variant_to_person gives
carriers directly, so we never read exome-MT genotypes (that OOM'd the Spark cluster). Run: %run aou-cope-r85h/burden_mask_v9.py
Phase 1 (pysam, streaming): VAT -> MASK-A (HC pLoF) + MASK-B (+REVEL>0.5 missense), rare (gvs_afr_af<0.001), per gene.
Phase 2 (BigQuery): union carriers of the mask vids via cb_variant_to_person -> ★ R85H × innate-LoF double-carriers by ancestry.
DESCRIPTIVE / L0. Ends 'run complete' / 'run failed'.
"""
import os, sys, subprocess, gzip, json
from collections import defaultdict
import pandas as pd
CDR = "wb-silky-artichoke-2408.C2025Q4R6"; PROJ, DS = CDR.split(".", 1)
VAT = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/vat/vat_complete.bgz.tsv.gz")
ANC = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv")
GENES = {'IRAK1':('chrX',154005501,154025650), 'IRAK3':('chr12',66183995,66259622),
         'IRAK4':('chr12',43753938,43803307), 'MYD88':('chr3',38133552,38148024),
         'TIRAP':('chr11',126277497,126303845), 'TLR2':('chr4',153679050,153713537),
         'TLR4':('chr9',117699170,117729735), 'TLR7':('chrX',12861994,12895361),
         'TLR9':('chr3',52216080,52230651), 'NFKB1':('chr4',102495911,102622302),
         'TNFAIP3':('chr6',137861383,137890836)}
R85H_VID = '19-18911007-C-T'
S = {}
def fnum(x):
    try: return float(x)
    except (TypeError, ValueError): return 0.0

try:
    # ---- Phase 1: VAT -> mask vids (streaming, low-memory) ----
    subprocess.run([sys.executable, '-m', 'pip', 'install', '-q', 'pysam'], capture_output=True, text=True)
    import pysam
    with gzip.open(VAT, 'rt') as fh: COLS = fh.readline().rstrip("\n").split("\t")
    IX = {c: COLS.index(c) for c in ['vid', 'gene_symbol', 'consequence', 'gvs_afr_af', 'revel', 'LoF']}
    tbx = pysam.TabixFile(VAT)
    A, B, pergA = set(), set(), defaultdict(set)
    print("== Phase 1: VAT mask (streaming; rare gvs_afr_af<0.1%) ==")
    for sym, (c, s, e) in GENES.items():
        seen = 0
        for line in tbx.fetch(c, s, e):
            f = line.split("\t")
            if f[IX['gene_symbol']] != sym: continue
            if fnum(f[IX['gvs_afr_af']]) >= 0.001: continue
            seen += 1
            vid, lof, cons, rev = f[IX['vid']], f[IX['LoF']], f[IX['consequence']], f[IX['revel']]
            isA = (lof == 'HC')
            if isA: A.add(vid); pergA[sym].add(vid)
            if isA or ('missense' in cons and fnum(rev) > 0.5): B.add(vid)
        print(f"   {sym:8s} rare rows {seen:>6} | HC-pLoF vids {len(pergA[sym])}")
    A_vids, B_vids = sorted(A), sorted(B)
    S['maskA_n'], S['maskB_n'] = len(A_vids), len(B_vids)
    S['maskA_per_gene'] = {g: len(v) for g, v in pergA.items()}
    print(f"   MASK-A (HC pLoF) unique vids = {len(A_vids)} | MASK-B = {len(B_vids)}")

    # ---- Phase 2: carriers via BigQuery cb_variant_to_person (cheap, clustered on vid; no Hail) ----
    from google.cloud import bigquery
    bq = bigquery.Client()
    T = f"`{PROJ}.{DS}.cb_variant_to_person`"
    def carriers_of(vids):                     # union of carriers over vids (batch <1000 keeps clustering pruning)
        ids = set()
        for i in range(0, len(vids), 900):
            inl = ",".join("'" + v + "'" for v in vids[i:i+900])
            ids |= set(int(r.pid) for r in bq.query(
                f"SELECT DISTINCT e.element pid FROM {T}, UNNEST(person_ids.list) e WHERE vid IN ({inl})"))
        return ids
    print("== Phase 2: carriers via cb_variant_to_person (BigQuery) ==")
    bcar = carriers_of(A_vids)
    rcar = carriers_of([R85H_VID])
    anc = pd.read_csv(ANC, sep="\t", usecols=['research_id', 'ancestry_pred'])
    anc['research_id'] = anc.research_id.astype('int64')
    amap = anc.set_index('research_id').ancestry_pred
    def by_anc(ids): return amap.reindex(list(ids)).value_counts().to_dict()
    dbl = bcar & rcar
    S['innate_burden_carriers'] = len(bcar); S['innate_burden_by_anc'] = by_anc(bcar)
    S['R85H_carriers'] = len(rcar)
    S['R85H_x_innateLoF_double'] = len(dbl); S['double_by_anc'] = by_anc(dbl)
    print(f"   innate-LoF (MASK-A) burden carriers = {len(bcar):,} by anc {S['innate_burden_by_anc']}")
    print(f"   R85H carriers = {len(rcar):,}")
    print(f"   ★ R85H × innate-LoF DOUBLE-CARRIERS = {len(dbl)} by anc {S['double_by_anc']}")

    print("\n===== BURDEN SUMMARY (paste back) =====")
    print(json.dumps(S, indent=1, default=str))
    print("run complete")
except Exception as e:
    import traceback
    print("run failed:", type(e).__name__, str(e)[:300])
    print(traceback.format_exc()[-1000:])
