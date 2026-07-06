#!/usr/bin/env python
"""AoU v9 COPE-R85H / IRAK3 — innate-pathway rare-LoF BURDEN MASK (Q2 second-hit). Dataproc/Spark app.
Run in a FRESH kernel: %run aou-cope-r85h/burden_mask_v9.py
Phase 1 (driver, pysam, STREAM-FILTERED to stay low-memory): fetch VAT for 11 innate genes -> MASK-A (HC pLoF)
  + MASK-B (+REVEL>0.5 missense), rare (gvs_afr_af<0.001), deduped by vid, per-gene counts.
Phase 2 (Hail): genotype mask vids in the exome MT -> per-person burden -> ★ R85H × innate-LoF DOUBLE-CARRIERS
  by ancestry (the Q2 keystone). DESCRIPTIVE / L0. Ends 'run complete' / 'run failed'.
Gene coordinates = Ensembl GRCh38 (padded ±5kb); gene_symbol filter corrects imprecision.
"""
import os, sys, subprocess, gzip, json
from collections import defaultdict
PROJECT = os.environ['GOOGLE_PROJECT']
GS   = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel"
EXOME_MT = f"{GS}/exome/splitMT/hail.mt"
ANC_GS   = f"{GS}/aux/ancestry/ancestry_preds.tsv"
VAT  = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/vat/vat_complete.bgz.tsv.gz")
GENES = {'IRAK1':('chrX',154005501,154025650), 'IRAK3':('chr12',66183995,66259622),
         'IRAK4':('chr12',43753938,43803307), 'MYD88':('chr3',38133552,38148024),
         'TIRAP':('chr11',126277497,126303845), 'TLR2':('chr4',153679050,153713537),
         'TLR4':('chr9',117699170,117729735), 'TLR7':('chrX',12861994,12895361),
         'TLR9':('chr3',52216080,52230651), 'NFKB1':('chr4',102495911,102622302),
         'TNFAIP3':('chr6',137861383,137890836)}
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
        for line in tbx.fetch(c, s, e):                 # generator: one row at a time, never the whole region
            f = line.split("\t")
            if f[IX['gene_symbol']] != sym: continue
            if fnum(f[IX['gvs_afr_af']]) >= 0.001: continue
            seen += 1
            vid, lof, cons, rev = f[IX['vid']], f[IX['LoF']], f[IX['consequence']], f[IX['revel']]
            is_A = (lof == 'HC')
            if is_A: A.add(vid); pergA[sym].add(vid)
            if is_A or ('missense' in cons and fnum(rev) > 0.5): B.add(vid)
        print(f"   {sym:8s} rare rows: {seen:>5} | HC-pLoF vids: {len(pergA[sym])}")
    A_vids, B_vids = sorted(A), sorted(B)
    S['maskA_n'], S['maskB_n'] = len(A_vids), len(B_vids)
    S['maskA_per_gene'] = {g: len(v) for g, v in pergA.items()}
    print(f"   MASK-A (HC pLoF) unique vids = {len(A_vids)} | per gene: {S['maskA_per_gene']}")
    print(f"   MASK-B (+REVEL>0.5 missense) unique vids = {len(B_vids)}")

    # ---- Phase 2: genotype mask vids in the exome MT -> burden -> R85H cross ----
    import hail as hl
    try:
        hl.init(default_reference='GRCh38', gcs_requester_pays_configuration=PROJECT, quiet=True)
    except Exception as e:
        print("   (hail already initialized:", str(e)[:60], ")")
    mt = hl.read_matrix_table(EXOME_MT)
    anc = hl.import_table(ANC_GS, types={'research_id': 'str'}).key_by('research_id')

    def keyed_ht(vids):
        t = hl.Table.parallelize([{'v': v} for v in vids], hl.tstruct(v=hl.tstr))
        t = t.annotate(pv=hl.parse_variant('chr' + t.v.replace('-', ':'), reference_genome='GRCh38'))
        return t.key_by(locus=t.pv.locus, alleles=t.pv.alleles).select()

    print("== Phase 2: per-person burden (MASK-A) + R85H cross ==")
    mA = mt.semi_join_rows(keyed_ht(A_vids))
    print("   MASK-A variants matched in exome MT:", mA.count_rows(), "of", len(A_vids))
    burden = mA.annotate_cols(bden=hl.agg.any(mA.GT.n_alt_alleles() > 0)).cols()
    r = hl.filter_intervals(mt, [hl.locus_interval('chr19', 18911007, 18911008, reference_genome='GRCh38')])
    r = r.filter_rows((r.alleles[0] == 'C') & (r.alleles[1] == 'T'))
    r85h = r.annotate_cols(r85h=hl.agg.any(r.GT.n_alt_alleles() > 0)).cols()
    comb = r85h.annotate(bden=hl.coalesce(burden[r85h.s].bden, False), anc=anc[r85h.s].ancestry_pred)
    ct = comb.aggregate(hl.struct(
        burden_carriers=hl.agg.count_where(comb.bden),
        burden_by_anc=hl.agg.filter(comb.bden, hl.agg.counter(comb.anc)),
        double=hl.agg.count_where(comb.r85h & comb.bden),
        double_by_anc=hl.agg.filter(comb.r85h & comb.bden, hl.agg.counter(comb.anc))))
    S['innate_burden_carriers'] = int(ct.burden_carriers)
    S['innate_burden_by_anc'] = {k: int(v) for k, v in ct.burden_by_anc.items() if k}
    S['R85H_x_innateLoF_double'] = int(ct.double)
    S['double_by_anc'] = {k: int(v) for k, v in ct.double_by_anc.items() if k}
    print(f"   innate-LoF (MASK-A) burden carriers = {S['innate_burden_carriers']:,} by anc {S['innate_burden_by_anc']}")
    print(f"   ★ R85H × innate-LoF DOUBLE-CARRIERS = {S['R85H_x_innateLoF_double']} by anc {S['double_by_anc']}")

    print("\n===== BURDEN SUMMARY (paste back) =====")
    print(json.dumps(S, indent=1, default=str))
    print("run complete")
except Exception as e:
    import traceback
    print("run failed:", type(e).__name__, str(e)[:300])
    print(traceback.format_exc()[-1200:])
