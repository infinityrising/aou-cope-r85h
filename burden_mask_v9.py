#!/usr/bin/env python
"""AoU v9 COPE-R85H / IRAK3 — innate-pathway rare-LoF BURDEN MASK (Q2 second-hit). Dataproc/Spark app.
Run in a FRESH kernel: %run aou-cope-r85h/burden_mask_v9.py
Phase 1 (driver, pysam): fetch VAT for 11 innate genes -> MASK-A (HC pLoF) + MASK-B (+REVEL missense), rare
  (gvs_afr_af<0.001), deduped by vid, per-gene counts. Phase 2 (Hail): genotype mask vids in the exome MT ->
  per-person burden -> ★ R85H × innate-LoF DOUBLE-CARRIERS by ancestry (the Q2 keystone gating number).
DESCRIPTIVE / L0. Ends 'run complete' / 'run failed'.
"""
import os, sys, subprocess, gzip, json
import pandas as pd
PROJECT = os.environ['GOOGLE_PROJECT']
GS   = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel"
EXOME_MT = f"{GS}/exome/splitMT/hail.mt"
ANC_GS   = f"{GS}/aux/ancestry/ancestry_preds.tsv"
VAT  = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/vat/vat_complete.bgz.tsv.gz")
# GRCh38 gene bodies, padded ±~20kb; gene_symbol filter corrects any imprecision (per-gene counts flag a miss)
GENES = {'IRAK3':('chr12',66150000,66245000), 'IRAK1':('chrX',153986000,154038000),
         'IRAK4':('chr12',44104000,44187000), 'MYD88':('chr3',38118000,38163000),
         'TIRAP':('chr11',126263000,126319000), 'TLR2':('chr4',153664000,153728000),
         'TLR4':('chr9',117684000,117745000), 'TLR7':('chrX',12847000,12933000),
         'TLR9':('chr3',52235000,52281000), 'NFKB1':('chr4',102481000,102638000),
         'TNFAIP3':('chr6',137847000,137904000)}
S = {}

try:
    # ---- Phase 1: VAT -> mask vids ----
    subprocess.run([sys.executable,'-m','pip','install','-q','pysam'], capture_output=True, text=True)
    import pysam
    with gzip.open(VAT,'rt') as fh: COLS = fh.readline().rstrip("\n").split("\t")
    tbx = pysam.TabixFile(VAT)
    frames = []
    print("== Phase 1: VAT mask (per-gene VAT rows) ==")
    for sym,(c,s,e) in GENES.items():
        rows = [r.split("\t") for r in tbx.fetch(c, s, e)]
        g = pd.DataFrame(rows, columns=COLS)
        g = g[g.gene_symbol == sym]
        frames.append(g); print(f"   {sym:8s} {len(g):>6} rows")
    vat = pd.concat(frames, ignore_index=True)
    for col in ['gvs_afr_af','revel']: vat[col] = pd.to_numeric(vat[col], errors='coerce')
    rare = vat.gvs_afr_af.fillna(0) < 0.001
    is_missense = vat.consequence.fillna('').str.contains('missense')
    maskA = vat[(vat.LoF == 'HC') & rare]
    maskB = vat[((vat.LoF == 'HC') | (is_missense & (vat.revel > 0.5))) & rare]
    A_vids, B_vids = sorted(set(maskA.vid)), sorted(set(maskB.vid))
    S['maskA_n'], S['maskB_n'] = len(A_vids), len(B_vids)
    S['maskA_per_gene'] = maskA.groupby('gene_symbol').vid.nunique().to_dict()
    print(f"   MASK-A (HC pLoF) unique vids = {len(A_vids)} | per gene: {S['maskA_per_gene']}")
    print(f"   MASK-B (+REVEL>0.5 missense) unique vids = {len(B_vids)}")

    # ---- Phase 2: genotype mask vids in the exome MT -> burden -> R85H cross ----
    import hail as hl
    try:
        hl.init(default_reference='GRCh38', gcs_requester_pays_configuration=PROJECT, quiet=True)
    except Exception as e:
        print("   (hail already initialized:", str(e)[:60], ")")
    mt = hl.read_matrix_table(EXOME_MT)
    anc = hl.import_table(ANC_GS, types={'research_id':'str'}).key_by('research_id')

    def vids_to_keyed_ht(vids):
        t = hl.Table.parallelize([{'v': v} for v in vids], hl.tstruct(v=hl.tstr))
        t = t.annotate(pv=hl.parse_variant('chr' + t.v.replace('-', ':'), reference_genome='GRCh38'))
        return t.key_by(locus=t.pv.locus, alleles=t.pv.alleles).select()

    print("== Phase 2: per-person burden (MASK-A) + R85H cross ==")
    htA = vids_to_keyed_ht(A_vids)
    mA = mt.semi_join_rows(htA)
    print("   MASK-A variants matched in exome MT:", mA.count_rows(), "of", len(A_vids))
    burden = mA.annotate_cols(bden=hl.agg.any(mA.GT.n_alt_alleles() > 0)).cols()
    # R85H genotype per person
    r = hl.filter_intervals(mt, [hl.locus_interval('chr19', 18911007, 18911008, reference_genome='GRCh38')])
    r = r.filter_rows((r.alleles[0] == 'C') & (r.alleles[1] == 'T'))
    r85h = r.annotate_cols(r85h=hl.agg.any(r.GT.n_alt_alleles() > 0)).cols()
    # combine per person: r85h, burden, ancestry
    comb = r85h.annotate(bden=burden[r85h.s].bden, anc=anc[r85h.s].ancestry_pred)
    comb = comb.annotate(bden=hl.coalesce(comb.bden, False))
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
