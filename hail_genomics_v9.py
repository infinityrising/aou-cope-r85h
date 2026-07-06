#!/usr/bin/env python
"""AoU v9 COPE-R85H / IRAK3 — Hail genomics (Phase 1, DESCRIPTIVE / L0).
REQUIRES the 'Hail Genomic Analysis' app (Dataproc + Spark) — NOT standard JupyterLab.
  Spin up: Environments -> delete current -> create -> App = 'JupyterLab Spark cluster' (Hail).
  Preflight (terminal):  python -c "import hail, pyspark; print(hail.__version__, pyspark.__version__)"
Run:  %run aou-cope-r85h/hail_genomics_v9.py
Does: (A) validate Hail + discover gs:// genomic paths, (B) R85H zygosity + per-ancestry + HWE,
      (C) STING unphased diplotypes, (D) innate-pathway burden mask [scaffold].
DESCRIPTIVE ONLY — NO confirmatory association (that is SHA-locked + run once, see PREREGISTRATION.md).
NOTE: Hail on Dataproc reads gs:// (workers can't see the FUSE mount); pandas/driver reads use gs:// via import_table.
Ends 'run complete' / 'run failed'.
"""
import subprocess, json
BUCKET = "gs://vwb-aou-datasets-controlled"
GS  = f"{BUCKET}/v9/wgs/short_read/snpindel"
ANC_GS = f"{GS}/aux/ancestry/ancestry_preds.tsv"
R85H = ('chr19', 18911007, 'C', 'T')
STING = {'R71H': ('chr5', 139481493, 'C', 'T'), 'G230A': ('chr5', 139478340, 'C', 'G'),
         'R293Q': ('chr5', 139477397, 'C', 'T'), 'R220H': ('chr5', 139478370, 'C', 'T')}
INNATE = ['IRAK1', 'IRAK3', 'IRAK4', 'MYD88', 'TIRAP', 'TLR2', 'TLR4', 'TLR7', 'TLR9', 'NFKB1', 'TNFAIP3']
S = {}

def gsls(p): return subprocess.run(['bash', '-lc', f'gsutil ls "{p}" 2>/dev/null | head -20'],
                                   capture_output=True, text=True).stdout.strip().splitlines()

try:
    import hail as hl
    hl.init(default_reference='GRCh38', log='/tmp/hail.log', quiet=True)
    print(f"== Hail {hl.version()} up ==")

    # --- A. discover gs:// genomic paths (first-run validation) ---
    print("== A. gs:// path discovery ==")
    for lab, p in [('exome/splitMT', f"{GS}/exome/splitMT/"), ('acaf/splitMT', f"{GS}/acaf_threshold/splitMT/"),
                   ('vds', f"{GS}/vds/"), ('aux/vat', f"{GS}/aux/vat/")]:
        print(f"  {lab}: {gsls(p)[:3]}")
    EXOME_MT = next((l.rstrip('/') for l in gsls(f"{GS}/exome/splitMT/") if l.rstrip('/').endswith('.mt')), None)
    print("  -> EXOME_MT =", EXOME_MT)
    if not EXOME_MT:
        print("run failed: exome splitMT .mt not found under", f"{GS}/exome/splitMT/ — inspect the listing above"); raise SystemExit

    # ancestry keyed by sample id (string)
    anc = hl.import_table(ANC_GS, types={'research_id': 'str'}, min_partitions=8).key_by('research_id')

    # --- B. R85H zygosity + per-ancestry + HWE ---
    print("== B. R85H zygosity + per-ancestry + HWE ==")
    mt = hl.read_matrix_table(EXOME_MT)
    iv = hl.locus_interval(R85H[0], R85H[1], R85H[1] + 1, reference_genome='GRCh38')
    mt = hl.filter_intervals(mt, [iv])
    mt = mt.filter_rows((mt.alleles[0] == R85H[2]) & (mt.alleles[1] == R85H[3]))
    mt = mt.annotate_cols(anc=anc[mt.s].ancestry_pred)
    zyg = mt.aggregate_entries(hl.agg.filter(hl.is_defined(mt.GT),
              hl.agg.group_by(mt.anc, hl.agg.counter(mt.GT.n_alt_alleles()))))
    zyg = {a: {int(k): int(v) for k, v in c.items()} for a, c in zyg.items() if a is not None}
    S['R85H_zyg_by_anc'] = zyg
    afr = zyg.get('afr', {})
    print(f"  per-ancestry {{0:homref,1:het,2:homalt}}: {zyg}")
    print(f"  ★ AFR: het={afr.get(1,0)} homalt={afr.get(2,0)}  (homozygotes = the 'healthy homozygotes exist' number)")
    hwe = mt.filter_cols(mt.anc == 'afr').annotate_rows(
              hwe=hl.agg.hardy_weinberg_test(mt.GT)).rows().collect()
    if hwe:
        S['R85H_AFR_HWE'] = {'het_freq_hwe': hwe[0].hwe.het_freq_hwe, 'p': hwe[0].hwe.p_value}
        print(f"  AFR HWE: p={hwe[0].hwe.p_value:.3g} (informative — R85H common enough for HWE)")

    # --- C. STING unphased diplotypes (per-person alt-allele counts) ---
    print("== C. STING unphased genotypes (phasing = long-read, sting_haplotypes step) ==")
    ivs = [hl.locus_interval(c, p, p + 1, reference_genome='GRCh38') for (c, p, _, _) in STING.values()]
    st = hl.filter_intervals(hl.read_matrix_table(EXOME_MT), ivs) if EXOME_MT else None
    # note: STING is chr5; if not in the exome MT (non-coding positions), switch to acaf/splitMT — flagged for iteration
    print("  [ADJUST] confirm STING chr5 positions are in the exome MT; else read acaf_threshold/splitMT. Rows found:",
          st.count_rows() if st is not None else 'NA')

    # --- D. innate-pathway burden mask [SCAFFOLD — needs the VAT annotations] ---
    print("== D. innate-pathway burden mask [scaffold] ==")
    print(f"  genes: {INNATE}")
    print("  [ADJUST] read exome splitMT + join aux/vat (gnomAD-AFR AF, LOFTEE, consequence);")
    print("           MASK-A = HC pLoF (LOFTEE HC, AF<0.001); MASK-B = +REVEL>0.5 missense; per-person carrier flag.")
    print("           (Built here on Hail; feeds the epistasis arm. IRAK3 alone = 6 carriers, so pathway burden is the test.)")

    print("\n===== HAIL SUMMARY (paste back) =====")
    print(json.dumps(S, indent=1, default=str))
    print("run complete")
except SystemExit:
    pass
except Exception as e:
    import traceback
    print("run failed:", type(e).__name__, str(e)[:300])
    print(traceback.format_exc()[-800:])
