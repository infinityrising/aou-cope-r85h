#!/usr/bin/env python
"""AoU v9 COPE-R85H / IRAK3 — Hail genomics (Phase 1, DESCRIPTIVE / L0).
REQUIRES the Dataproc/Spark app ('JupyterLab Spark cluster'). Run: %run aou-cope-r85h/hail_genomics_v9.py
The controlled bucket is REQUESTER-PAYS -> hl.init is given the billing project; Hail reads gs:// (workers
can't see the FUSE mount). Does: (A) open exome MT, (B) R85H zygosity + per-ancestry + AFR HWE,
(C) STING per-variant genotypes (cross-check), (D) VAT discovery for the innate-burden mask.
DESCRIPTIVE ONLY — no confirmatory association (see PREREGISTRATION.md). Ends 'run complete' / 'run failed'.
"""
import os, subprocess, json
PROJECT = os.environ['GOOGLE_PROJECT']                       # requester-pays billing
GS  = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel"
EXOME_MT = f"{GS}/exome/splitMT/hail.mt"
ANC_GS   = f"{GS}/aux/ancestry/ancestry_preds.tsv"
R85H = ('chr19', 18911007, 'C', 'T')
STING = {'R71H': ('chr5', 139481493, 'C', 'T'), 'G230A': ('chr5', 139478340, 'C', 'G'),
         'R293Q': ('chr5', 139477397, 'C', 'T'), 'R220H': ('chr5', 139478370, 'C', 'T')}
INNATE = ['IRAK1', 'IRAK3', 'IRAK4', 'MYD88', 'TIRAP', 'TLR2', 'TLR4', 'TLR7', 'TLR9', 'NFKB1', 'TNFAIP3']
S = {}

def one_locus(mt, chrom, pos, ref, alt):
    iv = hl.locus_interval(chrom, pos, pos + 1, reference_genome='GRCh38')
    m = hl.filter_intervals(mt, [iv])
    return m.filter_rows((m.alleles[0] == ref) & (m.alleles[1] == alt))

try:
    import hail as hl
    hl.init(gcs_requester_pays_configuration=PROJECT, quiet=True)
    hl.default_reference('GRCh38')
    print(f"== Hail {hl.version()} | requester-pays billing = {PROJECT} ==")

    # --- A. open the exome MatrixTable ---
    mt = hl.read_matrix_table(EXOME_MT)
    print("== A. exome MT ==")
    print("   col_key:", list(mt.col_key), "| n_samples:", mt.count_cols(), "| row fields:", list(mt.row)[:6])
    anc = hl.import_table(ANC_GS, types={'research_id': 'str'}).key_by('research_id')

    # --- B. R85H zygosity + per-ancestry + AFR HWE (Hail counts ALT=T=R85H directly) ---
    print("== B. R85H zygosity + HWE ==")
    r = one_locus(mt, *R85H).annotate_cols(anc=anc[mt.s].ancestry_pred)
    print("   R85H rows found:", r.count_rows())
    zyg = r.aggregate_entries(hl.agg.filter(hl.is_defined(r.GT),
              hl.agg.group_by(r.anc, hl.agg.counter(r.GT.n_alt_alleles()))))
    zyg = {a: {int(k): int(v) for k, v in c.items()} for a, c in zyg.items() if a}
    S['R85H_zyg_by_anc'] = zyg
    afr = zyg.get('afr', {})
    print(f"   per-anc {{0:homref,1:het,2:homalt}}: {zyg}")
    print(f"   ★ AFR het={afr.get(1,0)} homalt={afr.get(2,0)}  (cross-check plink2: 2539 het / 21 hom)")
    ra = r.filter_cols(r.anc == 'afr')
    hwe = ra.annotate_rows(h=hl.agg.hardy_weinberg_test(ra.GT)).rows().collect()
    if hwe:
        S['R85H_AFR_HWE'] = {'het_freq_hwe': float(hwe[0].h.het_freq_hwe), 'p': float(hwe[0].h.p_value)}
        print(f"   AFR HWE: het_freq_hwe={hwe[0].h.het_freq_hwe:.4g}  p={hwe[0].h.p_value:.3g}")

    # --- C. STING per-variant genotype counts (unphased; phasing = long-read step) ---
    print("== C. STING per-variant genotypes (cross-check BQ carrier counts) ==")
    sc = {}
    for name, (c, p, rf, al) in STING.items():
        s1 = one_locus(mt, c, p, rf, al)
        cnt = s1.aggregate_entries(hl.agg.filter(hl.is_defined(s1.GT), hl.agg.counter(s1.GT.n_alt_alleles())))
        sc[name] = {int(k): int(v) for k, v in cnt.items()}
    S['STING_geno'] = sc
    print("   {0:homref,1:het,2:homalt} per STING variant:", sc)

    # --- D. VAT discovery for the innate-pathway burden mask [next iteration] ---
    print("== D. VAT structure (for the burden mask) ==")
    vat = subprocess.run(['bash', '-lc', f'gsutil -u "{PROJECT}" ls "{GS}/aux/vat/" 2>&1 | head -8'],
                         capture_output=True, text=True).stdout
    print(vat or "  (empty)")
    print(f"  [next] join VAT (gnomAD-AFR AF, LOFTEE, consequence) to exome MT over {INNATE};")
    print("         MASK-A = HC pLoF AF<0.001; MASK-B = +REVEL>0.5 missense; per-person carrier flag.")

    print("\n===== HAIL SUMMARY (paste back) =====")
    print(json.dumps(S, indent=1, default=str))
    print("run complete")
except Exception as e:
    import traceback
    print("run failed:", type(e).__name__, str(e)[:300])
    print(traceback.format_exc()[-1000:])
