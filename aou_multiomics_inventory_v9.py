#!/usr/bin/env python
"""AoU v9 -- INVENTORY of available data modalities: are PROTEOMICS / metabolomics present, and where? Lists the mount's
v9 tree, flags proteomics/olink/somascan/metabolomics dirs, and for any tabular omics prints the header + a sample/feature
count. Goal: know exactly what juice AoU v9 has (beyond the RNA-seq + EMR we've used) before designing the squeeze.
Read-only. Standard app. Ends 'run complete' / 'run failed'.
"""
import os, subprocess
MNT=os.path.expanduser("~/workspace/vwb-aou-datasets-controlled-v9")
def sh(c): return subprocess.run(['bash','-lc',c],capture_output=True,text=True).stdout
try:
    print("== v9 top-level ==")
    print(sh(f'ls -la "{MNT}/v9" 2>/dev/null'))
    print("== v9/multiomics tree (2 levels) ==")
    print(sh(f'find "{MNT}/v9/multiomics" -maxdepth 2 2>/dev/null | head -80'))
    print("== hunt: proteomics / olink / somascan / metabolomics / plasma anywhere in the mount ==")
    print(sh(f'find "{MNT}" -maxdepth 5 -type d \\( -iname "*prote*" -o -iname "*olink*" -o -iname "*soma*" -o -iname "*metabol*" -o -iname "*plasma*" -o -iname "*proteom*" \\) 2>/dev/null | head -30'))
    print(sh(f'find "{MNT}" -maxdepth 6 \\( -iname "*prote*" -o -iname "*olink*" -o -iname "*soma*" -o -iname "*metabolom*" \\) 2>/dev/null | grep -iE "\\.(tsv|csv|txt|parquet|gz)" | head -30'))
    print("== RNA-seq (known modality) inventory ==")
    print(sh(f'ls -la "{MNT}/v9/multiomics/rnaseq" 2>/dev/null; echo "-- rsem --"; ls -la "{MNT}/v9/multiomics/rnaseq/rsem" 2>/dev/null | head'))
    print("== any manifest/README describing modalities ==")
    print(sh(f'find "{MNT}/v9" -maxdepth 3 \\( -iname "*manifest*" -o -iname "*readme*" -o -iname "*data_dictionary*" \\) 2>/dev/null | head -20'))
    print("run complete")
except Exception as e:
    print("run failed:",type(e).__name__,str(e)[:200])
