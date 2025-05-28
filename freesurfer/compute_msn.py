#!/usr/bin/env python3

mris_anatomical_stats -a /Volumes/LA_4TB/datasets/mrgfus/derivatives/freesurfer/sub-001/ses-01/label/lh.500.aparc.annot  -b sub-001/ses-01 lh

"""
compute_msn.py

Compute Morphometric Similarity Networks (MSNs) from FreeSurfer stats.

Usage examples:

  # match lh.500.aparc.stats & rh.500.aparc.stats
  
  python compute_msn.py \
    --root data/ \
    --outdir out/ \
    --stats 500.aparc

"""

import os
import re
import argparse
import pandas as pd
import numpy as np
from scipy.stats import zscore

METRICS = [
    "SurfArea", "GrayVol", "ThickAvg", "ThickStd",
    "MeanCurv", "GausCurv", "FoldInd", "CurvInd",
]

def parse_stats_file(fname):
    df = pd.read_csv(
        fname,
        delim_whitespace=True,
        comment="#",
        header=None,
        names=["Region", "NumVert"] + METRICS
    )
    return df[["Region"] + METRICS]

def compute_msn_for_session(stats_dir, file_pattern, out_csv):
    regex = re.compile(file_pattern)
    dfs = []
    for fn in sorted(os.listdir(stats_dir)):
        if regex.match(fn):
            dfs.append(parse_stats_file(os.path.join(stats_dir, fn)))

    if not dfs:
        raise RuntimeError(f"No files matching `{file_pattern}` in {stats_dir!r}")

    df = pd.concat(dfs, ignore_index=True).set_index("Region")
    mat = df.apply(lambda col: zscore(col, nan_policy="omit"), axis=0).values
    corr = np.corrcoef(mat)
    np.fill_diagonal(corr, 0.0)
    pd.DataFrame(corr).to_csv(out_csv, header=False, index=False)

def main(root, outdir, stats_name):
    # decide which pattern to use
    if stats_name:
        # escape any regex metachars in stats_name
        core = re.escape(stats_name)
        pattern = rf'^(?:lh|rh)\.{core}\.stats$'
    else:
        pattern = file_regex

    for subj in sorted(os.listdir(root)):
        subjdir = os.path.join(root, subj)
        if not os.path.isdir(subjdir):
            continue
        for ses in sorted(os.listdir(subjdir)):
            statsdir = os.path.join(subjdir, ses, "stats")
            if not os.path.isdir(statsdir):
                continue

            od = os.path.join(outdir, subj)
            os.makedirs(od, exist_ok=True)
            out_csv = os.path.join(od, f"{ses}_msn.csv")
            print(f"→ {subj}/{ses} → {out_csv}")
            compute_msn_for_session(statsdir, pattern, out_csv)

if __name__ == "__main__":
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--root",     required=True,
                   help="Root dir containing sub-*/ses-*/stats")
    p.add_argument("--outdir",   required=True,
                   help="Where to write per-session MSN CSVs")
    p.add_argument("--stats",
                   help="Core stats filename (e.g. '500.aparc' to match lh.500.aparc.stats & rh.500.aparc.stats)")
    args = p.parse_args()
    main(args.root, args.outdir, args.stats)

