 aparcstats2table --hemi lh --subjects 004 008 --parc 500.aparc.annot --meas meancurv --tablefile lh.a2005s.meancurv.txt

# To use:
#python compute_msn.py \
#  --root /path/to/data \
#  --outdir /path/to/output


#!/usr/bin/env python3
"""
compute_msn.py

Compute Morphometric Similarity Networks (MSNs) from FreeSurfer stats.

For each subject/session:
 - reads lh.*.stats and rh.*.stats in ses-*/stats
 - extracts SurfArea, GrayVol, ThickAvg, ThickStd, MeanCurv, GausCurv, FoldInd, CurvInd
 - z-scores each metric across regions
 - builds region×region Pearson correlation matrix (diagonal set to 0)
 - writes to CSV (no headers, no index)
"""

import os
import argparse
import pandas as pd
import numpy as np
from scipy.stats import zscore

# the eight metrics we want from the stats files
METRICS = [
    "SurfArea",
    "GrayVol",
    "ThickAvg",
    "ThickStd",
    "MeanCurv",
    "GausCurv",
    "FoldInd",
    "CurvInd",
]

def parse_stats_file(fname):
    """
    Read a FreeSurfer .stats file and return a DataFrame with
    columns: Region, SurfArea, GrayVol, ..., CurvInd
    """
    # Most .stats files have a header line like:
    # # ColHeaders  StructName  NumVert  SurfArea  GrayVol  ThickAvg  ...
    # and data lines:
    # bankssts  12345  45.67  890.1  2.345  0.123  ...
    df = pd.read_csv(
        fname,
        delim_whitespace=True,
        comment="#",
        header=None,
        names=["Region", "NumVert"] + METRICS
    )
    # drop NumVert, keep Region + our metrics
    return df[["Region"] + METRICS]

def compute_msn_for_session(stats_dir, out_csv):
    """
    stats_dir: path to a freesurfer session stats/ folder
    out_csv: where to write the N×N matrix
    """
    # find all *.stats
    files = os.listdir(stats_dir)
    dfs = []
    for f in files:
        if not f.endswith(".stats"):
            continue
        # only take left/right hem files
        if ("lh." in f) or ("rh." in f):
            full = os.path.join(stats_dir, f)
            dfs.append(parse_stats_file(full))

    if not dfs:
        raise RuntimeError(f"No lh/rh .stats found in {stats_dir}")

    # stack left + right
    df = pd.concat(dfs, ignore_index=True)
    df = df.set_index("Region")

    # z-score each metric across regions (axis=0)
    # nan_policy='omit' in case of missing values
    mat = df.apply(lambda col: zscore(col, nan_policy="omit"), axis=0).values

    # compute Pearson correlation between rows (regions)
    corr = np.corrcoef(mat)

    # zero out the diagonal
    np.fill_diagonal(corr, 0.0)

    # write to CSV without headers/index
    pd.DataFrame(corr).to_csv(out_csv, header=False, index=False)

def main(root, outdir):
    for subj in sorted(os.listdir(root)):
        subjdir = os.path.join(root, subj)
        if not os.path.isdir(subjdir):
            continue

        for ses in sorted(os.listdir(subjdir)):
            sesdir = os.path.join(subjdir, ses)
            statsdir = os.path.join(sesdir, "stats")
            if not os.path.isdir(statsdir):
                continue

            # prepare output path
            od = os.path.join(outdir, subj)
            os.makedirs(od, exist_ok=True)
            out_csv = os.path.join(od, f"{ses}_msn.csv")

            print(f"→ {subj}/{ses} → {out_csv}")
            compute_msn_for_session(statsdir, out_csv)

if __name__ == "__main__":
    p = argparse.ArgumentParser(__doc__)
    p.add_argument(
        "--root",
        required=True,
        help="Root directory containing sub-***/ses-*/stats"
    )
    p.add_argument(
        "--outdir",
        required=True,
        help="Where to save each subj/ses MSN CSV"
    )
    args = p.parse_args()
    main(args.root, args.outdir)
