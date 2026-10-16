#!/usr/bin/env python3
"""Compute marker-set WPS TSS summary features for one cfDNA BAM."""

import argparse
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from snyder_wps import Region, fetch_region_fragments, parse_bed, wps_arrays_for_region  # noqa: E402


def parse_marker_name(name):
    parts = str(name).split("|")
    if len(parts) >= 3:
        return parts[0], parts[1], "|".join(parts[2:])
    if len(parts) == 2:
        return "marker", parts[0], parts[1]
    return "marker", "all", str(name)


def tss_for_bed6(region):
    return region.end - 1 if region.strand == "-" else region.start


def safe_mean(values):
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(arr.mean()) if len(arr) else np.nan


def tss_summary(bam, item, args):
    tss = tss_for_bed6(item)
    start = max(0, tss - args.flank)
    end = tss + args.flank + 1
    region = Region(item.chrom, start, end, item.name, item.score, item.strand)
    fragments = fetch_region_fragments(
        bam,
        region,
        args.window_size,
        args.min_mapq,
        args.min_frag,
        args.max_frag,
        args.include_duplicates,
    )
    wps, _, _ = wps_arrays_for_region(region, fragments, args.window_size)
    rel = np.arange(region.start, region.end) - tss
    if item.strand == "-":
        rel = -rel
    center = np.abs(rel) <= args.center_bp
    left_peak = (rel >= -args.peak_outer) & (rel <= -args.peak_inner)
    right_peak = (rel >= args.peak_inner) & (rel <= args.peak_outer)
    full = np.ones_like(rel, dtype=bool)

    center_wps = safe_mean(wps[center])
    left_wps = safe_mean(wps[left_peak])
    right_wps = safe_mean(wps[right_peak])
    flank_wps = safe_mean([left_wps, right_wps])
    peak_wps = np.nanmax([left_wps, right_wps]) if np.isfinite([left_wps, right_wps]).any() else np.nan
    return {
        "center_wps": center_wps,
        "left_peak_wps": left_wps,
        "right_peak_wps": right_wps,
        "flank_wps": flank_wps,
        "center_depletion": flank_wps - center_wps if np.isfinite(flank_wps) and np.isfinite(center_wps) else np.nan,
        "phasing_amplitude": peak_wps - center_wps if np.isfinite(peak_wps) and np.isfinite(center_wps) else np.nan,
        "mean_wps": safe_mean(wps[full]),
        "n_positions": int(len(wps)),
        "n_fragments": int(len(fragments)),
    }


def summarize_group(sample_id, source, marker_set, rows):
    out = {"sample_id": sample_id, "source": source, "marker_set": marker_set, "n_tss": len(rows)}
    metrics = [
        "center_wps",
        "left_peak_wps",
        "right_peak_wps",
        "flank_wps",
        "center_depletion",
        "phasing_amplitude",
        "mean_wps",
    ]
    for metric in metrics:
        out[metric] = safe_mean([r[metric] for r in rows])
    out["n_positions"] = int(np.sum([r["n_positions"] for r in rows]))
    out["n_fragments"] = int(np.sum([r["n_fragments"] for r in rows]))
    out["mean_fragments_per_tss"] = out["n_fragments"] / out["n_tss"] if out["n_tss"] else np.nan
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--tss-bed", required=True)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--window-size", type=int, default=120)
    ap.add_argument("--flank", type=int, default=1000)
    ap.add_argument("--center-bp", type=int, default=50)
    ap.add_argument("--peak-inner", type=int, default=120)
    ap.add_argument("--peak-outer", type=int, default=250)
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--min-frag", type=int, default=50)
    ap.add_argument("--max-frag", type=int, default=1000)
    ap.add_argument("--include-duplicates", action="store_true")
    args = ap.parse_args()

    if not os.path.exists(args.bam):
        raise FileNotFoundError(args.bam)
    if not os.path.exists(args.tss_bed):
        raise FileNotFoundError(args.tss_bed)

    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError("pysam is required for BAM input") from exc

    groups = defaultdict(list)
    items = list(parse_bed(args.tss_bed))
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for item in items:
            source, marker_set, gene = parse_marker_name(item.name)
            rec = tss_summary(bam, item, args)
            rec.update({"gene_symbol": gene})
            groups[(source, marker_set)].append(rec)

    rows = [summarize_group(args.sample_id, source, marker_set, vals) for (source, marker_set), vals in sorted(groups.items())]
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False, float_format="%.6g")
    print(f"wrote {len(rows)} marker-set rows -> {out_path}")


if __name__ == "__main__":
    main()
