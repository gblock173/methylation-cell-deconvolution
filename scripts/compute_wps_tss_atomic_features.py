#!/usr/bin/env python3
"""Compute per-marker atomic WPS-TSS features plus weighted/rank-bin summaries."""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from snyder_wps import Region, fetch_region_fragments, wps_arrays_for_region  # noqa: E402

METRICS = [
    'center_wps', 'left_peak_wps', 'right_peak_wps', 'flank_wps',
    'center_depletion', 'phasing_amplitude', 'mean_wps',
    'mean_fragments_per_tss',
]


def safe_mean(values):
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(arr.mean()) if len(arr) else np.nan


def safe_weighted_mean(values, weights):
    v = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    mask = np.isfinite(v) & np.isfinite(w) & (w > 0)
    if not mask.any():
        return np.nan
    return float(np.average(v[mask], weights=w[mask]))


def direction_from_marker_set(marker_set):
    text = str(marker_set)
    if text.endswith('__hyper'):
        return 'hyper'
    if text.endswith('__hypo'):
        return 'hypo'
    return 'unspecified'


def rank_bin(rank):
    rank = int(rank)
    if rank <= 10:
        return 'rank01_10'
    if rank <= 25:
        return 'rank11_25'
    return 'rank26_50'


def load_marker_metadata(path):
    df = pd.read_csv(path, sep='\t')
    if 'primary_contig' in df.columns:
        df = df[df['primary_contig'].astype(bool)].copy()
    if 'mapped' in df.columns:
        df = df[df['mapped'].astype(bool)].copy()
    required = ['source', 'marker_set', 'gene_symbol_marker', 'rank', 'weight', 'chrom', 'tss0', 'tss1', 'strand']
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f'metadata missing required columns: {missing}')
    df['rank'] = pd.to_numeric(df['rank'], errors='raise').astype(int)
    df['weight'] = pd.to_numeric(df['weight'], errors='coerce').fillna(0.0).astype(float)
    df['tss0'] = pd.to_numeric(df['tss0'], errors='raise').astype(int)
    df['tss1'] = pd.to_numeric(df['tss1'], errors='raise').astype(int)
    df['direction'] = df['marker_set'].map(direction_from_marker_set)
    df['rank_bin'] = df['rank'].map(rank_bin)
    return df.sort_values(['source', 'marker_set', 'rank', 'gene_symbol_marker']).reset_index(drop=True)


def tss_for_row(row):
    return int(row.tss1) - 1 if row.strand == '-' else int(row.tss0)


def tss_summary(bam, row, args):
    tss = tss_for_row(row)
    start = max(0, tss - args.flank)
    end = tss + args.flank + 1
    region = Region(str(row.chrom), start, end, str(row.gene_symbol_marker), int(round(row.weight)), str(row.strand))
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
    if row.strand == '-':
        rel = -rel
    center = np.abs(rel) <= args.center_bp
    left_peak = (rel >= -args.peak_outer) & (rel <= -args.peak_inner)
    right_peak = (rel >= args.peak_inner) & (rel <= args.peak_outer)

    center_wps = safe_mean(wps[center])
    left_wps = safe_mean(wps[left_peak])
    right_wps = safe_mean(wps[right_peak])
    flank_wps = safe_mean([left_wps, right_wps])
    peak_wps = np.nanmax([left_wps, right_wps]) if np.isfinite([left_wps, right_wps]).any() else np.nan
    n_fragments = int(len(fragments))
    return {
        'center_wps': center_wps,
        'left_peak_wps': left_wps,
        'right_peak_wps': right_wps,
        'flank_wps': flank_wps,
        'center_depletion': flank_wps - center_wps if np.isfinite(flank_wps) and np.isfinite(center_wps) else np.nan,
        'phasing_amplitude': peak_wps - center_wps if np.isfinite(peak_wps) and np.isfinite(center_wps) else np.nan,
        'mean_wps': safe_mean(wps),
        'n_positions': int(len(wps)),
        'n_fragments': n_fragments,
        'mean_fragments_per_tss': float(n_fragments),
    }


def summarize_group(sample_id, rows, aggregation, rank_bin_label):
    first = rows.iloc[0]
    out = {
        'sample_id': sample_id,
        'source': first['source'],
        'marker_set': first['marker_set'],
        'direction': first['direction'],
        'aggregation': aggregation,
        'rank_bin': rank_bin_label,
        'n_tss': int(len(rows)),
        'weight_sum': float(rows['weight'].sum()),
    }
    for metric in METRICS:
        if aggregation == 'weighted':
            out[metric] = safe_weighted_mean(rows[metric], rows['weight'])
        else:
            out[metric] = safe_mean(rows[metric])
    out['n_positions'] = int(rows['n_positions'].sum())
    out['n_fragments'] = int(rows['n_fragments'].sum())
    return out


def build_summary(sample_id, atomic):
    rows = []
    group_cols = ['source', 'marker_set']
    for _, group in atomic.groupby(group_cols, sort=True):
        for aggregation in ['unweighted', 'weighted']:
            rows.append(summarize_group(sample_id, group, aggregation, 'all'))
            for rb, sub in group.groupby('rank_bin', sort=True):
                rows.append(summarize_group(sample_id, sub, aggregation, rb))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--bam', required=True)
    ap.add_argument('--marker-metadata', required=True)
    ap.add_argument('--sample-id', required=True)
    ap.add_argument('--atomic-output', required=True)
    ap.add_argument('--summary-output', required=True)
    ap.add_argument('--window-size', type=int, default=120)
    ap.add_argument('--flank', type=int, default=1000)
    ap.add_argument('--center-bp', type=int, default=50)
    ap.add_argument('--peak-inner', type=int, default=120)
    ap.add_argument('--peak-outer', type=int, default=250)
    ap.add_argument('--min-mapq', type=int, default=30)
    ap.add_argument('--min-frag', type=int, default=50)
    ap.add_argument('--max-frag', type=int, default=1000)
    ap.add_argument('--include-duplicates', action='store_true')
    args = ap.parse_args()

    if not os.path.exists(args.bam):
        raise FileNotFoundError(args.bam)
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError('pysam is required for BAM input') from exc

    markers = load_marker_metadata(args.marker_metadata)
    atomic_rows = []
    with pysam.AlignmentFile(args.bam, 'rb') as bam:
        for row in markers.itertuples(index=False):
            rec = tss_summary(bam, row, args)
            rec.update({
                'sample_id': args.sample_id,
                'source': row.source,
                'marker_set': row.marker_set,
                'direction': row.direction,
                'rank_bin': row.rank_bin,
                'gene_symbol': row.gene_symbol_marker,
                'rank': int(row.rank),
                'weight': float(row.weight),
                'chrom': row.chrom,
                'tss0': int(row.tss0),
                'tss1': int(row.tss1),
                'strand': row.strand,
            })
            atomic_rows.append(rec)
    atomic = pd.DataFrame(atomic_rows)
    atomic_cols = [
        'sample_id', 'source', 'marker_set', 'direction', 'rank_bin', 'gene_symbol',
        'rank', 'weight', 'chrom', 'tss0', 'tss1', 'strand',
        'center_wps', 'left_peak_wps', 'right_peak_wps', 'flank_wps',
        'center_depletion', 'phasing_amplitude', 'mean_wps', 'n_positions',
        'n_fragments', 'mean_fragments_per_tss',
    ]
    atomic = atomic[atomic_cols]
    summary = build_summary(args.sample_id, atomic)

    atomic_path = Path(args.atomic_output)
    summary_path = Path(args.summary_output)
    atomic_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    atomic.to_csv(atomic_path, sep='\t', index=False, float_format='%.6g', compression='infer')
    summary.to_csv(summary_path, sep='\t', index=False, float_format='%.6g')
    print(f'wrote atomic {atomic.shape[0]} rows -> {atomic_path}')
    print(f'wrote summary {summary.shape[0]} rows -> {summary_path}')


if __name__ == '__main__':
    main()
