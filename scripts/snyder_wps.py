#!/usr/bin/env python3
"""
Snyder et al. cfDNA nucleosome footprinting utilities.

Implements the core Windowed Protection Score (WPS) from paired-end cfDNA BAMs:
for each genomic coordinate, WPS is the number of fragments spanning a centered
window minus the number of fragment endpoints inside that window.

This script is intentionally region-oriented. Whole-genome base-resolution WPS is
large; start with BED panels such as TSS windows, DHS regions, marker blocks, or
other loci of interest.
"""

import argparse
import csv
import math
import os
import sys
from collections import defaultdict
from dataclasses import dataclass

import numpy as np

try:
    import pysam
except ImportError:  # pragma: no cover - exercised only on missing dependency
    pysam = None


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int
    end: int
    name: str = "."
    score: str = "."
    strand: str = "."


@dataclass(frozen=True)
class TssWindow:
    region: Region
    tss: int


@dataclass(frozen=True)
class Fragment:
    start: int
    end: int


def parse_bed(path):
    """Yield BED3+ records as Region objects."""
    with open(path, "r", newline="") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_no}: expected at least 3 BED columns")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            if end <= start:
                raise ValueError(f"{path}:{line_no}: BED end must be greater than start")
            name = fields[3] if len(fields) > 3 else f"{chrom}:{start}-{end}"
            score = fields[4] if len(fields) > 4 else "."
            strand = fields[5] if len(fields) > 5 else "."
            yield Region(chrom, start, end, name, score, strand)


def parse_region_string(region):
    """Parse chr:start-end or chr:start-end:name into a Region."""
    label = "."
    if region.count(":") > 1:
        region, label = region.rsplit(":", 1)
    chrom, coords = region.split(":", 1)
    start_s, end_s = coords.replace(",", "").split("-", 1)
    start = int(start_s)
    end = int(end_s)
    if end <= start:
        raise ValueError("region end must be greater than start")
    return Region(chrom, start, end, label if label != "." else region)


def add_diff_range(diff, start, end, delta):
    """Add delta to an inclusive zero-based range in a difference array."""
    if end < start:
        return
    start = max(0, start)
    end = min(len(diff) - 2, end)
    if end < start:
        return
    diff[start] += delta
    diff[end + 1] -= delta


def wps_arrays_for_region(region, fragments, window_size=120):
    """
    Compute per-base WPS arrays for one region from fragment intervals.

    Coordinates follow BED convention. Output index i corresponds to genomic
    coordinate region.start + i, interpreted as the center of the WPS window.
    For even window sizes, the centered window is [x-left, x+right), where
    left = window_size // 2 and right = window_size - left.
    """
    if window_size <= 0:
        raise ValueError("window_size must be positive")
    length = region.end - region.start
    left = window_size // 2
    right = window_size - left
    span_diff = np.zeros(length + 1, dtype=np.int32)
    endpoint_diff = np.zeros(length + 1, dtype=np.int32)

    for frag in fragments:
        if frag.end <= frag.start:
            continue

        # Fragment spans centered window [x-left, x+right).
        span_center_start = frag.start + left
        span_center_end = frag.end - right
        add_diff_range(
            span_diff,
            span_center_start - region.start,
            span_center_end - region.start,
            1,
        )

        # Endpoint coordinates are the first and last covered bases.
        for endpoint in (frag.start, frag.end - 1):
            hit_center_start = endpoint - right + 1
            hit_center_end = endpoint + left
            add_diff_range(
                endpoint_diff,
                hit_center_start - region.start,
                hit_center_end - region.start,
                1,
            )

    span = np.cumsum(span_diff[:-1])
    endpoints = np.cumsum(endpoint_diff[:-1])
    return span - endpoints, span, endpoints


def iter_fragments_from_bam(
    bam,
    chrom,
    start,
    end,
    min_mapq=30,
    min_frag=50,
    max_frag=1000,
    include_duplicates=False,
):
    """Yield one Fragment per paired-end template from an indexed BAM."""
    fetch_start = max(0, start)
    fetch_end = max(fetch_start + 1, end)
    try:
        reads = bam.fetch(chrom, fetch_start, fetch_end)
    except ValueError as exc:
        raise ValueError(
            f"cannot fetch {chrom}:{fetch_start}-{fetch_end}; check BAM index and chromosome names"
        ) from exc

    for read in reads:
        if read.is_unmapped or read.mate_is_unmapped:
            continue
        if read.is_secondary or read.is_supplementary or read.is_qcfail:
            continue
        if read.is_duplicate and not include_duplicates:
            continue
        if read.mapping_quality < min_mapq:
            continue
        if read.template_length <= 0:
            continue
        frag_len = abs(read.template_length)
        if frag_len < min_frag or frag_len > max_frag:
            continue
        frag_start = read.reference_start
        frag_end = frag_start + frag_len
        yield Fragment(frag_start, frag_end)


def fetch_region_fragments(bam, region, window_size, min_mapq, min_frag, max_frag, include_duplicates):
    """Fetch fragments whose spans or endpoints can affect a region's WPS."""
    pad = window_size + max_frag
    fetch_start = max(0, region.start - pad)
    fetch_end = region.end + pad
    return list(
        iter_fragments_from_bam(
            bam,
            region.chrom,
            fetch_start,
            fetch_end,
            min_mapq=min_mapq,
            min_frag=min_frag,
            max_frag=max_frag,
            include_duplicates=include_duplicates,
        )
    )


def write_wps(args):
    if pysam is None:
        raise RuntimeError("pysam is required for BAM input")
    regions = collect_regions(args)
    with pysam.AlignmentFile(args.bam, "rb") as bam, open(args.output, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(
            [
                "chrom",
                "start",
                "end",
                "region_name",
                "wps",
                "spanning_fragments",
                "endpoint_fragments",
            ]
        )
        for region in regions:
            fragments = fetch_region_fragments(
                bam,
                region,
                args.window_size,
                args.min_mapq,
                args.min_frag,
                args.max_frag,
                args.include_duplicates,
            )
            wps, span, endpoints = wps_arrays_for_region(region, fragments, args.window_size)
            for offset, value in enumerate(wps):
                pos = region.start + offset
                writer.writerow(
                    [
                        region.chrom,
                        pos,
                        pos + 1,
                        region.name,
                        int(value),
                        int(span[offset]),
                        int(endpoints[offset]),
                    ]
                )


def summarize_regions(args):
    if pysam is None:
        raise RuntimeError("pysam is required for BAM input")
    regions = collect_regions(args)
    with pysam.AlignmentFile(args.bam, "rb") as bam, open(args.output, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(
            [
                "chrom",
                "start",
                "end",
                "name",
                "n_positions",
                "n_fragments_fetched",
                "mean_wps",
                "median_wps",
                "min_wps",
                "max_wps",
                "mean_spanning_fragments",
                "mean_endpoint_fragments",
            ]
        )
        for region in regions:
            fragments = fetch_region_fragments(
                bam,
                region,
                args.window_size,
                args.min_mapq,
                args.min_frag,
                args.max_frag,
                args.include_duplicates,
            )
            wps, span, endpoints = wps_arrays_for_region(region, fragments, args.window_size)
            writer.writerow(
                [
                    region.chrom,
                    region.start,
                    region.end,
                    region.name,
                    len(wps),
                    len(fragments),
                    fmt_float(np.mean(wps)),
                    fmt_float(np.median(wps)),
                    int(np.min(wps)) if len(wps) else "nan",
                    int(np.max(wps)) if len(wps) else "nan",
                    fmt_float(np.mean(span)),
                    fmt_float(np.mean(endpoints)),
                ]
            )


def tss_profile(args):
    if pysam is None:
        raise RuntimeError("pysam is required for BAM input")
    tss_windows = []
    for item in parse_bed(args.tss_bed):
        if item.strand == "-":
            tss = item.end - 1
        else:
            tss = item.start
        start = max(0, tss - args.flank)
        end = tss + args.flank + 1
        region = Region(item.chrom, start, end, item.name, item.score, item.strand)
        tss_windows.append(TssWindow(region, tss))

    sums = defaultdict(float)
    counts = defaultdict(int)
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for window in tss_windows:
            region = window.region
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
            for offset, value in enumerate(wps):
                pos = region.start + offset
                rel = pos - window.tss
                if args.strand_aware and region.strand == "-":
                    rel = -rel
                if -args.flank <= rel <= args.flank:
                    sums[rel] += float(value)
                    counts[rel] += 1

    with open(args.output, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["relative_position", "mean_wps", "n_tss"])
        for rel in range(-args.flank, args.flank + 1):
            if counts[rel]:
                writer.writerow([rel, fmt_float(sums[rel] / counts[rel]), counts[rel]])
            else:
                writer.writerow([rel, "nan", 0])


def fmt_float(value):
    if value is None or not math.isfinite(float(value)):
        return "nan"
    return f"{float(value):.6g}"


def collect_regions(args):
    regions = []
    if getattr(args, "regions_bed", None):
        regions.extend(parse_bed(args.regions_bed))
    if getattr(args, "region", None):
        regions.extend(parse_region_string(item) for item in args.region)
    if not regions:
        raise ValueError("provide --regions-bed and/or --region")
    return regions


def add_bam_filters(parser):
    parser.add_argument("--bam", required=True, help="Coordinate-sorted, indexed paired-end BAM")
    parser.add_argument("--window-size", type=int, default=120, help="Centered WPS window size")
    parser.add_argument("--min-mapq", type=int, default=30)
    parser.add_argument("--min-frag", type=int, default=50)
    parser.add_argument("--max-frag", type=int, default=1000)
    parser.add_argument("--include-duplicates", action="store_true")
    parser.add_argument("-o", "--output", required=True)


def build_parser():
    parser = argparse.ArgumentParser(
        description="Snyder et al. cfDNA Windowed Protection Score utilities"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_wps = sub.add_parser("wps", help="write base-resolution WPS over BED/region loci")
    add_bam_filters(p_wps)
    p_wps.add_argument("--regions-bed", help="BED loci to score")
    p_wps.add_argument("--region", action="append", help="Single chr:start-end locus; repeatable")
    p_wps.set_defaults(func=write_wps)

    p_summary = sub.add_parser("region-summary", help="write per-region WPS summary features")
    add_bam_filters(p_summary)
    p_summary.add_argument("--regions-bed", help="BED loci to summarize")
    p_summary.add_argument("--region", action="append", help="Single chr:start-end locus; repeatable")
    p_summary.set_defaults(func=summarize_regions)

    p_tss = sub.add_parser("tss-profile", help="aggregate WPS around BED6 transcription starts")
    add_bam_filters(p_tss)
    p_tss.add_argument("--tss-bed", required=True, help="BED6 genes/promoters; strand is used")
    p_tss.add_argument("--flank", type=int, default=1000)
    p_tss.add_argument(
        "--no-strand-aware",
        dest="strand_aware",
        action="store_false",
        help="Do not reverse minus-strand loci",
    )
    p_tss.set_defaults(func=tss_profile, strand_aware=True)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not os.path.exists(args.bam):
        parser.error(f"BAM does not exist: {args.bam}")
    try:
        args.func(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
