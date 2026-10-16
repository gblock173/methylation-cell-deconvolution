"""Build a reference atlas using explicit, user-supplied manifests.

Baseline numeric and CpG slice conventions are preserved. This module uses
module-level settings and must be run as one process per atlas build.
"""
import argparse
import gzip
import random
from pathlib import Path
import numpy as np
import pandas as pd

N_CPG = 0
CPG_BED = None
CELLTYPES = []
CELLTYPE_MAP = {}
COVERAGE_FILES = []

def load_cpg_index():
    chrom_list, pos_list = [], []
    with (gzip.open(CPG_BED, 'rt') if str(CPG_BED).endswith('.gz') else open(CPG_BED, 'rt')) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            chrom_list.append(parts[0])
            pos_list.append(int(parts[1]))
    return chrom_list, pos_list


def read_beta_file(path):
    return np.memmap(path, dtype=np.uint8, mode='r').reshape(-1, 2)


def step1_filter_cpgs(file_idx, min_range=0.3):
    print(f"\n{'='*60}")
    print(f"STEP 1: Filter CpGs (range threshold = {min_range})")
    print(f"{'='*60}")

    beta_max = np.full(N_CPG, -1.0, dtype=np.float32)
    beta_min = np.full(N_CPG, 2.0, dtype=np.float32)
    valid_count = np.zeros(N_CPG, dtype=np.uint8)

    for ci, ct in enumerate(CELLTYPES):
        reference_groups = CELLTYPE_MAP[ct]
        rep_files = []
        for ln in reference_groups:
            rep_files.extend(file_idx.get(ln, []))
        if not rep_files:
            print(f"  WARNING: no .beta files for {ct}")
            continue

        print(f"  [{ci+1:2d}/{len(CELLTYPES)}] {ct}: {len(rep_files)} replicates ...", end='', flush=True)
        meth = np.zeros(N_CPG, dtype=np.uint32)
        cov = np.zeros(N_CPG, dtype=np.uint32)
        for bf in rep_files:
            data = read_beta_file(bf)
            n = min(len(data), N_CPG)
            meth[:n] += data[:n, 0].astype(np.uint32)
            cov[:n] += data[:n, 1].astype(np.uint32)

        has_cov = cov > 0
        ct_beta = np.where(has_cov, meth.astype(np.float32) / cov.astype(np.float32), np.nan)
        valid_count[has_cov] += 1
        np.fmax(beta_max, ct_beta, out=beta_max, where=has_cov)
        np.fmin(beta_min, ct_beta, out=beta_min, where=has_cov)
        del meth, cov, ct_beta
        print(" done")

    valid_all = valid_count == len(CELLTYPES)
    beta_range = beta_max - beta_min
    passed = valid_all & (beta_range >= min_range)

    n_valid = int(valid_all.sum())
    n_passed = int(passed.sum())
    print(f"\n  CpGs with coverage in all {len(CELLTYPES)} types: {n_valid:,}")
    print(f"  CpGs with range >= {min_range}: {n_passed:,}")
    print(f"  Filtered out: {n_valid - n_passed:,}")

    del beta_max, beta_min, valid_count
    return passed


def step2_build_blocks(passed, min_cpgs=4):
    print(f"\n{'='*60}")
    print(f"STEP 2: Build contiguous blocks (min {min_cpgs} CpGs)")
    print(f"{'='*60}")

    print("  Loading CpG index ...")
    chroms, positions = load_cpg_index()

    blocks = []
    cur_chr = None
    run_start_idx = None
    run_count = 0

    for i in range(N_CPG):
        c = chroms[i]
        if c != cur_chr:
            if run_count >= min_cpgs and cur_chr is not None:
                blocks.append((cur_chr, positions[run_start_idx], positions[i-1] + 1,
                               run_start_idx + 1, i, run_count))
            cur_chr = c
            run_start_idx = None
            run_count = 0

        if passed[i]:
            if run_start_idx is None:
                run_start_idx = i
                run_count = 1
            else:
                run_count += 1
        else:
            if run_count >= min_cpgs:
                blocks.append((cur_chr, positions[run_start_idx], positions[i-1] + 1,
                               run_start_idx + 1, i, run_count))
            run_start_idx = None
            run_count = 0

    if run_count >= min_cpgs and cur_chr is not None:
        blocks.append((cur_chr, positions[run_start_idx], positions[N_CPG-1] + 1,
                       run_start_idx + 1, N_CPG, run_count))

    blocks_df = pd.DataFrame(blocks, columns=['chr', 'start', 'end', 'startCpG', 'endCpG', 'n_cpgs'])
    print(f"  Total blocks (>= {min_cpgs} CpGs): {len(blocks_df):,}")
    print(f"  Total CpGs in blocks: {blocks_df['n_cpgs'].sum():,}")
    print(f"  Block size distribution: "
          f"min={blocks_df['n_cpgs'].min()}, "
          f"median={blocks_df['n_cpgs'].median():.0f}, "
          f"max={blocks_df['n_cpgs'].max()}")

    del chroms, positions
    return blocks_df


def step3_coverage_filter(blocks_df, n_samples=10, min_depth=10, seed=42):
    print(f"\n{'='*60}")
    print(f"STEP 3: cfDNA coverage filter ({n_samples} samples, min depth {min_depth})")
    print(f"{'='*60}")

    cfdna_files = list(COVERAGE_FILES)
    print(f"  Found {len(cfdna_files)} cfDNA .beta files")

    random.seed(seed)
    selected = random.sample(cfdna_files, min(n_samples, len(cfdna_files)))
    print(f"  Selected {len(selected)} samples:")

    starts = blocks_df['startCpG'].values.astype(int)
    ends = blocks_df['endCpG'].values.astype(int)
    n_blocks = len(blocks_df)
    all_pass = np.ones(n_blocks, dtype=bool)

    for fi, bf in enumerate(selected):
        print(f"  [{fi+1}/{len(selected)}] checking coverage ...", end='', flush=True)
        data = read_beta_file(bf)
        dlen = len(data)
        for bi in range(n_blocks):
            if not all_pass[bi]:
                continue
            s, e = starts[bi], ends[bi]
            if e > dlen:
                all_pass[bi] = False
                continue
            total_cov = data[s:e, 1].astype(np.int64).sum()
            if total_cov < min_depth:
                all_pass[bi] = False
        n_surviving = int(all_pass.sum())
        print(f" {n_surviving}/{n_blocks} blocks pass")

    blocks_filtered = blocks_df[all_pass].reset_index(drop=True)
    print(f"\n  Blocks after coverage filter: {len(blocks_filtered):,} / {n_blocks:,}")
    return blocks_filtered


def step4_build_reference(blocks_df, file_idx):
    print(f"\n{'='*60}")
    print(f"STEP 4: Build AMF reference atlas")
    print(f"{'='*60}")

    starts = blocks_df['startCpG'].values.astype(int)
    ends = blocks_df['endCpG'].values.astype(int)
    n_blocks = len(blocks_df)

    ref = np.full((n_blocks, len(CELLTYPES)), np.nan)

    for ci, ct in enumerate(CELLTYPES):
        reference_groups = CELLTYPE_MAP[ct]
        rep_files = []
        for ln in reference_groups:
            rep_files.extend(file_idx.get(ln, []))
        if not rep_files:
            print(f"  WARNING: no files for {ct}")
            continue

        print(f"  [{ci+1:2d}/{len(CELLTYPES)}] {ct}: {len(rep_files)} reps ...", end='', flush=True)
        agg_meth = np.zeros(n_blocks, dtype=np.float64)
        agg_cov = np.zeros(n_blocks, dtype=np.float64)

        for bf in rep_files:
            data = read_beta_file(bf)
            dlen = len(data)
            for bi in range(n_blocks):
                s, e = starts[bi], ends[bi]
                if e <= dlen:
                    agg_meth[bi] += data[s:e, 0].astype(np.float64).sum()
                    agg_cov[bi] += data[s:e, 1].astype(np.float64).sum()

        has_cov = agg_cov > 0
        ref[has_cov, ci] = agg_meth[has_cov] / agg_cov[has_cov]
        n_nan = int((~has_cov).sum())
        print(f" done (NaN blocks: {n_nan})")

    atlas_df = blocks_df.copy()
    for ci, ct in enumerate(CELLTYPES):
        atlas_df[ct] = ref[:, ci]

    return atlas_df


def main():
    global N_CPG, CPG_BED, CELLTYPES, CELLTYPE_MAP, COVERAGE_FILES
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-manifest", required=True,
                        help="TSV columns: cell_type, beta_file")
    parser.add_argument("--coverage-manifest", required=True,
                        help="TSV column: beta_file")
    parser.add_argument("--cpg-index", required=True,
                        help="CpG chromosome and coordinate TSV, optionally gzip")
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-range", type=float, default=0.3)
    parser.add_argument("--min-cpgs", type=int, default=4)
    parser.add_argument("--min-depth", type=int, default=10)
    parser.add_argument("--coverage-samples", type=int, default=10)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    if not 0 <= args.min_range <= 1:
        parser.error("--min-range must be between zero and one")
    if min(args.min_cpgs, args.min_depth, args.coverage_samples) < 1:
        parser.error("Count and depth settings must be positive")
    refs = pd.read_csv(args.reference_manifest, sep="\t")
    coverage = pd.read_csv(args.coverage_manifest, sep="\t")
    if not {"cell_type", "beta_file"} <= set(refs.columns):
        parser.error("Reference manifest requires cell_type and beta_file")
    if "beta_file" not in coverage.columns or coverage.empty:
        parser.error("Coverage manifest requires at least one beta_file")
    if refs.empty or refs[["cell_type", "beta_file"]].isna().any().any():
        parser.error("Reference manifest must have complete entries")
    if coverage["beta_file"].isna().any():
        parser.error("Coverage manifest must have complete entries")
    CELLTYPES = sorted(refs.cell_type.unique())
    if not 2 <= len(CELLTYPES) <= 255:
        parser.error("Expected between 2 and 255 cell types")
    CELLTYPE_MAP = {ct: [ct] for ct in CELLTYPES}
    file_idx = {ct: list(group.beta_file) for ct, group in refs.groupby("cell_type")}
    COVERAGE_FILES = list(coverage.beta_file)
    CPG_BED = args.cpg_index
    chroms, positions = load_cpg_index()
    N_CPG = len(chroms)
    del chroms, positions
    if N_CPG < args.min_cpgs:
        parser.error("CpG index is too short")
    for beta_file in list(refs.beta_file) + COVERAGE_FILES:
        if Path(beta_file).stat().st_size != N_CPG * 2:
            parser.error("Every beta input must have one uint8 pair per indexed CpG")
    passed = step1_filter_cpgs(file_idx, args.min_range)
    blocks = step2_build_blocks(passed, args.min_cpgs)
    if blocks.empty:
        parser.error("No blocks passed the CpG filters")
    blocks = step3_coverage_filter(blocks, args.coverage_samples, args.min_depth, args.seed)
    if blocks.empty:
        parser.error("No blocks passed the coverage filter")
    step4_build_reference(blocks, file_idx).to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
