"""Coverage-weighted methylation deconvolution.

The extraction, weighting, NNLS, normalization, R2 and minimum-block logic
retain the baseline implementation. Inputs are supplied by the caller.
Author: Yujin Jeon.
"""
import numpy as np
import pandas as pd
from scipy import optimize

META_COLS = ["chr", "start", "end", "startCpG", "endCpG", "n_cpgs"]
BLOOD_TYPES = [
    "Blood-B", "Blood-Granul", "Blood-Mono+Macro",
    "Blood-NK", "Blood-T", "Eryth-prog",
]

def extract_block_amf_from_bedgraph(bg_path, blocks_df):
    """Single-pointer sweep intersection of per-CpG bedGraph against blocks.
    Returns (betas, covs) per block. AMF = sum(meth) / sum(meth+unmeth)."""
    n_blocks = len(blocks_df)
    block_meth = np.zeros(n_blocks, dtype=np.float64)
    block_total = np.zeros(n_blocks, dtype=np.float64)
    blk_chr = blocks_df['chr'].values
    blk_start = blocks_df['start'].values
    blk_end = blocks_df['end'].values

    chr_blocks = {}
    for i in range(n_blocks):
        c = blk_chr[i]
        if c not in chr_blocks:
            chr_blocks[c] = []
        chr_blocks[c].append(i)
    for c in chr_blocks:
        indices = chr_blocks[c]
        indices.sort(key=lambda i: blk_start[i])
        chr_blocks[c] = indices

    with open(bg_path, 'r') as f:
        cur_chr = None
        cur_indices = []
        cur_ptr = 0
        for line in f:
            if line.startswith('track'):
                continue
            parts = line.rstrip('\n').split('\t')
            chrom = parts[0]
            pos = int(parts[1])
            m = int(parts[4])
            u = int(parts[5])
            total = m + u
            if chrom != cur_chr:
                cur_chr = chrom
                cur_indices = chr_blocks.get(chrom, [])
                cur_ptr = 0
            while cur_ptr < len(cur_indices) and blk_end[cur_indices[cur_ptr]] <= pos:
                cur_ptr += 1
            if cur_ptr >= len(cur_indices):
                continue
            bi = cur_indices[cur_ptr]
            if pos >= blk_start[bi] and pos < blk_end[bi]:
                block_meth[bi] += m
                block_total[bi] += total

    betas = np.full(n_blocks, np.nan)
    covs = np.zeros(n_blocks)
    has_cov = block_total > 0
    betas[has_cov] = block_meth[has_cov] / block_total[has_cov]
    covs[has_cov] = block_total[has_cov]
    return betas, covs


def _nnls_sqrt_cov(A, cov, s):
    w = np.sqrt(cov)
    coef, _ = optimize.nnls(A * w[:, None], s * w)
    return coef


def _r2_sqrt_cov(A, cov, s, coef):
    w = np.sqrt(cov)
    sw = s * w
    pred = (A * w[:, None]) @ coef
    ssr = float(np.sum((sw - pred) ** 2))
    sst = float(np.sum((sw - sw.mean()) ** 2))
    return 1.0 - ssr / sst if sst > 0 else 0.0


class V1Model:
    def __init__(self, atlas_path, config=None):
        self.cfg = dict(config or {})
        self.min_cov = int(self.cfg.get("min_cov", 3))
        if self.min_cov < 1:
            raise ValueError("min_cov must be positive")
        self.blood_types = self.cfg.get("blood_types", BLOOD_TYPES)
        meta_cols = self.cfg.get("meta_cols", META_COLS)
        df = pd.read_csv(atlas_path, sep='\t')
        self.celltypes = [c for c in df.columns if c not in meta_cols]
        # sanity: config-pinned order must match atlas header
        if "celltypes" in self.cfg and self.cfg["celltypes"] != self.celltypes:
            raise ValueError("atlas celltype order does not match config-pinned order")
        if not self.celltypes:
            raise ValueError("Atlas has no cell-type columns")
        self.A = df[self.celltypes].values.astype(float)
        self.blocks = df[meta_cols].copy()
        self._finite = np.all(np.isfinite(self.A), axis=1)

    def deconv(self, betas, covs):
        mask = (np.isfinite(betas) & (covs >= self.min_cov) & self._finite)
        n_blocks = int(mask.sum())
        if n_blocks < 50:
            return None
        s = betas[mask]
        c = covs[mask]
        Am = self.A[mask]
        coef = _nnls_sqrt_cov(Am, c, s)
        raw_sum = float(coef.sum())
        frac = coef / raw_sum if raw_sum > 0 else coef
        r2 = _r2_sqrt_cov(Am, c, s, coef)
        out = {ct: float(f) for ct, f in zip(self.celltypes, frac)}
        out['blood_total'] = float(sum(out.get(b, 0.0) for b in self.blood_types))
        out['raw_sum'] = raw_sum
        out['R2'] = r2
        out['n_blocks'] = n_blocks
        out['mean_cov'] = float(c.mean())
        return out

    def deconv_bedgraph(self, bg_path):
        betas, covs = extract_block_amf_from_bedgraph(bg_path, self.blocks)
        return self.deconv(betas, covs)
