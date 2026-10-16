"""Balanced marker selection by one-versus-rest specificity."""
import argparse
import numpy as np
import pandas as pd
from .v1 import META_COLS

def compute_specificity_scores(A, celltypes):
    """One-vs-rest z-score for each (block, cell_type) pair."""
    n_blocks, n_types = A.shape
    scores = np.zeros_like(A)

    for t in range(n_types):
        rest_mask = np.ones(n_types, dtype=bool)
        rest_mask[t] = False
        rest = A[:, rest_mask]
        rest_mean = rest.mean(axis=1)
        rest_std = rest.std(axis=1, ddof=0)
        rest_std[rest_std < 1e-6] = 1e-6
        scores[:, t] = (A[:, t] - rest_mean) / rest_std

    best_type = np.argmax(np.abs(scores), axis=1)
    best_score = np.max(np.abs(scores), axis=1)
    return scores, best_type, best_score


def select_blocks(best_type, best_score, n_types, k):
    """Select top-K blocks per cell type, balanced."""
    selected = []
    for t in range(n_types):
        mask = best_type == t
        indices = np.where(mask)[0]
        type_scores = best_score[indices]
        top_k = indices[np.argsort(-type_scores)[:k]]
        selected.extend(top_k)
    return np.sort(np.array(selected))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--per-cell-type", type=int, default=100)
    args = parser.parse_args()
    if args.per_cell_type < 1:
        parser.error("--per-cell-type must be positive")
    frame = pd.read_csv(args.atlas, sep="\t")
    celltypes = [c for c in frame.columns if c not in META_COLS]
    if len(celltypes) < 2:
        parser.error("At least two cell types are required")
    values = frame[celltypes].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        parser.error("Reference values must be finite")
    _, owner, score = compute_specificity_scores(values, celltypes)
    selected = select_blocks(owner, score, len(celltypes), args.per_cell_type)
    frame.iloc[selected].to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
