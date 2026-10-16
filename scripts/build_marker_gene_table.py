#!/usr/bin/env python3
"""Build marker_gene TSV from a gene-expression matrix and sample metadata.

Input expression matrix is wide: one row per gene, one column per sample.
The script performs one-vs-rest scoring within the provided cohort/source.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def choose_gene_col(df, requested):
    if requested:
        if requested not in df.columns:
            raise ValueError(f"gene column not found: {requested}")
        return requested
    for col in ["gene_symbol", "gene", "symbol", "Gene", "GeneSymbol", "gene_name"]:
        if col in df.columns:
            return col
    return df.columns[0]


def weighted_score(log2fc, group_mean, t_stat):
    expr_term = np.log2(group_mean + 1.0)
    t_term = np.nan_to_num(t_stat, nan=0.0, posinf=0.0, neginf=0.0)
    raw = np.maximum(log2fc, 0.0) * np.maximum(expr_term, 0.0) * np.maximum(t_term, 0.0)
    if len(raw) == 0 or np.nanmax(raw) <= 0:
        return np.zeros_like(raw)
    return 1000.0 * raw / np.nanmax(raw)


def welch_t(group, rest):
    g_mean = group.mean(axis=1)
    r_mean = rest.mean(axis=1)
    g_var = group.var(axis=1, ddof=1) if group.shape[1] > 1 else pd.Series(0, index=group.index)
    r_var = rest.var(axis=1, ddof=1) if rest.shape[1] > 1 else pd.Series(0, index=rest.index)
    denom = np.sqrt(g_var / max(group.shape[1], 1) + r_var / max(rest.shape[1], 1))
    return (g_mean - r_mean) / denom.replace(0, np.nan)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--expression", required=True, help="wide gene x sample expression TSV/CSV")
    ap.add_argument("--metadata", required=True, help="sample metadata TSV/CSV")
    ap.add_argument("--source", required=True, help="GTEx, TCGA, or another source label")
    ap.add_argument("--group-col", required=True, help="tissue/cancer-type column in metadata")
    ap.add_argument("--sample-id-col", default="sample_id")
    ap.add_argument("--gene-col", default=None)
    ap.add_argument("--sep", default="\t")
    ap.add_argument("--metadata-sep", default="\t")
    ap.add_argument("--top-n", type=int, default=100)
    ap.add_argument("--min-group-mean", type=float, default=1.0)
    ap.add_argument("--min-log2fc", type=float, default=1.0)
    ap.add_argument("--min-samples", type=int, default=3)
    ap.add_argument("--pseudocount", type=float, default=0.25)
    ap.add_argument("--allowed-groups", default=None, help="comma-separated subset of groups")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    expr = pd.read_csv(args.expression, sep=args.sep)
    meta = pd.read_csv(args.metadata, sep=args.metadata_sep, dtype=str)
    gene_col = choose_gene_col(expr, args.gene_col)
    genes = expr[gene_col].astype(str).str.strip()
    expr = expr.drop(columns=[gene_col]).apply(pd.to_numeric, errors="coerce")
    expr.index = genes
    expr = expr[~expr.index.duplicated(keep="first")]

    if args.sample_id_col not in meta.columns:
        raise ValueError(f"sample id column not found in metadata: {args.sample_id_col}")
    if args.group_col not in meta.columns:
        raise ValueError(f"group column not found in metadata: {args.group_col}")
    meta[args.sample_id_col] = meta[args.sample_id_col].astype(str).str.strip()
    meta[args.group_col] = meta[args.group_col].astype(str).str.strip()

    samples = [s for s in meta[args.sample_id_col] if s in expr.columns]
    if not samples:
        raise ValueError("no metadata samples matched expression matrix columns")
    expr = expr[samples]
    meta = meta[meta[args.sample_id_col].isin(samples)].copy()

    groups = sorted(meta[args.group_col].dropna().unique())
    if args.allowed_groups:
        allowed = set(x.strip() for x in args.allowed_groups.split(",") if x.strip())
        groups = [g for g in groups if g in allowed]

    rows = []
    for group_name in groups:
        group_samples = meta.loc[meta[args.group_col] == group_name, args.sample_id_col].tolist()
        rest_samples = meta.loc[meta[args.group_col] != group_name, args.sample_id_col].tolist()
        if len(group_samples) < args.min_samples or len(rest_samples) < args.min_samples:
            continue
        group_expr = expr[group_samples]
        rest_expr = expr[rest_samples]
        group_mean = group_expr.mean(axis=1)
        rest_mean = rest_expr.mean(axis=1)
        log2fc = np.log2((group_mean + args.pseudocount) / (rest_mean + args.pseudocount))
        t_stat = welch_t(group_expr, rest_expr)
        candidates = pd.DataFrame({
            "source": args.source,
            "marker_set": str(group_name).strip(),
            "gene_symbol": expr.index,
            "group_mean": group_mean.values,
            "rest_mean": rest_mean.values,
            "log2fc": log2fc.values,
            "t_stat": t_stat.values,
            "n_group": len(group_samples),
            "n_rest": len(rest_samples),
        })
        candidates = candidates[
            (candidates["group_mean"] >= args.min_group_mean)
            & (candidates["log2fc"] >= args.min_log2fc)
        ].copy()
        if candidates.empty:
            continue
        candidates["weight"] = weighted_score(candidates["log2fc"].values, candidates["group_mean"].values, candidates["t_stat"].values)
        candidates = candidates.sort_values(["log2fc", "group_mean", "t_stat"], ascending=[False, False, False]).head(args.top_n)
        candidates.insert(3, "rank", np.arange(1, len(candidates) + 1))
        rows.append(candidates)

    if not rows:
        raise ValueError("no markers passed filters")
    out = pd.concat(rows, ignore_index=True)
    cols = ["source", "marker_set", "gene_symbol", "rank", "weight", "log2fc", "group_mean", "rest_mean", "t_stat", "n_group", "n_rest"]
    out = out[cols]
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False, float_format="%.6g")
    print(f"wrote {len(out)} markers across {out[['source','marker_set']].drop_duplicates().shape[0]} marker sets -> {out_path}")
    print(out.groupby(["source", "marker_set"]).size().head(50).to_string())


if __name__ == "__main__":
    main()
