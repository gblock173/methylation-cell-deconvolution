#!/usr/bin/env python3
"""Convert GTEx/TCGA marker genes into hg38 TSS BED6 from GTF/GFF3."""

import argparse
import gzip
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

ATTR_RE = re.compile(r'(\S+) \"([^\"]*)\"')


def open_maybe_gzip(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def parse_attrs(attr_text):
    if "=" in attr_text and not ATTR_RE.search(attr_text):
        attrs = {}
        for item in attr_text.strip().strip(";").split(";"):
            if not item or "=" not in item:
                continue
            key, value = item.split("=", 1)
            attrs[key] = value
        return attrs
    return dict(ATTR_RE.findall(attr_text))


def infer_format(path):
    suffixes = "".join(Path(str(path)).suffixes).lower()
    if "gff" in suffixes:
        return "gff3"
    if "gtf" in suffixes:
        return "gtf"
    return "auto"


def feature_aliases(feature, annotation_format):
    if feature == "transcript" and annotation_format == "gff3":
        return {"transcript", "mRNA"}
    return {feature}



def iter_candidate_lines(annotation_path, gene_keys, annotation_format):
    """Yield annotation lines, using grep to prefilter target gene names when possible."""
    path_s = str(annotation_path)
    can_grep = bool(gene_keys) and not path_s.endswith(".gz") and shutil.which("grep")
    if not can_grep:
        with open_maybe_gzip(annotation_path) as handle:
            for line in handle:
                yield line
        return

    patterns = []
    for gene in sorted(gene_keys):
        if annotation_format == "gff3":
            patterns.extend([f"gene_name={gene}", f"Name={gene}"])
        else:
            patterns.extend([f'gene_name "{gene}"', f'gene "{gene}"'])
    tmp_name = None
    try:
        with tempfile.NamedTemporaryFile("w", delete=False) as tmp:
            tmp_name = tmp.name
            tmp.write("\n".join(patterns))
            tmp.write("\n")
        proc = subprocess.Popen(
            ["grep", "-F", "-f", tmp_name, path_s],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            yield line
        _, stderr = proc.communicate()
        if proc.returncode not in (0, 1):
            raise RuntimeError(f"grep prefilter failed: {stderr.strip()}")
    finally:
        if tmp_name and os.path.exists(tmp_name):
            os.unlink(tmp_name)


def normalize_gene_symbol(value):
    if pd.isna(value):
        return None
    value = str(value).strip()
    return value.upper() if value else None


def marker_table(path, default_source):
    df = pd.read_csv(path, sep="\t", dtype=str)
    cols = {c.lower(): c for c in df.columns}
    gene_col = cols.get("gene_symbol") or cols.get("gene") or cols.get("symbol")
    set_col = cols.get("marker_set") or cols.get("tissue") or cols.get("cancer_type") or cols.get("type")
    source_col = cols.get("source")
    if gene_col is None or set_col is None:
        raise ValueError("marker TSV needs gene_symbol/gene/symbol and marker_set/tissue/cancer_type/type")
    out = pd.DataFrame({
        "source": df[source_col] if source_col else default_source,
        "marker_set": df[set_col],
        "gene_symbol": df[gene_col],
    })
    for opt in ["rank", "weight", "score", "pvalue", "qvalue"]:
        if opt in cols:
            out[opt] = df[cols[opt]]
    out["source"] = out["source"].fillna(default_source).astype(str).str.strip()
    out["marker_set"] = out["marker_set"].astype(str).str.strip()
    out["gene_symbol"] = out["gene_symbol"].astype(str).str.strip()
    out["gene_key"] = out["gene_symbol"].map(normalize_gene_symbol)
    out = out.dropna(subset=["gene_key"])
    out = out[out["marker_set"] != ""]
    return out.drop_duplicates(["source", "marker_set", "gene_key"])


def annotation_tss_table(annotation_path, feature="transcript", annotation_format="auto", gene_keys=None):
    if annotation_format == "auto":
        annotation_format = infer_format(annotation_path)
    wanted = feature_aliases(feature, annotation_format)
    gene_keys = set(gene_keys or [])
    rows = []
    for line in iter_candidate_lines(annotation_path, gene_keys, annotation_format):
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9 or parts[2] not in wanted:
            continue
        chrom, _, feat, start_s, end_s, _, strand, _, attrs_s = parts
        if strand not in {"+", "-"}:
            continue
        attrs = parse_attrs(attrs_s)
        gene_name = attrs.get("gene_name") or attrs.get("gene") or attrs.get("Name")
        gene_id = attrs.get("gene_id") or attrs.get("Parent") or attrs.get("ID")
        gene_key = normalize_gene_symbol(gene_name or gene_id)
        if gene_keys and gene_key not in gene_keys:
            continue
        transcript_id = attrs.get("transcript_id") or attrs.get("ID") or "."
        if gene_id and gene_id.startswith("gene:"):
            gene_id = gene_id.split(":", 1)[1]
        if transcript_id and transcript_id.startswith("transcript:"):
            transcript_id = transcript_id.split(":", 1)[1]
        if not gene_name and not gene_id:
            continue
        start = int(start_s)
        end = int(end_s)
        tss0 = start - 1 if strand == "+" else end - 1
        tags = attrs.get("tag", "")
        canonical = int(
            "MANE_Select" in tags
            or "Ensembl_canonical" in tags
            or attrs.get("transcript_support_level") == "1"
        )
        rows.append({
            "chrom": chrom,
            "tss0": tss0,
            "tss1": tss0 + 1,
            "strand": strand,
            "annotation_feature": feat,
            "gene_symbol": gene_name or gene_id,
            "gene_id": gene_id,
            "gene_key": gene_key,
            "transcript_id": transcript_id,
            "gene_type": attrs.get("gene_type") or attrs.get("gene_biotype") or attrs.get("biotype", "."),
            "canonical_score": canonical,
            "span": end - start + 1,
        })
    if not rows:
        raise ValueError(f"no {feature} records parsed from {annotation_path}")
    return pd.DataFrame(rows)


def choose_representative_tss(gtf_df):
    gtf_df = gtf_df.sort_values(
        ["gene_key", "canonical_score", "span", "transcript_id"],
        ascending=[True, False, False, True],
    )
    return gtf_df.drop_duplicates("gene_key", keep="first")


def bed_safe(value):
    return str(value).replace(" ", "_").replace("|", "_").replace("/", "_")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--markers", required=True, help="GTEx/TCGA marker gene TSV")
    ap.add_argument("--gtf", dest="annotation", help="backward-compatible alias for --annotation")
    ap.add_argument("--annotation", help="hg38 GTF/GFF3, optionally gzipped")
    ap.add_argument("--annotation-format", default="auto", choices=["auto", "gtf", "gff3"])
    ap.add_argument("--feature", default="transcript", choices=["transcript", "gene", "mRNA"])
    ap.add_argument("--default-source", default="marker")
    ap.add_argument("--all-transcripts", action="store_true")
    ap.add_argument("--output-bed", required=True)
    ap.add_argument("--output-metadata", required=True)
    args = ap.parse_args()

    annotation = args.annotation
    if annotation is None:
        raise ValueError("provide --annotation or --gtf")
    markers = marker_table(args.markers, args.default_source)
    annotation_df = annotation_tss_table(annotation, args.feature, args.annotation_format, set(markers["gene_key"]))
    if not args.all_transcripts:
        annotation_df = choose_representative_tss(annotation_df)
    merged = markers.merge(annotation_df, on="gene_key", how="left", suffixes=("_marker", "_annotation"))
    merged["mapped"] = merged["chrom"].notna()

    bed_df = merged[merged["mapped"]].copy()
    bed_df["name"] = bed_df.apply(
        lambda r: "|".join([bed_safe(r["source"]), bed_safe(r["marker_set"]), bed_safe(r["gene_symbol_marker"])]),
        axis=1,
    )
    bed_df["tss0"] = pd.to_numeric(bed_df["tss0"], errors="raise").astype(int)
    bed_df["tss1"] = pd.to_numeric(bed_df["tss1"], errors="raise").astype(int)
    if "weight" in bed_df.columns:
        bed_df["score_bed"] = pd.to_numeric(bed_df["weight"], errors="coerce").fillna(0).round().clip(0, 1000).astype(int)
    elif "score" in bed_df.columns:
        bed_df["score_bed"] = pd.to_numeric(bed_df["score"], errors="coerce").fillna(0).round().clip(0, 1000).astype(int)
    else:
        bed_df["score_bed"] = 0

    bed_cols = ["chrom", "tss0", "tss1", "name", "score_bed", "strand"]
    out_bed = Path(args.output_bed)
    out_meta = Path(args.output_metadata)
    out_bed.parent.mkdir(parents=True, exist_ok=True)
    out_meta.parent.mkdir(parents=True, exist_ok=True)
    bed_df[bed_cols].to_csv(out_bed, sep="\t", index=False, header=False)
    merged.to_csv(out_meta, sep="\t", index=False)

    print(f"markers={len(markers)} mapped={int(merged['mapped'].sum())} unmapped={int((~merged['mapped']).sum())}")
    print(f"wrote BED6 -> {out_bed}")
    print(f"wrote metadata -> {out_meta}")
    if (~merged["mapped"]).any():
        print("top unmapped genes:")
        print(merged.loc[~merged["mapped"], ["source", "marker_set", "gene_symbol_marker"]].head(30).to_string(index=False))


if __name__ == "__main__":
    main()
