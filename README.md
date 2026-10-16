# Methylation Cell Deconvolution

Author: **Yujin Jeon**  
License: **MIT**

Coverage-weighted non-negative least squares (NNLS) for methylation cell-type
fractions, reference-atlas construction, marker selection, and WPS/TSS tools.

This code release includes an entirely synthetic example generator. Biological
inputs, fitted atlases, sample metadata, cohort results, figures, credentials
and the original execution history are excluded.

## Install and run

~~bash
python -m pip install .
python examples/generate_synthetic.py --output-dir demo-output
python -m cell_deconvolution --atlas demo-output/atlas.tsv --bedgraph demo-output/mixture.bedGraph --output demo-output/fractions.json
~~

The artificial example has expected fractions 0.25 and 0.75. All its coordinates
and counts are generated mathematically.

~~python
from cell_deconvolution import V1Model
model = V1Model("demo-output/atlas.tsv")
result = model.deconv_bedgraph("demo-output/mixture.bedGraph")
~~

Use --config with a private JSON to override meta_cols, blood_types, or to
pin celltypes order. Use --min-cov to override the default of 3.

## Input and numeric conventions

The atlas is TSV: chr, start, end, startCpG, endCpG, n_cpgs, then numeric
cell-type columns. BedGraph has six tab-separated columns: chromosome, start,
end, percentage, methylated count, unmethylated count. Percentage is ignored;
a track header is optional.

Supply chromosome-grouped, coordinate-sorted bedGraph and non-overlapping atlas
blocks. The baseline single-pointer reader assigns each CpG by start coordinate
to [start, end). Overlapping blocks and unsorted inputs are unsupported.

Block AMF = sum(methylated counts) / sum(total counts). NNLS uses sqrt(coverage)
weights, coverage >= 3, finite reference/sample values, and at least 50 blocks.
Coefficients are normalized by their sum. raw_sum is the unnormalized sum.
blood_total sums configured blood cell types.

R2 retains the baseline definition on the square-root-coverage-transformed
response, centered by that transformed response's mean. Extraction, fit,
normalization and reporting preserve the baseline numeric implementation.

## Atlas construction and selection

~~bash
python -m cell_deconvolution.atlas --help
python -m cell_deconvolution.selection --help
~~

The builder takes a reference TSV (cell_type, beta_file), a coverage TSV
(beta_file), and a CpG index. None are distributed. Beta input is uint8
methylated/total pairs, one pair per indexed CpG. File references in manifests
resolve from the current working directory.

The builder retains baseline range, contiguous-run, coverage and pooled-count
calculations. Its legacy startCpG/endCpG values are used directly as Python
slice bounds, distinct from genomic start/end. This is compatibility behavior,
not a corrected indexing convention or newly validated atlas. Use a compatible,
independently checked reference. Whole-genome atlas building needs substantial
memory. Selection uses absolute one-versus-rest z-scores and up to K owned
blocks per cell type.

## Optional fragmentomics

~~bash
python -m pip install ".[fragmentomics]"
python scripts/snyder_wps.py --help
python scripts/compute_wps_tss_features.py --help
python scripts/compute_wps_tss_atomic_features.py --help
python scripts/build_marker_gene_table.py --help
python scripts/build_marker_tss_bed.py --help
~~

Run optional scripts from a checkout. Their BAM, expression, marker and
annotation inputs are user-supplied. Generated outputs can inherit input labels
and file references; keep them outside the release.

## Verification and privacy

~~bash
python -m unittest discover -s tests -v
python tools/privacy_check.py
git config core.hooksPath .githooks
~~

Only .release-files allowlisted files may be committed; .gitignore denies new
files by default. The checker inspects the Git index, reachable historical trees,
and commit/tag metadata. Commit/push hooks and CI run the checks.

Keep all data, reports, logs and runtime configuration outside Git. New allowlist
entries need review. Automated matching cannot prove arbitrary future content
is non-identifying. See RELEASE_AUDIT.md for this release's checks.
