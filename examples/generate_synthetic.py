"""Generate an entirely artificial two-component example."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def generate(directory):
    destination = Path(directory)
    destination.mkdir(parents=True, exist_ok=True)
    count = 60
    first = np.tile([0.0, 1.0], count // 2)
    second = 1.0 - first
    frame = pd.DataFrame({
        "chr": ["synthetic"] * count,
        "start": np.arange(count) * 10,
        "end": np.arange(count) * 10 + 5,
        "startCpG": np.arange(count),
        "endCpG": np.arange(count) + 1,
        "n_cpgs": np.ones(count, dtype=int),
        "component_a": first,
        "component_b": second,
    })
    frame.to_csv(destination / "atlas.tsv", sep="\t", index=False)
    observed = 0.25 * first + 0.75 * second
    meth = np.rint(observed * 100).astype(int)
    bedgraph = pd.DataFrame({
        "chr": frame["chr"], "start": frame.start, "end": frame.start + 1,
        "pct": observed * 100, "meth": meth, "unmeth": 100 - meth,
    })
    bedgraph.to_csv(destination / "mixture.bedGraph", sep="\t", index=False, header=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True)
    generate(parser.parse_args().output_dir)
