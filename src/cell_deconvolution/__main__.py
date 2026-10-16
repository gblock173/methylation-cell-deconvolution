"""Single-input command line interface."""
import argparse
import json
import sys
from pathlib import Path

from .v1 import V1Model


def main():
    parser = argparse.ArgumentParser(description="Coverage-weighted NNLS deconvolution")
    parser.add_argument("--atlas", required=True)
    parser.add_argument("--bedgraph", required=True)
    parser.add_argument("--output", help="JSON output; defaults to standard output")
    parser.add_argument("--config", help="Optional JSON with metadata columns and cell-type order")
    parser.add_argument("--min-cov", type=int, default=None)
    args = parser.parse_args()
    try:
        config = json.loads(Path(args.config).read_text()) if args.config else {}
        if args.min_cov is not None:
            config["min_cov"] = args.min_cov
        result = V1Model(args.atlas, config=config).deconv_bedgraph(args.bedgraph)
        if result is None:
            print("Insufficient coverage: at least 50 usable blocks are required.", file=sys.stderr)
            return 2
        payload = json.dumps(result, indent=2, allow_nan=False) + "\n"
        if args.output:
            Path(args.output).write_text(payload, encoding="utf-8")
        else:
            print(payload, end="")
        return 0
    except (OSError, ValueError, KeyError, IndexError, RuntimeError):
        print("Cannot process inputs. Check file access, schema, ordering and numeric values.", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
