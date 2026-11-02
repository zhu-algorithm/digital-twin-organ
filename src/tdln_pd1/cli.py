from __future__ import annotations

import argparse
import json
from pathlib import Path
from .data import load_cells, load_labels
from .features import aggregate_cells
from .model import train_evaluate


def main() -> None:
    parser = argparse.ArgumentParser(description="TDLN PD-1 response baseline")
    sub = parser.add_subparsers(dest="command", required=True)
    featurize = sub.add_parser("featurize")
    featurize.add_argument("--cells", required=True)
    featurize.add_argument("--output", required=True)
    train = sub.add_parser("train")
    train.add_argument("--features", required=True)
    train.add_argument("--labels", required=True)
    train.add_argument("--output", required=True)
    train.add_argument("--folds", type=int, default=5)
    args = parser.parse_args()
    if args.command == "featurize":
        aggregate_cells(load_cells(args.cells)).to_csv(args.output, index=False)
    else:
        import pandas as pd
        metrics = train_evaluate(pd.read_csv(args.features), load_labels(args.labels), Path(args.output), args.folds)
        print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
