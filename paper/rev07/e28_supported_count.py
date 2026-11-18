"""Agreement of consistent neighbor evidence by the number of supported similar reference compounds.

Answers reviewer question R2 Q3: how many consistent labels rest on one supported reference compound,
and how agreement and chance agreement change with the number of supported reference compounds
(Tanimoto >= 0.6 among the five most similar, ties kept). Post hoc and exploratory; two-sided tests.
Input: compounds_annotated_<sample>.csv.gz written by e28_posthoc.py (in E28_OUT).
Usage: E28_OUT=<folder> python paper/rev07/e28_supported_count.py [primary|replication]
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from e28_posthoc import n1_dataset, n1_test  # noqa: E402

OUT = Path(os.environ.get("E28_OUT", "e28_posthoc_output"))
LEVELS = (("1", lambda n: n == 1), ("2", lambda n: n == 2), ("3", lambda n: n == 3),
          ("4 or more", lambda n: n >= 4), ("2 or more", lambda n: n >= 2))


def main(sample="primary"):
    comp = pd.read_csv(OUT / f"compounds_annotated_{sample}.csv.gz", low_memory=False)
    cons = comp[comp.group.eq("consistent")]
    rows = []
    for ep, e in cons.groupby("endpoint"):
        for label, keep in LEVELS:
            t = e[e.supported_reference_count.map(keep)]
            parts = [n1_dataset(g) for _, g in t.groupby("dataset")]
            parts = [(n, o, x) for d, n, o, x in parts if n]
            chance = sum(n * x for n, o, x in parts) / sum(n for n, o, x in parts)
            test = n1_test(t, min_n=10)
            rows.append(dict(endpoint=ep, supported_reference_compounds=label, compounds=int(len(t)),
                             share=len(t) / len(e), agreement=float(t.measured_class.eq(t.ev).mean()),
                             pooled_chance=chance, **{f"n1_{k}": v for k, v in test.items()}))
    table = pd.DataFrame(rows)
    table.to_csv(OUT / f"q3_supported_count_{sample}.csv", index=False)
    print(table.round(3).to_string(index=False))


if __name__ == "__main__":
    main(*sys.argv[1:2])
