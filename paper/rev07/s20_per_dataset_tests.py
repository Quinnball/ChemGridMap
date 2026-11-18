"""S20 (post hoc): Wilcoxon signed-rank test of the per-dataset class-standardised differences.

Reads the s20_mh.json written by s20_mh.py and adds, for each design and endpoint, the number of
datasets, how many differences are positive or zero, their median and a two-sided Wilcoxon P value
(zeros dropped; at least five non-zero differences).
Usage: python paper/rev07/s20_per_dataset_tests.py <folder with s20_mh.json>
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats


def main(folder):
    path = Path(folder) / "s20_mh.json"
    rows = json.loads(path.read_text())
    out = {}
    for row in rows:
        key = f"{row['design']}_{row['endpoint']}"
        x = np.array(row["per_dataset"], dtype=float)
        nz = x[x != 0]
        p = float(stats.wilcoxon(nz, alternative="two-sided").pvalue) if len(nz) >= 5 else None
        out[key] = dict(row, pd_n=int(len(x)), pd_positive=int((x > 0).sum()), pd_zero=int((x == 0).sum()),
                        pd_median=float(np.median(x)) if len(x) else None, wilcoxon_two_sided=p)
        print(key, row["mh_rd"], (row["boot_lo"], row["boot_hi"]), "crude", row["crude"], "| per-dataset n", len(x),
              "positive", int((x > 0).sum()), "P", p)
    (Path(folder) / "s20_tests.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main(sys.argv[1])
