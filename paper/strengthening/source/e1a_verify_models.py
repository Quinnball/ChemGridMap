"""只读复核 E1-A 模型输出；不调用待验的模型构建函数。"""
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / "paper"), str(ROOT / "src")]

import numpy as np
import pandas as pd
from scipy.stats import norm

from reliability_activity import entry_folds
from chemgridmap.evidence import document_values, document_pair_offsets

S = ROOT / "paper/strengthening"
TARGETS = ("chembl205", "chembl204", "chembl240")
KEY = ["target", "method", "fold", "entry", "held_out_unit"]


def pooled_tau(units):
    moments = units.groupby("entry").document_value.agg(["size", "var"])
    multi = moments[moments["size"] >= 2]
    degrees = multi["size"] - 1
    tau = float(np.sqrt((degrees * multi["var"]).sum() / degrees.sum())) if len(multi) >= 5 else .64
    return tau, multi


def normal_probabilities(values, tau):
    center = values.mean()
    scale = tau * np.sqrt(1 + 1 / len(values))
    low, high = norm.cdf((6 - center) / scale), norm.sf((7 - center) / scale)
    return np.array([low, 1 - low - high, high])


def main():
    pred = pd.read_csv(S / "output/E1A/predictions.csv", float_precision="round_trip")
    fits = pd.read_csv(S / "qa/E1A/fold_fits.csv", float_precision="round_trip")
    frozen = pd.read_csv(ROOT / "paper/output/reliability_v1/activity/lodo_predictions.csv")
    base = pred[pred.population == "all"].copy()
    assert not base.duplicated(KEY).any()
    assert set(base.method) == {"B0", "B1", "G", "G-record", "G-fixed", "X"}
    assert base.groupby("method").size().eq(820).all()
    for method in ("B0", "G", "X"):
        expected = frozen[(frozen.method == method) & (frozen.evaluation_set == "merged")].set_index(KEY).sort_index()
        actual = base[base.method == method].set_index(KEY).sort_index()
        pd.testing.assert_frame_equal(actual[expected.columns], expected, check_exact=True)
        print(f"{method}: 820 rows, all frozen fields exactly preserved after CSV round-trip")

    errors = {"B1": 0., "G-record": 0., "G-fixed": 0.}
    fold_count, record_count, merged_held_count = 0, 0, 0
    for target in TARGETS:
        source = pd.read_csv(ROOT / f"chemgridmap/paper/output/revision_v6/{target}/{target}_retained_records.csv")
        doc = document_values(source)
        rec = document_values(source, document="activity_id", merge_identical=False)
        assignments = entry_folds(doc, seed=20261024)
        pairs = document_pair_offsets(source)
        flagged = pairs[pairs.systematic_offset]
        documents = set(flagged.document_a) | set(flagged.document_b)
        # 从原始记录直接识别整条目剔除集合，独立于合并单元的构成文献实现。
        excluded = set(source.loc[source.document_chembl_id.isin(documents), "molecule_identity_key"])
        selected = base[base.target == target]
        assert np.array_equal(selected.flag_clean, ~selected.entry.isin(excluded))
        assert selected.groupby("entry").flag_clean.nunique().eq(1).all()
        expected_clean = selected[selected.flag_clean].set_index(KEY).sort_index()
        actual_clean = pred[(pred.target == target) & (pred.population == "flag-clean")].set_index(KEY).sort_index()
        pd.testing.assert_frame_equal(actual_clean.drop(columns="population"),
                                      expected_clean.drop(columns="population"), check_exact=True)
        clean_n = len(selected[(selected.method == "G") & selected.flag_clean])
        print(f"{target}: flagged pairs={len(flagged)}, excluded entries={len(excluded)}, clean G predictions={clean_n}")
        b0 = selected[selected.method == "B0"]
        sources = dict(tuple(source.groupby("molecule_identity_key")))
        units = dict(tuple(doc.groupby("entry")))
        for fold in range(10):
            test = {entry for entry, assigned in assignments.items() if assigned == fold}
            train = rec[~rec.entry.isin(test)]
            assert set(train.entry).isdisjoint(test)
            tau, multi = pooled_tau(train)
            fitted = fits[(fits.target == target) & (fits.fold == fold) & (fits.method == "G-record")].iloc[0]
            assert abs(fitted.tau - tau) < 1e-14
            assert fitted.n_entries == len(multi) and fitted.n_units == multi["size"].sum()
            assert fitted.training_entries == train.entry.nunique()
            out_of_fold = b0[~b0.entry.isin(test)]
            constant = float((out_of_fold.observed == out_of_fold.displayed).mean())
            b1 = selected[(selected.method == "B1") & (selected.fold == fold)]
            assert all(assignments[entry] == fold for entry in b1.entry)
            for row in b1.itertuples(index=False):
                expected = np.full(3, (1 - constant) / 2)
                expected[row.displayed] = constant
                error = np.max(np.abs(expected - np.array([row.p0, row.p1, row.p2])))
                errors["B1"] = max(errors["B1"], error)
                assert error < 1e-14
            for method in ("G-record", "G-fixed"):
                for row in selected[(selected.method == method) & (selected.fold == fold)].itertuples(index=False):
                    held = set(row.held_out_unit.split(";"))
                    raw = sources[row.entry]
                    remaining = raw.loc[~raw.document_chembl_id.astype(str).isin(held), "pchembl_value"].to_numpy(float)
                    rest_units = units[row.entry].loc[units[row.entry].documents != row.held_out_unit, "document_value"].to_numpy(float)
                    values, used_tau = (remaining, tau) if method == "G-record" else (rest_units, .64)
                    expected = normal_probabilities(values, used_tau)
                    error = np.max(np.abs(expected - np.array([row.p0, row.p1, row.p2])))
                    errors[method] = max(errors[method], error)
                    assert error < 1e-14 and abs(row.tau - used_tau) < 1e-14
                    assert row.remaining_records == len(remaining)
                    assert row.remaining_units == len(rest_units)
                    assert row.stratum == ("m=1" if len(rest_units) == 1 else "m>=2")
                    if method == "G-record":
                        assert row.used_units == len(remaining)
                        record_count += 1
                        merged_held_count += len(held) > 1
            fold_count += 1
    print(f"B1: {fold_count} folds independently recomputed; max probability error={errors['B1']:.17g}")
    print(f"G-record: {fold_count} fold variances and {record_count} predictions independently recomputed; "
          f"merged held units={merged_held_count}; max probability error={errors['G-record']:.17g}")
    print(f"G-fixed: 820 predictions independently recomputed; max probability error={errors['G-fixed']:.17g}")
    print("MODEL_OUTPUT_AUDIT PASS")


if __name__ == "__main__":
    main()
