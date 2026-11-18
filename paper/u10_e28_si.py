"""Build the version-0.6 addendum to Additional file 1 (sections S23-S25) from the E28 outputs.

Every number is formatted from the output files; nothing is typed by hand. Usage::

    python paper/u10_e28_si.py --e28 /path/e28_out --cores /path/e28_out/cores \
        --detailed-cores paper/output/u10_cores/summary.json --checks checks.json --output addendum.md
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

SET_LABELS = {"primary": "Primary (50 sampled datasets)", "all52": "Secondary (all 52 evaluation datasets)",
              "replication": "Replication (23 new datasets)", "primary_cutoff2010": "Primary, cutoff 2010 (sensitivity)"}
CATEGORY_LABELS = {"none": "No similar reference compound", "insufficient": "Insufficient", "mixed": "Mixed",
                   "active": "Consistent: active", "medium": "Consistent: medium", "inactive": "Consistent: inactive"}
RULE_LABELS = {"most_frequent": "Most frequent scaffold", "coverage": "Coverage target",
               "coverage_brics": "Coverage target + BRICS", "coverage_recap": "Coverage target + RECAP"}


def pct(x, digits=1):
    return "—" if x is None or pd.isna(x) else f"{100 * x:.{digits}f}%"


def num(x, digits=3):
    return "—" if x is None or pd.isna(x) else f"{x:.{digits}f}".replace("-", "−")


def pval(p):
    if p is None or pd.isna(p):
        return "—"
    if p < 0.001:
        mantissa, exponent = f"{p:.1e}".split("e")
        return f"{mantissa} × 10<sup>{int(exponent)}</sup>".replace("-", "−")
    return f"{p:.3f}"


def table(header, rows):
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return "\n".join(lines)


def section_e28(root):
    report = json.loads((root / "report.json").read_text(encoding="utf-8"))
    out = ["## S23 Neighbor evidence for later compounds (E28)", ""]
    out += ["The sealed plan, its translation and the deviation record are in Additional file 2. "
            "Records were dated by their document year; compounds first reported after the cutoff (2015, or 2010 in the "
            "sensitivity analysis) were inserted into a reference map built from the earlier records (PCA, no refinement, "
            "default curation, class support at 0.8). Datasets entered with at least 100 reference entries and 20 later compounds. "
            "Tests: one-sided Wilcoxon signed-rank tests across datasets, Holm-corrected within each endpoint; tests with fewer "
            "than five eligible datasets were not performed.", ""]
    out += ["### S23.1 Deviations", "",
            "D1 (CHEMBL4235 IC50): two curated stereoisomers share one canonical structure without stereochemistry, and "
            "default map building stops on repeated structures; the reference map of this dataset was built with every entry kept "
            "(`duplicate_policy=\"keep\"`). D2 (CHEMBL1980 IC50): a batch of 1,000 later compounds exceeded the empty cells of the "
            "reference grid; batches were limited to the number of empty cells. Neither change alters the neighbor evidence of "
            "any compound, and neither applies to datasets without these errors.", ""]
    out += ["### S23.2 Inclusion", ""]
    rows = []
    for name, entry in report.items():
        inclusion = pd.read_csv(root / name / "dataset_inclusion.csv")
        excluded = "; ".join(f"{r['dataset']}" for r in entry["excluded"]) or "—"
        rows.append([SET_LABELS.get(name, name), entry["datasets"], entry["included"], excluded,
                     f"{entry['undated_raw_records']:,}",
                     ", ".join(f"{k} {v:,}" for k, v in sorted(entry["insertion"]["later_status_counts"].items()))])
    out += [table(["Set", "Datasets", "Included", "Excluded (fewer than 100 reference entries or 20 later compounds)",
                   "Undated raw records dropped", "Later compounds by insertion status"], rows), ""]
    out += ["`in_reference`: later compound identical in structure to a reference compound (not inserted); "
            "`duplicate_query`: repeated structure among the later compounds of one batch (not inserted).", ""]
    out += ["### S23.3 Predeclared tests", ""]
    rows = []
    for name, entry in report.items():
        for test in entry["tests"]:
            rows.append([SET_LABELS.get(name, name), test["endpoint"], test["test"], test["datasets"], test["positive"],
                         num(test["median_difference"]), pval(test["p_one_sided"]), pval(test.get("p_holm"))])
    out += [table(["Set", "Endpoint", "Test", "Datasets", "Positive", "Median difference", "P (one-sided)", "P (Holm)"], rows), ""]
    out += ["N1: agreement of consistent evidence minus chance expectation (datasets with at least ten consistent compounds). "
            "N2: nearest-reference agreement of consistent minus unsettled (insufficient or mixed) compounds. N3: nearest-reference "
            "agreement of consistent minus `none` compounds. For N2 and N3, at least five compounds in each group.", ""]
    out += ["### S23.4 Categories and agreement", ""]
    rows = []
    for name, entry in report.items():
        for endpoint, e in entry["by_endpoint"].items():
            for category in CATEGORY_LABELS:
                agreement = next((r for r in e["per_class"] if r["neighbor_evidence"] == category), None)
                rows.append([SET_LABELS.get(name, name), endpoint, CATEGORY_LABELS[category], f"{e[f'{category}_n']:,}",
                             pct(e[f"{category}_pooled_fraction"]), pct(e[f"{category}_median_fraction"]),
                             (f"{pct(agreement['agreement'])} ({pct(agreement['ci_low'])}–{pct(agreement['ci_high'])}; "
                              f"{agreement['datasets']} datasets)") if agreement else "—"])
    out += [table(["Set", "Endpoint", "Neighbor evidence", "Later compounds", "Pooled fraction",
                   "Median fraction per dataset", "Agreement with measured class (Wilson 95% CI)"], rows), ""]
    rows = []
    for name, entry in report.items():
        for endpoint, e in entry["by_endpoint"].items():
            md = e["measured_distribution"]
            rows.append([SET_LABELS.get(name, name), endpoint, f"{e['later_compounds']:,}", f"{md['inactive']:,}",
                         f"{md['medium']:,}", f"{md['active']:,}", pct(e["n1_agreement_median"]), pct(e["n1_chance_median"]),
                         f"{pct(e['nearest_consistent_pooled'])} ({e['nearest_consistent_n']:,})",
                         f"{pct(e['nearest_unsettled_pooled'])} ({e['nearest_unsettled_n']:,})",
                         f"{pct(e['nearest_none_pooled'])} ({e['nearest_none_n']:,})"])
    out += [table(["Set", "Endpoint", "Later compounds", "Measured inactive", "Measured medium", "Measured active",
                   "N1 agreement, median per dataset", "N1 chance, median per dataset",
                   "Nearest-reference agreement: consistent (n)", "Unsettled (n)", "None (n)"], rows), ""]
    out += ["Nearest-reference agreement is pooled over compounds; for consistent and unsettled compounds the most similar "
            "reference compound has a Tanimoto similarity of at least 0.6, for `none` compounds below 0.6.", ""]
    out += ["### S23.5 Per-dataset quantities (primary set)", ""]
    metrics = pd.read_csv(root / "primary" / "metrics_by_dataset.csv")
    rows = []
    for _, r in metrics.iterrows():
        rows.append([r.target, r.endpoint, f"{int(r.n):,}", pct(r.n_none / r.n), pct((r.n_insufficient + r.n_mixed) / r.n),
                     int(r.n_consistent), pct(r.get("consistent_agreement")), pct(r.get("consistent_chance")),
                     pct(r.nearest_consistent_agreement), pct(r.nearest_unsettled_agreement), pct(r.nearest_none_agreement)])
    out += [table(["Target", "Endpoint", "Later compounds", "None", "Unsettled", "Consistent (n)", "Agreement",
                   "Chance", "Nearest ref.: consistent", "Unsettled", "None"], rows), ""]
    ins = report["primary"]["insertion"]
    out += ["### S23.6 Insertion time and placement", "",
            f"In the primary set, {ins['batches']} insertion batches placed {ins['placement_distance_n']:,} later compounds. "
            f"Batches of at least 900 compounds took a median of {ins['batch_900_plus_seconds_median']:.1f} s "
            f"(range {ins['batch_900_plus_seconds_range'][0]:.1f}–{ins['batch_900_plus_seconds_range'][1]:.1f} s) to the tables and the "
            "map, with four analyses running in parallel on a four-core Linux container. In these unrefined PCA maps, new compounds "
            f"were placed a median of {ins['placement_distance_median']:.1f} cells (interquartile range "
            f"{ins['placement_distance_q25_q75'][0]:.1f}–{ins['placement_distance_q25_q75'][1]:.1f}) from their most similar reference compound.", ""]
    return "\n".join(out)


def section_cores(cores, detailed):
    summary = json.loads((cores / "summary.json").read_text(encoding="utf-8"))
    out = ["## S24 Core rules with fragment candidates", ""]
    out += ["Post hoc on the 480 groups of the three detailed targets, and as the predeclared descriptive analysis of the E28 plan "
            "(section 6) on the groups with at least two members of the 50 sampled datasets (default curation; for datasets "
            "with more than 1,200 entries, an activity-independent subset of 1,200 drawn with seed 42 after sorting by the structure identity "
            "key; Butina groups at 0.6). Times come from one run per dataset with four datasets in parallel on a four-core "
            "Linux container, and are therefore only indicative.", ""]
    blocks = [("Three detailed targets (480 groups, post hoc)", detailed)]
    blocks += [(f"Sampled datasets, all ({summary['all']['groups']:,} groups in {summary['all']['datasets']} datasets)", summary["all"])]
    blocks += [(f"Sampled datasets, {k} ({v['groups']:,} groups)", v) for k, v in summary["by_endpoint"].items()]
    rows = []
    for label, s in blocks:
        for rule, r in s["by_rule"].items():
            if rule not in RULE_LABELS:
                continue
            rows.append([label, RULE_LABELS[rule], f"{r['target_met']:,}", f"{r['strict_majority']:,}", r["mean_heavy_atoms"],
                         f"{r['benzene_cores']:,}", r["total_seconds"], r["median_ms"]])
    out += [table(["Groups", "Rule", "Target met (≥ 0.8)", "Strict majority", "Mean heavy atoms", "Benzene cores",
                   "Total time (s)", "Median per group (ms)"], rows), ""]
    rows = []
    for label, s in blocks:
        c = s["changes"]
        g = s["generic"]
        rows.append([label, f"{c['coverage_brics_vs_coverage']['changed']:,} / {c['coverage_brics_vs_coverage']['larger']:,}",
                     f"{c['coverage_recap_vs_coverage']['changed']:,} / {c['coverage_recap_vs_coverage']['larger']:,}",
                     f"{c['brics_equals_recap']:,}", f"{c['coverage_recap_vs_coverage_brics']['smaller']:,}",
                     f"{g['generic_majority']:,}", f"{g['typed_exact_majority']:,}", f"{g['generic_groups_more_than_typed']:,}"])
    out += [table(["Groups", "BRICS: changed / larger", "RECAP: changed / larger", "BRICS = RECAP", "RECAP smaller than BRICS",
                   "Generic scaffold exact in a strict majority", "Most frequent typed scaffold exact in a strict majority",
                   "Generic groups more members than typed"], rows), ""]
    rows = []
    for label, s in blocks[1:2]:
        for rule, classes in s["by_size_class"].items():
            if rule not in RULE_LABELS:
                continue
            rows.append([RULE_LABELS[rule]] + [f"{classes[k]['mean_heavy_atoms']} ({classes[k]['groups']:,})" if k in classes else "—"
                                               for k in ("2", "3-4", "5-9", ">=10")])
    out += ["Mean heavy atoms of the core by group size, sampled datasets (number of groups):", "",
            table(["Rule", "2", "3–4", "5–9", "≥ 10"], rows), ""]
    if summary["all"].get("recap_groups_with_skipped_members"):
        out += [f"RECAP skipped members with more than 80 heavy atoms in {summary['all']['recap_groups_with_skipped_members']} groups; "
                "these members still count in coverage.", ""]
    return "\n".join(out)


def section_checks(checks):
    data = json.loads(Path(checks).read_text(encoding="utf-8"))
    out = ["## S25 Software checks for version 0.6", ""]
    out += [f"- **Tests:** {data['tests']}", f"- **Unchanged default outputs:** {data['bytes']}",
            f"- **Browser checks:** {data['browser']}", ""]
    return "\n".join(out)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--e28", type=Path, required=True)
    parser.add_argument("--cores", type=Path, required=True)
    parser.add_argument("--detailed-cores", type=Path, required=True)
    parser.add_argument("--checks", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    detailed = json.loads(args.detailed_cores.read_text(encoding="utf-8"))
    parts = ["# Additional file 1, addendum for version 0.6", "",
             "Sections S23–S25 extend Additional file 1 with the evaluation of the features added in version 0.6. "
             "They are generated by `paper/u10_e28_si.py` from the E28 outputs.", "",
             section_e28(args.e28), section_cores(args.cores, detailed), section_checks(args.checks)]
    args.output.write_text("\n".join(parts).rstrip() + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
