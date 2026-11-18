"""Post hoc sensitivity analysis for parent-structure writing variants."""
from __future__ import annotations

from datetime import datetime
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
S = ROOT / "paper/strengthening"
SOURCE = S / "source"
sys.path[:0] = [str(SOURCE), str(ROOT / "paper"), str(ROOT / "src")]

import numpy as np
import pandas as pd
from rdkit import Chem

from chemgridmap.chembl import curate_chembl_activity_data
from chemgridmap.provenance import environment_versions
from chemgridmap.structures import parent_structure
from e2_analysis import (
    ENDPOINTS,
    METRICS,
    MODELS,
    POPULATIONS,
    STRATA,
    comparator_aucs,
    point_scores,
)
from e2_models import (
    FOLD_SEED,
    build_global_fits,
    build_target_predictions,
    prepare_target,
)
from e2b_analysis import (
    BOOT_N,
    BOOTSTRAP_SEED,
    SUPPORT_THRESHOLD,
    agreement_by_display_status,
    h1_h2_summary,
    hypothesis_tests,
)
from e11_analysis import predict_entry_from_records, validation_tests

MAIN = S / "data/chembl37_api"
VALIDATION = S / "data/chembl37_api_validation"
E2 = S / "output/E2"
E11 = S / "output/E11"
OUT = S / "output/E16"
MODELS_E16 = tuple(name for name in MODELS if name != "B1")

MAIN_NEW = (
    ("CHEMBL4979", "IC50"),
    ("CHEMBL2487", "IC50"),
    ("CHEMBL2820", "Ki"),
)
VALIDATION_NEW = (
    ("CHEMBL3884", "IC50"),
    ("CHEMBL1913", "IC50"),
    ("CHEMBL3759", "Ki"),
    ("CHEMBL255", "Ki"),
    ("CHEMBL3105", "Ki"),
)
AUDIT_ENDPOINTS = tuple(
    ("chembl37_api", target, endpoint)
    for target, endpoint in MAIN_NEW
) + tuple(
    ("chembl37_api_validation", target, endpoint)
    for target, endpoint in VALIDATION_NEW
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_safe(value):
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(json_safe(value), indent=2, ensure_ascii=False, allow_nan=False)
        + "\n",
        encoding="utf-8",
    )


def activity_path(dataset: str, target: str, endpoint: str) -> Path:
    base = MAIN if dataset == "chembl37_api" else VALIDATION
    return base / "activities" / f"{target}_{endpoint}.csv"


@lru_cache(maxsize=50000)
def identity_details(smiles: str):
    parent_smiles, _ = parent_structure(smiles)
    parent = Chem.RemoveHs(Chem.MolFromSmiles(parent_smiles))
    normalized = Chem.MolToSmiles(
        parent, canonical=True, isomericSmiles=True
    )
    standard_inchi = Chem.MolToInchi(parent)
    standard_inchikey = Chem.MolToInchiKey(parent)
    charges = tuple(sorted(
        atom.GetFormalCharge()
        for atom in parent.GetAtoms()
        if atom.GetFormalCharge()
    ))
    return (
        parent_smiles,
        normalized,
        standard_inchi,
        standard_inchikey,
        charges,
    )


def classify_variant_kind(parent_smiles, normalized, charges):
    """Classify the reviewed writing variants without changing the comparison."""
    if "[H]" in parent_smiles and parent_smiles != normalized:
        return "explicit_H"
    if len({item for item in charges}) or charges:
        return "resonance"
    return "tautomer"


def audit_endpoint(dataset: str, target: str, endpoint: str):
    raw = pd.read_csv(
        activity_path(dataset, target, endpoint), low_memory=False
    )
    rows = []
    grouped = raw.dropna(
        subset=["parent_molecule_chembl_id", "canonical_smiles"]
    ).groupby("parent_molecule_chembl_id", sort=True)
    for parent_id, group in grouped:
        writings = group.canonical_smiles.astype(str).value_counts()
        if len(writings) <= 1:
            continue
        unique_keys = {
            identity_details(smiles)[3] for smiles in writings.index
        }
        if len(unique_keys) != 1:
            continue
        details = [identity_details(smiles) for smiles in writings.index]
        kinds = {
            classify_variant_kind(item[0], item[1], item[4])
            for item in details
        }
        kind = "explicit_H" if "explicit_H" in kinds else (
            "resonance" if "resonance" in kinds else "tautomer"
        )
        for smiles, count in writings.items():
            parent_smiles, normalized, inchi, inchikey, _ = identity_details(smiles)
            rows.append({
                "dataset": dataset,
                "target": target,
                "endpoint": endpoint,
                "parent_molecule_chembl_id": str(parent_id),
                "source_smiles": smiles,
                "parent_smiles": parent_smiles,
                "remove_hs_smiles": normalized,
                "standard_inchi": inchi,
                "standard_inchikey": inchikey,
                "record_count": int(count),
                "variant_kind": kind,
            })
    return rows


def audit_all_downloads():
    files = sorted([*MAIN.glob("activities/*.csv"), *VALIDATION.glob("activities/*.csv")])
    parent_groups = 0
    variant_groups = 0
    conflicts = []
    for index, path in enumerate(files, 1):
        raw = pd.read_csv(
            path,
            usecols=["parent_molecule_chembl_id", "canonical_smiles"],
            low_memory=False,
        ).dropna(subset=["parent_molecule_chembl_id", "canonical_smiles"])
        for parent_id, group in raw.groupby(
            "parent_molecule_chembl_id", sort=True
        ):
            parent_groups += 1
            writings = group.canonical_smiles.astype(str).unique()
            if len(writings) <= 1:
                continue
            keys = {identity_details(smiles)[3] for smiles in writings}
            if len(keys) == 1:
                variant_groups += 1
            else:
                conflicts.append({
                    "file": path.relative_to(S).as_posix(),
                    "parent_molecule_chembl_id": str(parent_id),
                    "standard_inchikeys": ";".join(sorted(keys)),
                    "record_count": int(len(group)),
                })
        print(f"audit {index}/{len(files)} {path.name}", flush=True)
    return files, {
        "activity_files": len(files),
        "parent_groups": parent_groups,
        "same_inchikey_variant_groups": variant_groups,
        "different_inchikey_conflict_groups": len(conflicts),
    }, conflicts


def curate_dataset(dataset, target, endpoint):
    base = MAIN if dataset == "chembl37_api" else VALIDATION
    raw = pd.read_csv(
        base / "activities" / f"{target}_{endpoint}.csv", low_memory=False
    )
    result = curate_chembl_activity_data(
        raw, target_id=target, activity_type=endpoint
    )
    context = prepare_target(result.retained_records)
    entries = context["doc"].groupby("entry").size()
    multi = int((entries >= 2).sum())
    return raw, result, context, {
        "dataset": dataset,
        "target": target,
        "endpoint": endpoint,
        "records": int(len(raw)),
        "retained_records": int(len(result.retained_records)),
        "entries": int(len(entries)),
        "multi_unit_entries": multi,
        "included": bool(multi >= 30),
        "identity_variant_parents": int(
            result.structure_variants.molecule_identity_key.nunique()
            if len(result.structure_variants) else 0
        ),
    }


def metrics_for_predictions(predictions):
    rows = []
    for (target, endpoint), frame in predictions.groupby(
        ["target", "endpoint"], sort=True
    ):
        for population in POPULATIONS:
            population_frame = frame.loc[frame.population == population]
            for model in MODELS_E16:
                model_frame = population_frame.loc[
                    population_frame.method == model
                ]
                for stratum in STRATA:
                    selected = (
                        model_frame
                        if stratum == "all"
                        else model_frame.loc[model_frame.stratum == stratum]
                    )
                    scores = point_scores(selected)
                    rows.append({
                        "target": target,
                        "endpoint": endpoint,
                        "population": population,
                        "model": model,
                        "stratum": stratum,
                        "n_predictions": int(len(selected)),
                        "n_entries": int(selected.entry.nunique()),
                        **scores,
                    })
    return pd.DataFrame(rows)


def pooled_for_predictions(predictions):
    from e1a_statistics import cluster_weights, percentile_ci, weighted_scores

    rows = []
    for endpoint in ENDPOINTS:
        endpoint_frame = predictions.loc[predictions.endpoint == endpoint]
        for population in POPULATIONS:
            population_frame = endpoint_frame.loc[
                endpoint_frame.population == population
            ]
            for model in MODELS_E16:
                model_frame = population_frame.loc[
                    population_frame.method == model
                ]
                for stratum in STRATA:
                    selected = (
                        model_frame
                        if stratum == "all"
                        else model_frame.loc[model_frame.stratum == stratum]
                    )
                    point = point_scores(selected)
                    weights = cluster_weights(
                        selected,
                        np.random.default_rng(BOOTSTRAP_SEED),
                        2000,
                    )
                    bootstrap = weighted_scores(selected, weights)
                    for metric in METRICS:
                        low, high, valid = percentile_ci(bootstrap[metric])
                        rows.append({
                            "endpoint": endpoint,
                            "population": population,
                            "model": model,
                            "stratum": stratum,
                            "metric": metric,
                            "estimate": point[metric],
                            "ci_low": low,
                            "ci_high": high,
                            "n_bootstrap_valid": int(valid),
                        })
    return pd.DataFrame(rows)


def main_combined(original_inclusion, curated_new):
    keys = [
        (str(row.target), str(row.endpoint))
        for row in original_inclusion.loc[
            original_inclusion.included
        ].itertuples(index=False)
    ]
    new_keys = {
        key for key, item in curated_new.items() if item["included"]
    }
    if not new_keys:
        empty = pd.DataFrame()
        return {
            "metrics": pd.read_csv(E2 / "metrics_by_target.csv"),
            "pooled": pd.read_csv(E2 / "pooled_metrics.csv"),
            "comparators": pd.read_csv(E2 / "comparators.csv"),
            "h1": pd.read_csv(E2 / "h1_h2_summary.csv"),
            "tests": pd.read_csv(E2 / "hypothesis_tests.csv"),
            "support": pd.read_csv(E2 / "agreement_by_display_status.csv"),
            "new_support": empty.copy(),
            "new_dataset_metrics": empty.copy(),
            "records": {key: None for key in keys},
        }
    records = {}
    contexts = {}
    for index, (target, endpoint) in enumerate(keys, 1):
        _, result, context, _ = curate_dataset(
            "chembl37_api", target, endpoint
        )
        records[(target, endpoint)] = result.retained_records
        contexts[(target, endpoint)] = context
        print(f"main baseline {index}/{len(keys)} {target} {endpoint}", flush=True)
    for key, item in curated_new.items():
        if item["included"]:
            records[key] = item["result"].retained_records
            contexts[key] = item["context"]
    global_fits_by_endpoint = {}
    for endpoint in ENDPOINTS:
        endpoint_keys = [key for key in records if key[1] == endpoint]
        global_fits_by_endpoint[endpoint] = build_global_fits(
            {
                target: contexts[(target, endpoint)]["doc"]
                for target, _ in endpoint_keys
            },
            {
                target: contexts[(target, endpoint)]["clean_entries"]
                for target, _ in endpoint_keys
            },
        )
    frames = []
    for index, (target, endpoint) in enumerate(records, 1):
        frame, _ = build_target_predictions(
            records[(target, endpoint)],
            target,
            endpoint,
            global_fits=global_fits_by_endpoint[endpoint],
            models=MODELS_E16,
            context=contexts[(target, endpoint)],
        )
        frames.append(frame)
        print(f"main models {index}/{len(records)} {target} {endpoint}", flush=True)
    predictions = pd.concat(frames, ignore_index=True)
    metrics = metrics_for_predictions(predictions)
    pooled = pooled_for_predictions(predictions)
    comparators = comparator_aucs(predictions)
    h1 = h1_h2_summary(metrics, pooled)
    tests = hypothesis_tests(metrics, comparators, records)
    support = agreement_by_display_status(
        predictions,
        records,
    )
    new_predictions = predictions.loc[
        predictions.set_index(["target", "endpoint"]).index.isin(new_keys)
    ]
    new_support = agreement_by_display_status(
        new_predictions,
        {key: records[key] for key in new_keys},
    )
    new_metrics = metrics.loc[
        metrics.set_index(["target", "endpoint"]).index.isin(new_keys)
    ].copy()
    new_h1 = new_metrics.loc[
        (new_metrics.stratum == "all") & (new_metrics.model == "B0"),
        ["target", "endpoint", "population", "accuracy"],
    ].rename(columns={"accuracy": "b0_accuracy"})
    new_h1["h1"] = 1.0 - new_h1.pop("b0_accuracy")
    new_auc = new_metrics.loc[
        (new_metrics.stratum == "all") & (new_metrics.model == "G"),
        ["target", "endpoint", "population", "auc"],
    ].rename(columns={"auc": "g_auc"})
    new_auc["h3"] = new_auc.g_auc - 0.5
    new_compare = comparators.loc[
        comparators.set_index(["target", "endpoint"]).index.isin(new_keys)
    ].pivot(
        index=["target", "endpoint", "population"],
        columns="model",
        values="auc",
    ).reset_index()
    new_compare["h8"] = new_compare["G"] - new_compare["B-range"]
    new_dataset_metrics = (
        new_h1.merge(new_auc, on=["target", "endpoint", "population"])
        .merge(new_compare, on=["target", "endpoint", "population"])
    )
    return {
        "predictions": predictions,
        "metrics": metrics,
        "pooled": pooled,
        "comparators": comparators,
        "h1": h1,
        "tests": tests,
        "support": support,
        "new_support": new_support,
        "new_dataset_metrics": new_dataset_metrics,
        "records": records,
    }


def combine_report_rows(original, updated, hypotheses):
    original = original.loc[original.hypothesis.isin(hypotheses)].copy()
    updated = updated.loc[updated.hypothesis.isin(hypotheses)].copy()
    original["scenario"] = "original"
    updated["scenario"] = "with_reviewed_variants"
    columns = [
        "scenario", "hypothesis", "endpoint", "population", "n_targets",
        "median_difference", "median_ci_low", "median_ci_high",
        "wilcoxon_statistic", "raw_p", "holm_p",
    ]
    return pd.concat(
        [original[columns], updated[columns]], ignore_index=True
    )


def combine_h1_rows(original, updated):
    original = original.loc[original.hypothesis == "H1"].copy()
    updated = updated.loc[updated.hypothesis == "H1"].copy()
    original["scenario"] = "original"
    updated["scenario"] = "with_reviewed_variants"
    columns = [
        "scenario", "hypothesis", "endpoint", "population", "n_targets",
        "median", "q1", "q3", "pooled_estimate", "pooled_ci_low",
        "pooled_ci_high",
    ]
    return pd.concat(
        [original[columns], updated[columns]], ignore_index=True
    )


def validation_combined(curated_new):
    original = pd.read_csv(E11 / "predictions_e11.csv", low_memory=False)
    frames = [original]
    for index, (key, item) in enumerate(curated_new.items(), 1):
        if not item["included"]:
            continue
        frames.append(predict_entry_from_records(
            item["result"].retained_records,
            key[0],
            key[1],
            context=item["context"],
        ))
        print(f"validation models {index}/{len(curated_new)} {key[0]} {key[1]}", flush=True)
    predictions = pd.concat(frames, ignore_index=True)
    original_out = validation_tests.__globals__["OUT"]
    validation_tests.__globals__["OUT"] = OUT
    try:
        tests = validation_tests(predictions)
    finally:
        validation_tests.__globals__["OUT"] = original_out
    baseline = pd.read_csv(E11 / "validation_tests.csv")
    baseline["scenario"] = "original"
    tests["scenario"] = "with_reviewed_variants"
    columns = [
        "scenario", "endpoint", "test", "n_datasets", "median_difference",
        "wilcoxon_statistic", "p_value", "holm_adjusted_p",
        "holm_significant",
    ]
    return pd.concat(
        [baseline[columns], tests[columns]], ignore_index=True
    ), tests


def write_manifest(started, ended, inputs):
    output_hashes = {
        path.name: sha256(path)
        for path in sorted(OUT.iterdir())
        if path.is_file() and path.name != "manifest.json"
    }
    manifest = {
        "stage": "E16",
        "analysis_role": "post hoc sensitivity",
        "started_at": started,
        "ended_at": ended,
        "parameters": {
            "fold_seed": FOLD_SEED,
            "bootstrap_seed": BOOTSTRAP_SEED,
            "bootstrap_repetitions": 2000,
            "inclusion_rule": ">=30 entries with >=2 evidence units",
            "reference_identity": "RDKit standard InChIKey after RemoveHs",
        },
        "software_versions": environment_versions(),
        "inputs_sha256": {
            path.relative_to(ROOT).as_posix(): sha256(path)
            for path in inputs
        },
        "source_sha256": {
            "source/e16_identity_sensitivity.py": sha256(Path(__file__)),
            "source/test_e16.py": sha256(SOURCE / "test_e16.py"),
        },
        "outputs_sha256": output_hashes,
    }
    write_json(OUT / "manifest.json", manifest)


def run():
    started = datetime.now().astimezone().isoformat()
    OUT.mkdir(parents=True, exist_ok=True)
    audit_rows = []
    for item in AUDIT_ENDPOINTS:
        audit_rows.extend(audit_endpoint(*item))
    pd.DataFrame(audit_rows).to_csv(
        OUT / "reviewed_variants.csv", index=False
    )

    files, audit_summary, conflicts = audit_all_downloads()
    pd.DataFrame(conflicts).to_csv(
        OUT / "all_parent_conflicts.csv", index=False
    )
    write_json(OUT / "audit_summary.json", audit_summary)

    curated = {}
    for dataset, target, endpoint in AUDIT_ENDPOINTS:
        raw, result, context, inclusion = curate_dataset(
            dataset, target, endpoint
        )
        del raw
        curated[(target, endpoint)] = {
            "result": result,
            "context": context,
            **inclusion,
        }
    inclusion = pd.DataFrame([
        {key: value for key, value in item.items()
         if key not in {"result", "context"}}
        for item in curated.values()
    ])
    inclusion.to_csv(OUT / "identity_inclusion.csv", index=False)

    main_new = {
        key: item for key, item in curated.items() if key in MAIN_NEW
    }
    main = main_combined(
        pd.read_csv(E2 / "target_inclusion.csv"),
        main_new,
    )
    original_h1 = pd.read_csv(E2 / "h1_h2_summary.csv")
    original_tests = pd.read_csv(E2 / "hypothesis_tests.csv")
    h1_compare = combine_h1_rows(original_h1, main["h1"])
    test_compare = combine_report_rows(
        original_tests, main["tests"], ("H3", "H8")
    )
    pd.concat([h1_compare, test_compare], ignore_index=True).to_csv(
        OUT / "main_combined_results.csv", index=False
    )
    main["new_dataset_metrics"].to_csv(
        OUT / "main_new_dataset_metrics.csv", index=False
    )
    main["new_support"].to_csv(
        OUT / "main_new_support.csv", index=False
    )
    pd.DataFrame({
        "target": [key[0] for key in main["records"]],
        "endpoint": [key[1] for key in main["records"]],
    }).to_csv(OUT / "main_combined_datasets.csv", index=False)

    validation_new = {
        key: item for key, item in curated.items() if key in VALIDATION_NEW
    }
    validation_compare, _ = validation_combined(validation_new)
    validation_compare.to_csv(
        OUT / "validation_combined_tests.csv", index=False
    )

    inputs = files + [
        E2 / "target_inclusion.csv",
        E2 / "h1_h2_summary.csv",
        E2 / "hypothesis_tests.csv",
        E11 / "target_inclusion.csv",
        E11 / "predictions_e11.csv",
        E11 / "validation_tests.csv",
    ]
    ended = datetime.now().astimezone().isoformat()
    write_manifest(started, ended, inputs)
    print(
        json.dumps({
            "audit": audit_summary,
            "included": inclusion[["dataset", "target", "endpoint", "included"]]
            .to_dict("records"),
            "main_combined_datasets": len(main["records"]),
            "validation_combined_datasets": len(
                pd.read_csv(E11 / "target_inclusion.csv").query("included")
            ) + int(
                inclusion.set_index(["target", "endpoint"]).loc[
                    list(VALIDATION_NEW), "included"
                ].sum()
            ),
        }, ensure_ascii=False),
        flush=True,
    )


if __name__ == "__main__":
    run()
