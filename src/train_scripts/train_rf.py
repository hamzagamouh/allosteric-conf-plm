"""Cross-validation training with a Random Forest classifier.

Usage examples:
    python -m train_scripts.train_rf --model maccs
    python -m train_scripts.train_rf --model esm+dpocket --folds 0 1 2 3 4
    python -m train_scripts.train_rf --model dpocket --n-estimators 100 --output-dir results/
    python -m train_scripts.train_rf --model maccs --maccs-suffix bioplausible --output-dir results/

Available models: maccs | dpocket | esm | esm+maccs | esm+dpocket | maccs+dpocket | esm+dpocket+maccs

--maccs-suffix loads ligand-filtered fingerprints built by
`python -m train_scripts.data.build_maccs_filtered`.
"""
import argparse
import json
import os
import numpy as np
from collections import defaultdict

from train_scripts.config import FEATURES_DIR, N_FOLDS
from train_scripts.data.io import (
    load_fold_arrays, clean_dpocket_nans, compose_features, get_feature_names,
)
from train_scripts.models.rf import train_rf
from train_scripts.models.metrics import balance_subsample, build_val_arrays

AVAILABLE_MODELS = ["maccs", "dpocket", "esm", "esm+maccs", "esm+dpocket", "maccs+dpocket", "esm+dpocket+maccs"]
METRIC_NAMES = ["mcc", "f1", "accuracy", "precision", "recall", "roc_auc"]


def run(model, folds, features_dir, n_estimators, output_dir, maccs_suffix=None):
    fold_train_metrics = defaultdict(list)
    fold_val_metrics = defaultdict(list)
    fold_importances = []
    feature_names = None

    if maccs_suffix:
        print(f"Using ligand-filtered MACCS arrays: *_maccs_{maccs_suffix}_fold_*.npy")

    for fold in folds:
        print(f"\n--- fold {fold} [{model}] ---")
        arrays = clean_dpocket_nans(load_fold_arrays(fold, features_dir, maccs_suffix=maccs_suffix))

        if feature_names is None:
            feature_names = get_feature_names(model, arrays)

        train_allo = compose_features(arrays, model, "train", "allo")
        train_non_allo = compose_features(arrays, model, "train", "non_allo")
        val_allo = compose_features(arrays, model, "val", "allo")
        val_non_allo = compose_features(arrays, model, "val", "non_allo")

        train_feats, y_train = balance_subsample(train_allo, train_non_allo)
        val_feats, y_val = build_val_arrays(val_allo, val_non_allo)

        tm, vm, importances = train_rf(train_feats, y_train, val_feats, y_val, n_estimators=n_estimators)

        for metric in METRIC_NAMES:
            fold_train_metrics[metric].append(tm[metric])
            fold_val_metrics[metric].append(vm[metric])

        fold_importances.append(importances)

        print(f"  train: mcc={tm['mcc']:.4f}  f1={tm['f1']:.4f}  acc={tm['accuracy']:.4f}"
              f"  prec={tm['precision']:.4f}  rec={tm['recall']:.4f}  auc={tm['roc_auc']:.4f}")
        print(f"  val:   mcc={vm['mcc']:.4f}  f1={vm['f1']:.4f}  acc={vm['accuracy']:.4f}"
              f"  prec={vm['precision']:.4f}  rec={vm['recall']:.4f}  auc={vm['roc_auc']:.4f}")

    mean_imp = np.mean(fold_importances, axis=0)
    top10_idx = np.argsort(mean_imp)[::-1][:10]
    top10 = [
        {"rank": int(i + 1), "feature": feature_names[idx], "importance": float(mean_imp[idx])}
        for i, idx in enumerate(top10_idx)
    ]

    print(f"\n===== RESULTS [{model}] =====")
    for metric in METRIC_NAMES:
        tr = fold_train_metrics[metric]
        vl = fold_val_metrics[metric]
        print(f"  {metric:9s}  train={np.mean(tr):.3f}±{np.std(tr):.3f}  val={np.mean(vl):.3f}±{np.std(vl):.3f}")

    print("\nTop 10 features:")
    for entry in top10:
        print(f"  {entry['rank']:2d}. {entry['feature']:40s} {entry['importance']:.4f}")

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        slug = model.replace("+", "_")
        if maccs_suffix:
            slug = f"{slug}_maccs_{maccs_suffix}"

        results = {
            "model": model,
            "folds": folds,
            "n_estimators": n_estimators,
            "maccs_suffix": maccs_suffix,
        }
        for metric in METRIC_NAMES:
            tr = fold_train_metrics[metric]
            vl = fold_val_metrics[metric]
            results[f"train_{metric}_per_fold"] = tr
            results[f"val_{metric}_per_fold"] = vl
            results[f"train_{metric}_mean"] = float(np.mean(tr))
            results[f"train_{metric}_std"] = float(np.std(tr))
            results[f"val_{metric}_mean"] = float(np.mean(vl))
            results[f"val_{metric}_std"] = float(np.std(vl))

        json.dump(results, open(f"{output_dir}/{slug}_cv_results.json", "w"), indent=2)
        json.dump(top10, open(f"{output_dir}/{slug}_top10_features.json", "w"), indent=2)
        print(f"\nSaved results to {output_dir}/")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--model", required=True, choices=AVAILABLE_MODELS,
        help="Feature set to train on.",
    )
    parser.add_argument("--folds", nargs="+", type=int, default=list(range(N_FOLDS)))
    parser.add_argument("--n-estimators", type=int, default=10)
    parser.add_argument("--features-dir", default=FEATURES_DIR)
    parser.add_argument("--output-dir", default=None, help="Directory to save results JSON files.")
    parser.add_argument(
        "--maccs-suffix", default=None,
        help="Load ligand-filtered MACCS arrays '*_maccs_<suffix>_fold_*.npy' "
             "(produced by build_maccs_filtered.py, e.g. 'bioplausible') instead of "
             "the default MACCS arrays.",
    )
    args = parser.parse_args()

    run(args.model, args.folds, args.features_dir, args.n_estimators, args.output_dir,
        maccs_suffix=args.maccs_suffix)


if __name__ == "__main__":
    main()
