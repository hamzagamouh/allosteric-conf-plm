"""
Rebuild the per-residue MACCS fingerprint feature arrays (`*_maccs_*.npy`),
keeping only the ligands whose 3-letter PDB code belongs to a given set.

By default the set is the "biologically plausible" ligand list
(``asd_data/biologically_plausible_ligands.json``): endogenous metabolites,
cofactors, ions, hormones, neurotransmitters, signalling molecules, lipids,
sugars and microbial/plant natural products - i.e. everything *except*
synthesized medicinal-chemistry compounds / crystallisation additives.

For every binding-site residue the MACCS feature is the element-wise maximum
(bit-wise union) of the fingerprints of the *allowed* ligands seen at that
residue's pockets. Residues whose ligands are all filtered out get an all-zero
vector, so the output arrays stay **row-aligned** with the existing
``*_esm_*.npy`` / ``*_dpocket_*.npy`` / ``*_res_*.json`` files and can be dropped
straight into training (``train_rf.py --maccs-suffix ...``).

Inputs
------
    {json_dir}/{mode}/allosteric_ligands_{mode}_fold_{fold}.json   (asd_data/)
    dpocket outputs (DPOCKET_OUT_DIR)   - only to reproduce residue selection
    MACCS_FPS_PATH                      - {code: np.ndarray(167)}

Outputs (to --output-dir, default FEATURES_DIR)
----------------------------------------------
    {mode}_{allo|non_allo}_maccs_{suffix}_fold_{fold}.npy
    {mode}_{allo|non_allo}_res_{suffix}_fold_{fold}.json
    maccs_{suffix}_manifest.json

Examples
--------
    python -m train_scripts.data.build_maccs_filtered
    python -m train_scripts.data.build_maccs_filtered --ligand-set my_codes.json --suffix myset
    python -m train_scripts.data.build_maccs_filtered --folds 0 --modes train
"""
import argparse
import json
import os
import pickle

import numpy as np

try:
    from joblib import Parallel, delayed
    _HAVE_JOBLIB = True
except ImportError:  # serial fallback
    _HAVE_JOBLIB = False

from train_scripts.config import (
    FEATURES_DIR, MACCS_FPS_PATH, N_FOLDS, DPOCKET_OUT_DIR, DPOCKET_FEATURES,
)


def get_dpocket_feats(pockets):
    """Same behaviour as extract_features.get_dpocket_feats (kept local to avoid a
    hard joblib import). Returns None iff no pocket has a usable dpocket table -
    this is exactly the condition that makes extract_features skip a whole group,
    so residue selection stays identical."""
    import pandas as pd
    dfs, valid = [], []
    for pocket in pockets:
        path = f"{DPOCKET_OUT_DIR}/{pocket}/{pocket}.txt"
        if not os.path.exists(path):
            continue
        try:
            df = pd.read_csv(path, sep=r"\s+")[DPOCKET_FEATURES]
            if len(df) > 0:
                dfs.append(df)
                valid.append(pocket)
        except Exception as e:
            print(f"Problem with pocket {pocket}: {e}")
    if not dfs:
        return None
    result = pd.concat(dfs)
    result["pocket"] = valid
    return result

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
DEFAULT_LIGAND_SET = os.path.join(REPO_ROOT, "asd_data", "biologically_plausible_ligands.json")
DEFAULT_JSON_DIR = os.path.join(REPO_ROOT, "asd_data")

_GROUPS = [("Allosteric", "allo"), ("Non-Allosteric", "non_allo")]

maccs_fps = pickle.load(open(MACCS_FPS_PATH, "rb"))
_FP_DIM = next(iter(maccs_fps.values())).shape[0]


def ligand_code(pocket_id):
    """`1abc-A-ATP-401` -> `ATP` (matches extract_features convention)."""
    return pocket_id.split("-")[2]


def filtered_union_fp(selected_pockets, allowed):
    """Bit-wise union of MACCS fps for allowed ligands; zeros if none allowed."""
    fps = [
        maccs_fps[ligand_code(p)].reshape(1, -1)
        for p in selected_pockets
        if ligand_code(p) in allowed
    ]
    if not fps:
        return np.zeros((1, _FP_DIM), dtype=maccs_fps[next(iter(maccs_fps))].dtype), 0
    return np.concatenate(fps, axis=0).max(axis=0).reshape(1, -1), len(fps)


def _process_entry(results_data, unp, allowed):
    """Reproduce extract_features residue selection, emit filtered MACCS rows."""
    y = results_data[unp]
    out = {"allo": ([], [], 0, 0), "non_allo": ([], [], 0, 0)}

    for group, label in _GROUPS:
        maccs_rows, res_list = [], []
        n_kept_residues, n_zeroed = 0, 0

        try:
            pocket_df = get_dpocket_feats(y[f"{group} ligands"])
        except Exception:
            print(f"Problem getting dpocket feats for {group} pockets of {unp}")
            pocket_df = None
        if pocket_df is None:
            out[label] = (maccs_rows, res_list, n_kept_residues, n_zeroed)
            continue

        for residue in y[f"{group} residues"]:
            for res_name, selected_pockets in residue.items():
                if not selected_pockets:
                    continue
                try:
                    fp, n_allowed = filtered_union_fp(selected_pockets, allowed)
                except Exception:
                    print(f"Problem at residue {res_name} in {unp}")
                    continue
                maccs_rows.append(fp)
                res_list.append(f"{unp}-{res_name}")
                n_kept_residues += 1
                if n_allowed == 0:
                    n_zeroed += 1

        out[label] = (maccs_rows, res_list, n_kept_residues, n_zeroed)

    return out["allo"], out["non_allo"]


def build_for_fold(mode, fold, allowed, json_dir, output_dir, suffix):
    json_path = os.path.join(json_dir, mode, f"allosteric_ligands_{mode}_fold_{fold}.json")
    results_data = json.load(open(json_path))

    if _HAVE_JOBLIB:
        raw = Parallel(n_jobs=-1, backend="multiprocessing", verbose=1)(
            delayed(_process_entry)(results_data, unp, allowed) for unp in results_data
        )
    else:
        raw = [_process_entry(results_data, unp, allowed) for unp in results_data]

    stats = {}
    for label, idx in (("allo", 0), ("non_allo", 1)):
        rows, res, n_res, n_zero = [], [], 0, 0
        for entry in raw:
            m, r, k, z = entry[idx]
            rows.extend(m)
            res.extend(r)
            n_res += k
            n_zero += z

        maccs = np.concatenate(rows, axis=0) if rows else np.zeros((0, _FP_DIM))
        maccs_out = f"{output_dir}/{mode}_{label}_maccs_{suffix}_fold_{fold}.npy"
        res_out = f"{output_dir}/{mode}_{label}_res_{suffix}_fold_{fold}.json"
        np.save(maccs_out, maccs)
        json.dump(res, open(res_out, "w"), indent=1)

        # sanity: alignment with the default (unfiltered) arrays
        default_maccs = f"{output_dir}/{mode}_{label}_maccs_fold_{fold}.npy"
        aligned = None
        if os.path.exists(default_maccs):
            n_default = np.load(default_maccs, mmap_mode="r").shape[0]
            aligned = bool(n_default == maccs.shape[0])
            flag = "OK" if aligned else "!! MISMATCH"
            print(f"[{mode} {label} fold {fold}] rows={maccs.shape[0]} "
                  f"(default={n_default}) {flag}; zeroed={n_zero}/{n_res}")
            if not aligned:
                print("  WARNING: row count differs from the default MACCS array - "
                      "mixed feature sets (esm+maccs, ...) will be misaligned for this "
                      "group. Re-run extract_features so all arrays share one residue set.")
        else:
            print(f"[{mode} {label} fold {fold}] rows={maccs.shape[0]}; "
                  f"zeroed={n_zero}/{n_res} (no default array to compare)")

        stats[label] = {
            "rows": int(maccs.shape[0]),
            "residues_zeroed": int(n_zero),
            "aligned_with_default": aligned,
        }
    return stats


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--ligand-set", default=DEFAULT_LIGAND_SET,
                        help="JSON file with a top-level 'codes' list (or a bare JSON list) "
                             "of allowed 3-letter ligand codes. "
                             f"Default: {DEFAULT_LIGAND_SET}")
    parser.add_argument("--suffix", default=None,
                        help="Filename tag for the outputs. Default: 'bioplausible' for the "
                             "default set, otherwise the ligand-set file stem.")
    parser.add_argument("--json-dir", default=DEFAULT_JSON_DIR,
                        help="Directory holding {mode}/allosteric_ligands_{mode}_fold_{fold}.json")
    parser.add_argument("--output-dir", default=FEATURES_DIR)
    parser.add_argument("--folds", nargs="+", type=int, default=list(range(N_FOLDS)))
    parser.add_argument("--modes", nargs="+", default=["train", "val"])
    args = parser.parse_args()

    spec = json.load(open(args.ligand_set))
    codes = spec["codes"] if isinstance(spec, dict) else spec
    allowed = {c.strip().upper() for c in codes}

    if args.suffix:
        suffix = args.suffix
    elif os.path.abspath(args.ligand_set) == DEFAULT_LIGAND_SET:
        suffix = "bioplausible"
    else:
        suffix = os.path.splitext(os.path.basename(args.ligand_set))[0]

    os.makedirs(args.output_dir, exist_ok=True)
    print(f"Ligand filter: {len(allowed)} codes from {args.ligand_set}")
    print(f"Output suffix: '{suffix}'  ->  {args.output_dir}/{{mode}}_{{allo|non_allo}}_maccs_{suffix}_fold_*.npy\n")

    manifest = {
        "ligand_set": os.path.abspath(args.ligand_set),
        "n_allowed_codes": len(allowed),
        "allowed_codes": sorted(allowed),
        "suffix": suffix,
        "folds": {},
    }
    for fold in args.folds:
        for mode in args.modes:
            print(f"=== {mode.upper()} fold {fold} ===")
            manifest["folds"].setdefault(str(fold), {})[mode] = build_for_fold(
                mode, fold, allowed, args.json_dir, args.output_dir, suffix
            )

    manifest_path = f"{args.output_dir}/maccs_{suffix}_manifest.json"
    json.dump(manifest, open(manifest_path, "w"), indent=2)
    print(f"\nWrote {manifest_path}")


if __name__ == "__main__":
    main()
