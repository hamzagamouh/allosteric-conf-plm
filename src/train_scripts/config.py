import os
from pathlib import Path

# Repository root (…/allosteric-conf-plm)
REPO_ROOT = Path(__file__).resolve().parents[2]

# Root folder holding the raw data. Defaults to <repo>/data; override with the
# ALLOSTERIC_DATA_DIR environment variable or edit this line.
DATA_DIR = os.environ.get("ALLOSTERIC_DATA_DIR", str(REPO_ROOT / "data"))

# Raw data paths
PDB_FOLDER = f"{DATA_DIR}/pdb_files"
LIGAND_INFO_DIR = f"{DATA_DIR}/asd_processing/ligand_info"
MACCS_FPS_PATH = f"{DATA_DIR}/asd_processing/ligands_maccs.pkl"
ASD_ENTRIES_PATH = f"{DATA_DIR}/ASD_entries.json"

# dpocket
DPOCKET_OUT_DIR = f"{DATA_DIR}/dpocket_feats/outputs"
DPOCKET_INP_DIR = f"{DATA_DIR}/dpocket_feats/inputs"
DPOCKET_POCKET_PDBS_DIR = f"{DATA_DIR}/dpocket_feats/pocket_pdbs"
DPOCKET_EXEC = os.environ.get("DPOCKET_EXEC", "dpocket")  # assumes dpocket is on PATH

# ESM embeddings (ASD_<unp>_esm_embs.npy files)
ESM_EMB_DIR = f"{DATA_DIR}/esm_embeddings"

# Fold JSON files and pre-extracted feature arrays live here
# (allosteric_ligands_{mode}_fold_{fold}.json, {mode}_{label}_{feat}_fold_{fold}.npy, etc.)
FEATURES_DIR = os.environ.get("ALLOSTERIC_FEATURES_DIR", str(REPO_ROOT / "features"))

N_FOLDS = 5

DPOCKET_FEATURES = (
    "pock_vol nb_AS nb_AS_norm mean_as_ray mean_as_solv_acc apol_as_prop "
    "apol_as_prop_norm mean_loc_hyd_dens mean_loc_hyd_dens_norm hydrophobicity_score "
    "volume_score polarity_score polarity_score_norm charge_score flex prop_polar_atm "
    "as_density as_density_norm as_max_dst as_max_dst_norm drug_score convex_hull_volume "
    "surf_pol_vdw14 surf_pol_vdw22 surf_apol_vdw14 surf_apol_vdw22 n_abpa"
).split()

# Columns 27-53 (second stat-block) are systematically NaN-heavy; drop them by default.
# Each of 27 dpocket features has 7 stats → 189 total columns.
# valid_cols keeps [0..26] + [54..188]
DPOCKET_VALID_COLS = list(range(27)) + list(range(54, len(DPOCKET_FEATURES) * 7))
