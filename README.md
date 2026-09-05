# From Sequence to Structure: Ligand Chemistry and Conformational Ensembles as Predictors of Allosteric Residues

Benchmark codebase for the paper *"Accurate prediction of allosteric sites..."*  
Compares sequence-based (ESM-2), ligand-based (MACCS fingerprints), and structure-based (FPocket/dpocket conformational ensembles) features for residue-level binary classification of allosteric vs. non-allosteric binding-site residues.

---

## Repository layout

```
allosteric-conf-plm/
├── asd_data/               # Allosteric Database snapshots and CV fold splits
│   ├── ASD_2023.tsv        # Raw ASD export
│   ├── ASD_entries.json    # Per-UniProt allosteric site annotations
│   ├── ASD_chains.json     # PDB chain → UniProt mappings
│   ├── ASD_pdbs.json       # PDB accession list
│   ├── ASD_seqs.fasta      # Protein sequences (for MMseqs2 clustering)
│   ├── ASD_clusters_30.json  # 30% sequence-identity clusters
│   ├── train/              # allosteric_ligands_train_fold_{0-4}.json
│   └── val/                # allosteric_ligands_val_fold_{0-4}.json
├── features/               # Pre-extracted .npy arrays (gitignored — see below)
├── results/
│   ├── maccs_cv_results.json       # RF cross-validation scores (MACCS model)
│   ├── maccs_top10_features.json   # Top-10 MACCS feature importances
│   ├── dpocket_cv_results.json     # RF cross-validation scores (dpocket model)
│   └── dpocket_top10_features.json # Top-10 dpocket feature importances
├── metacentrum/            # Cluster job scripts and outputs (gitignored)
└── src/train_scripts/
    ├── config.py                   # All data paths and constants — edit this first
    ├── mmseqs_command.sh           # Sequence clustering command (30% identity)
    ├── data/
    │   ├── prepare_reports.py      # Step 1 – map residues to ligands per fold
    │   ├── extract_features.py     # Step 2 – extract dpocket / MACCS / ESM arrays
    │   ├── asd_data.py             # ASD data accessors
    │   └── io.py                   # Array loading and feature composition
    ├── models/
    │   ├── rf.py                   # Random Forest (primary model)
    │   ├── knn.py                  # k-NN baseline
    │   ├── mlp.py                  # MLP baseline (PyTorch)
    │   └── metrics.py              # Classification metrics helpers
    ├── train_rf.py                 # Step 3 – train RF and save results
    ├── train_knn.py                # Baseline: k-NN across all feature types
    └── train_mlp.py                # Baseline: MLP across all feature types
```

`features/` and `metacentrum/` are gitignored because they contain large binary files and environment-specific cluster scripts respectively.

---

## Prerequisites

External tools required before running the pipeline:

| Tool | Purpose |
|---|---|
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | Sequence clustering for CV fold construction |
| [fpocket / dpocket](https://github.com/Discngine/fpocket) | Pocket detection on PDB conformational ensembles |
| [ESM-2](https://github.com/facebookresearch/esm) | Protein language model embeddings |
| RDKit | MACCS fingerprint computation |

Python dependencies: `scikit-learn`, `numpy`, `pandas`, `rdkit`, `tqdm`, `joblib`, `torch` (for MLP only).

---

## Configuration

All file system paths are centralised in [src/train_scripts/config.py](src/train_scripts/config.py). Update these before running any script:

```python
HOME_FOLDER      = "/storage/praha1/home/hamzagamouh"

PDB_FOLDER       = f"{HOME_FOLDER}/allosteric/pdb_files"       # raw PDB structures
LIGAND_INFO_DIR  = f"{HOME_FOLDER}/allosteric/asd_processing/ligand_info"
MACCS_FPS_PATH   = f"{HOME_FOLDER}/allosteric/asd_processing/ligands_maccs.pkl"
ASD_ENTRIES_PATH = ".../ASD_entries.json"

DPOCKET_OUT_DIR  = ".../dpocket_feats/outputs"   # dpocket output tables
ESM_EMB_DIR      = ".../method_2"                # ASD_<unp>_esm_embs.npy files
FEATURES_DIR     = ".../ligand_analysis"         # fold JSON + .npy arrays land here
```

---

## Pipeline

### Step 0 — Sequence clustering (CV construction)

Cluster ASD sequences at 30% identity to build cross-validation folds that prevent train/val leakage:

```bash
bash src/train_scripts/mmseqs_command.sh
```

This produces `ASD_30_clusterRes` which is used to populate `asd_data/ASD_clusters_30.json` and subsequently the per-fold split files in `asd_data/train/` and `asd_data/val/`.

### Step 1 — Prepare residue–ligand reports

Map every binding-site residue to the ligands that contact it, and split them into allosteric vs. non-allosteric groups for each CV fold:

```bash
python -m train_scripts.data.prepare_reports \
    --features-dir /path/to/FEATURES_DIR \
    --folds 0 1 2 3 4 \
    --modes train val
```

**Inputs** (read from `FEATURES_DIR`):
- `{mode}_unps_fold_{fold}.json` — UniProt IDs for each protein in the fold
- `{mode}_y_fold_{fold}.pkl` — per-residue binary labels

**Outputs**:
- `allosteric_ligands_{mode}_fold_{fold}.json` — per-protein dict with allosteric/non-allosteric residue lists and their associated ligand identifiers

### Step 2 — Extract feature arrays

Convert the residue–ligand reports into numerical feature matrices (`.npy` files). The pre-extracted arrays used in the paper are not tracked in git due to their size; place them in `features/` (the path set as `FEATURES_DIR` in `config.py`).

To re-extract from raw data:

```bash
python -m train_scripts.data.extract_features \
    --output-dir features/ \
    --folds 0 1 2 3 4 \
    --modes train val
```

Optional flags to skip individual feature types:
- `--skip-dpocket-maccs` — skip dpocket and MACCS extraction
- `--skip-esm` — skip ESM embedding extraction

**Outputs written per fold and mode**:

| File | Description |
|---|---|
| `{mode}_{allo\|non_allo}_dpocket_fold_{fold}.npy` | FPocket/dpocket statistics aggregated across conformational ensemble (27 features × 7 stats = 162 columns after NaN filtering) |
| `{mode}_{allo\|non_allo}_maccs_fold_{fold}.npy` | Union MACCS fingerprint (167 bits) of all ligands contacting the residue |
| `{mode}_{allo\|non_allo}_esm_fold_{fold}.npy` | ESM-2 per-residue embedding (1280-dim) |
| `{mode}_{allo\|non_allo}_res_fold_{fold}.json` | Residue identifier list (`<unp>-<aa><pos>`) matching array row order |

**dpocket feature aggregation**: each residue is described by dpocket runs across all PDB conformers. The 27 raw pocket descriptors (`pock_vol`, `flex`, `drug_score`, surface area terms, etc.) are each summarised with 7 statistics (mean, std, min, 25%, 50%, 75%, max), then the NaN-heavy second stat block (columns 27–53) is dropped, leaving 162 valid columns.

**MACCS fingerprint**: when multiple ligands contact the same residue, their 167-bit MACCS keys are OR-combined (max pooling across ligands).

### Step 3 — Train classifiers

All classifiers are run from `src/` with the `features/` directory on the path set in `config.py`. Training data is class-balanced by downsampling the non-allosteric majority class to match allosteric count.

The primary classifier is a Random Forest evaluated with 5-fold cross-validation. The trained model reports six metrics per fold: MCC, F1, accuracy, precision, recall, and ROC-AUC. Results and top-10 feature importances are written to JSON.

```bash
cd src
python -m train_scripts.train_rf --model maccs --n-estimators 100 --output-dir ../results/
```

Available `--model` values:

| Value | Feature vector |
|---|---|
| `maccs` | 167-bit MACCS ligand fingerprint |
| `dpocket` | 162-column dpocket ensemble statistics |
| `esm` | 1280-dim ESM-2 sequence embedding |
| `esm+maccs` | ESM-2 + MACCS (concatenated) |
| `esm+dpocket` | ESM-2 + dpocket (concatenated) |
| `esm+dpocket+maccs` | ESM-2 + dpocket + MACCS (concatenated) |

Two additional baselines are also available: a k-NN classifier (`train_knn.py --n-neighbors 30`) and a two-hidden-layer MLP (`train_mlp.py --num-epochs 1000 --hidden-dim 512`); the MLP saves learning-curve plots per feature type.

---

## Results

### Cross-validation MCC (Random Forest, 100 estimators, 5 folds)

| Model | Val MCC (mean ± std) |
|---|---|
| `maccs` | **0.421 ± 0.026** |
| `dpocket` | 0.332 ± 0.057 |

Full per-fold breakdowns: [results/maccs_cv_results.json](results/maccs_cv_results.json), [results/dpocket_cv_results.json](results/dpocket_cv_results.json).

### Top feature importances

**MACCS model** — allosteric ligands are characterised by aromatic and ring-bearing chemistry:

| Rank | Feature | Importance |
|---|---|---|
| 1 | MACCS-165 RING | 0.0405 |
| 2 | MACCS-162 AROMATIC | 0.0317 |
| 3 | MACCS-125 AROMATIC RING > 1 | 0.0298 |
| 4 | MACCS-96 5M RING | 0.0292 |
| 5 | MACCS-121 N HETEROCYCLE | 0.0291 |

Full list: [results/maccs_top10_features.json](results/maccs_top10_features.json)

**dpocket model** — allosteric sites show distinctive structural variability across conformers:

| Rank | Feature | Importance |
|---|---|---|
| 1 | pock_vol-std | 0.0854 |
| 2 | surf_pol_vdw22-min | 0.0852 |
| 3 | surf_pol_vdw14-min | 0.0504 |
| 4 | n_abpa-std | 0.0372 |
| 5 | flex-50% | 0.0206 |

The dominance of **standard-deviation** and **min/max** statistics over mean-based features confirms that conformational heterogeneity — not mean pocket geometry — is the primary discriminator of allosteric sites.

Full list: [results/dpocket_top10_features.json](results/dpocket_top10_features.json)

---

## Output file format

`{model}_cv_results.json` — one entry per metric (MCC, F1, accuracy, precision, recall, ROC-AUC):
```json
{
  "model": "maccs",
  "folds": [0, 1, 2, 3, 4],
  "n_estimators": 100,
  "train_mcc_per_fold": [...],   "val_mcc_per_fold": [...],
  "train_mcc_mean": 0.922,       "train_mcc_std": 0.013,
  "val_mcc_mean": 0.421,         "val_mcc_std": 0.026,
  "train_f1_per_fold": [...],    "val_f1_per_fold": [...],
  "train_f1_mean": ...,          "val_f1_mean": ...,
  "train_accuracy_mean": ...,    "val_accuracy_mean": ...,
  "train_precision_mean": ...,   "val_precision_mean": ...,
  "train_recall_mean": ...,      "val_recall_mean": ...,
  "train_roc_auc_mean": ...,     "val_roc_auc_mean": ...
}
```

`{model}_top10_features.json`:
```json
[
  {"rank": 1, "feature": "MACCS-165 RING", "importance": 0.0405},
  ...
]
```
