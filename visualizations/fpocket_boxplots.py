"""Boxplots comparing top-10 dpocket features between allosteric and non-allosteric residues.

Uses fold-0 train + val arrays (each residue appears exactly once).
"""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJDIR = "/auto/vestec1-elixir/home/hamzagamouh/allosteric-conf-plm"
FEATURES_DIR = f"{PROJDIR}/features"
TOP10_PATH = f"{PROJDIR}/results/dpocket_top10_features.json"
OUT_DIR = f"{PROJDIR}/visualizations"

DPOCKET_FEATURES = (
    "pock_vol nb_AS nb_AS_norm mean_as_ray mean_as_solv_acc apol_as_prop "
    "apol_as_prop_norm mean_loc_hyd_dens mean_loc_hyd_dens_norm hydrophobicity_score "
    "volume_score polarity_score polarity_score_norm charge_score flex prop_polar_atm "
    "as_density as_density_norm as_max_dst as_max_dst_norm drug_score convex_hull_volume "
    "surf_pol_vdw14 surf_pol_vdw22 surf_apol_vdw14 surf_apol_vdw22 n_abpa"
).split()
STATS = ["mean", "std", "min", "25%", "50%", "75%", "max"]

all_dpocket_names = [f"{feat}-{stat}" for feat in DPOCKET_FEATURES for stat in STATS]
VALID_COLS = list(range(27)) + list(range(54, len(DPOCKET_FEATURES) * 7))
DPOCKET_FEATURE_NAMES = [all_dpocket_names[i] for i in VALID_COLS]

top10 = json.load(open(TOP10_PATH))
top10_names = [e["feature"] for e in top10]
top10_importances = [e["importance"] for e in top10]

col_indices = [DPOCKET_FEATURE_NAMES.index(name) for name in top10_names]

# Fold 0 train + val gives each residue exactly once
allo_parts, non_parts = [], []
for mode in ("train", "val"):
    raw_allo = np.load(f"{FEATURES_DIR}/{mode}_allo_dpocket_fold_0.npy")[:, VALID_COLS]
    raw_non = np.load(f"{FEATURES_DIR}/{mode}_non_allo_dpocket_fold_0.npy")[:, VALID_COLS]
    allo_parts.append(raw_allo)
    non_parts.append(raw_non)

allo_full = np.concatenate(allo_parts, axis=0)
non_full = np.concatenate(non_parts, axis=0)

allo_top10 = allo_full[:, col_indices]
non_top10 = non_full[:, col_indices]

n_allo = (~np.isnan(allo_top10).any(axis=1)).sum()
n_non = (~np.isnan(non_top10).any(axis=1)).sum()
print(f"Allosteric residues (no NaN row): {n_allo}")
print(f"Non-allosteric residues (no NaN row): {n_non}")

fig, axes = plt.subplots(2, 5, figsize=(22, 9))
axes = axes.flatten()

COLORS = {"Allosteric": "#e74c3c", "Non-allosteric": "#3498db"}

for i, (feat_name, ax) in enumerate(zip(top10_names, axes)):
    col_allo = allo_top10[:, i]
    col_non = non_top10[:, i]
    col_allo = col_allo[~np.isnan(col_allo)]
    col_non = col_non[~np.isnan(col_non)]

    bp = ax.boxplot(
        [col_allo, col_non],
        labels=["Allosteric", "Non-allo."],
        patch_artist=True,
        medianprops=dict(color="black", linewidth=2),
        whiskerprops=dict(linewidth=1.2),
        capprops=dict(linewidth=1.2),
        showfliers=False,
        widths=0.5,
    )
    bp["boxes"][0].set_facecolor(COLORS["Allosteric"])
    bp["boxes"][0].set_alpha(0.75)
    bp["boxes"][1].set_facecolor(COLORS["Non-allosteric"])
    bp["boxes"][1].set_alpha(0.75)

    ax.set_title(
        f"{feat_name}\n(imp={top10_importances[i]:.4f})",
        fontsize=8.5,
        fontweight="bold",
        pad=4,
    )
    ax.tick_params(axis="x", labelsize=8)
    ax.tick_params(axis="y", labelsize=7)
    ax.yaxis.grid(True, linestyle="--", alpha=0.5)
    ax.set_axisbelow(True)

from matplotlib.patches import Patch
legend_handles = [
    Patch(facecolor=COLORS["Allosteric"], alpha=0.75, label=f"Allosteric (n={n_allo})"),
    Patch(facecolor=COLORS["Non-allosteric"], alpha=0.75, label=f"Non-allosteric (n={n_non})"),
]
fig.legend(handles=legend_handles, loc="upper center", ncol=2, fontsize=11,
           bbox_to_anchor=(0.5, 1.01), frameon=True)

fig.suptitle(
    "Top 10 fpocket Features: Allosteric vs Non-Allosteric Residues",
    fontsize=14, fontweight="bold", y=1.05,
)
plt.tight_layout()
out_path = os.path.join(OUT_DIR, "fpocket_features_boxplots.png")
plt.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"Saved: {out_path}")
