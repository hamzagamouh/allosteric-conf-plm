"""Compare distribution of ligand types binding allosteric vs non-allosteric residues.

Uses fold-0 train + val JSON files (each protein appears exactly once).
Classifies ligand 3-letter codes via ASD_2023.tsv, with fallback heuristics.
"""
import csv
import json
import os
from collections import Counter, defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJDIR = "/auto/vestec1-elixir/home/hamzagamouh/allosteric-conf-plm"
ASD_TSV = f"{PROJDIR}/asd_data/ASD_2023.tsv"
TRAIN_JSON = f"{PROJDIR}/asd_data/train/allosteric_ligands_train_fold_0.json"
VAL_JSON = f"{PROJDIR}/asd_data/val/allosteric_ligands_val_fold_0.json"
OUT_DIR = f"{PROJDIR}/visualizations"

# ── Build alias → normalised class from ASD_2023.tsv ─────────────────────────

def is_ion(raw):
    return raw.strip().lower().split(";")[0].strip() in ("ion", "lon")

alias_is_ion = {}
with open(ASD_TSV) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        alias = row.get("modulator_alias", "").strip().upper()
        if alias:
            alias_is_ion[alias] = is_ion(row.get("modulator_class", ""))

# ── Fallback: known PDB ion codes not in ASD ─────────────────────────────────
KNOWN_IONS = {
    "MG", "ZN", "CA", "NA", "CL", "FE", "MN", "CU", "K", "SE", "CO", "NI",
    "CD", "PB", "HG", "BA", "SR", "RB", "CS", "LI", "BE", "AL", "CR", "MO",
    "W", "RE", "OS", "IR", "PT", "AU", "TL", "IN", "SN", "BI", "PO", "EU",
    "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "LA", "CE", "PR", "ND",
    "SM", "Y", "SC", "V", "TI", "GA", "GE", "AS", "RU", "RH", "PD",
    "AG", "SB", "TE", "I", "XE", "BR", "KR",
    "PO4", "SO4", "NO3", "CL4", "CL3",
}

def classify_code(code):
    code = code.upper()
    if code in alias_is_ion:
        return "Ion" if alias_is_ion[code] else "Non-Ion"
    return "Ion" if code in KNOWN_IONS else "Non-Ion"

def extract_code(lig_id):
    """Extract 3-letter code from 'pdb-chain-code-resnum'."""
    parts = lig_id.split("-")
    return parts[2].upper() if len(parts) >= 3 else lig_id.upper()

# ── Load fold-0 data (train + val → all proteins, each once) ─────────────────

def load_and_merge():
    merged = {}
    for path in (TRAIN_JSON, VAL_JSON):
        d = json.load(open(path))
        merged.update(d)
    return merged

data = load_and_merge()
print(f"Total proteins: {len(data)}")

# Per-protein unique ligand codes for allosteric and non-allosteric sites
allo_codes_all = []
non_codes_all = []

for unp, entry in data.items():
    allo_ligs = entry.get("Allosteric ligands", [])
    non_ligs = entry.get("Non-Allosteric ligands", [])

    allo_codes = {extract_code(l) for l in allo_ligs}
    non_codes = {extract_code(l) for l in non_ligs}

    allo_codes_all.extend(allo_codes)
    non_codes_all.extend(non_codes)

# ── Classify and count ────────────────────────────────────────────────────────
CATEGORIES = ["Ion", "Non-Ion"]

def count_categories(codes):
    cnt = Counter(classify_code(c) for c in codes)
    return {cat: cnt.get(cat, 0) for cat in CATEGORIES}

allo_counts = count_categories(allo_codes_all)
non_counts = count_categories(non_codes_all)

allo_total = sum(allo_counts.values())
non_total = sum(non_counts.values())

allo_pct = {k: 100 * v / allo_total for k, v in allo_counts.items()}
non_pct = {k: 100 * v / non_total for k, v in non_counts.items()}

print("\nAllosteric ligand type counts:", allo_counts, f"(total={allo_total})")
print("Non-allosteric ligand type counts:", non_counts, f"(total={non_total})")
print("\nAllosteric %:", {k: f"{v:.1f}" for k, v in allo_pct.items()})
print("Non-allosteric %:", {k: f"{v:.1f}" for k, v in non_pct.items()})

# ── Plot ──────────────────────────────────────────────────────────────────────
CAT_COLORS = {
    "Ion":     "#e74c3c",
    "Non-Ion": "#3498db",
}

fig, ax = plt.subplots(figsize=(7, 6))

# Grouped bar chart (percentages)
x = np.arange(len(CATEGORIES))
width = 0.35

bars_allo = ax.bar(x - width / 2,
                   [allo_pct[c] for c in CATEGORIES],
                   width, label=f"Allosteric (n={allo_total})",
                   color=[CAT_COLORS[c] for c in CATEGORIES],
                   alpha=0.85, edgecolor="white", linewidth=0.8)
bars_non = ax.bar(x + width / 2,
                  [non_pct[c] for c in CATEGORIES],
                  width, label=f"Non-allosteric (n={non_total})",
                  color=[CAT_COLORS[c] for c in CATEGORIES],
                  alpha=0.4, edgecolor="black", linewidth=0.8,
                  hatch="//")

ax.set_xticks(x)
ax.set_xticklabels(CATEGORIES, fontsize=11)
ax.set_ylabel("Percentage (%)", fontsize=12)
ax.set_title("Ligand Type Distribution: Allosteric vs Non-Allosteric Residues", fontsize=12, fontweight="bold")
ax.yaxis.grid(True, linestyle="--", alpha=0.5)
ax.set_axisbelow(True)
ax.legend(fontsize=10)

# Add value labels on bars
for bar in bars_allo:
    h = bar.get_height()
    if h > 1:
        ax.text(bar.get_x() + bar.get_width() / 2, h + 0.5,
                f"{h:.1f}%", ha="center", va="bottom", fontsize=8, fontweight="bold")
for bar in bars_non:
    h = bar.get_height()
    if h > 1:
        ax.text(bar.get_x() + bar.get_width() / 2, h + 0.5,
                f"{h:.1f}%", ha="center", va="bottom", fontsize=8)

plt.tight_layout()
out_path = os.path.join(OUT_DIR, "ligand_type_comparison.png")
plt.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"\nSaved: {out_path}")
