"""MACCS fingerprint t-SNE embedding: allosteric vs non-allosteric ligands.

Two plots side by side: training set (fold 0) and validation set (fold 0).
Large splits are stratified-subsampled to keep t-SNE tractable.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

PROJDIR = "/auto/vestec1-elixir/home/hamzagamouh/allosteric-conf-plm"
FEATURES_DIR = f"{PROJDIR}/features"
OUT_DIR = f"{PROJDIR}/visualizations"

MAX_SAMPLES = 10_000   # cap per split to keep t-SNE runtime manageable
PERPLEXITY   = 40
RANDOM_STATE = 42


def load_maccs(mode, fold=0):
    allo     = np.load(f"{FEATURES_DIR}/{mode}_allo_maccs_fold_{fold}.npy")
    non_allo = np.load(f"{FEATURES_DIR}/{mode}_non_allo_maccs_fold_{fold}.npy")
    X = np.concatenate([allo, non_allo], axis=0).astype(float)
    y = np.array([1] * len(allo) + [0] * len(non_allo))
    return X, y


def stratified_subsample(X, y, max_n, rng):
    if len(y) <= max_n:
        return X, y
    classes, counts = np.unique(y, return_counts=True)
    fractions = counts / counts.sum()
    indices = []
    for cls, frac in zip(classes, fractions):
        cls_idx = np.where(y == cls)[0]
        n_take  = int(round(max_n * frac))
        n_take  = min(n_take, len(cls_idx))
        indices.append(rng.choice(cls_idx, size=n_take, replace=False))
    idx = np.concatenate(indices)
    rng.shuffle(idx)
    return X[idx], y[idx]


def embed_tsne(X, rng):
    # PCA first to reduce noise and speed up t-SNE
    n_pca = min(50, X.shape[1], X.shape[0] - 1)
    pca   = PCA(n_components=n_pca, random_state=RANDOM_STATE)
    X_pca = pca.fit_transform(X)
    print(f"    PCA variance explained ({n_pca} components): "
          f"{pca.explained_variance_ratio_.sum():.1%}")
    tsne = TSNE(
        n_components=2,
        perplexity=PERPLEXITY,
        random_state=RANDOM_STATE,
        n_iter=1000,
        init="pca",
        learning_rate=200.0,
    )
    return tsne.fit_transform(X_pca)


rng = np.random.default_rng(RANDOM_STATE)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, mode, title in zip(axes,
                            ("train", "val"),
                            ("Training set", "Validation set")):
    print(f"\n── {title} ──")
    X, y = load_maccs(mode)
    print(f"  Full: allo={y.sum()}, non-allo={(y==0).sum()}, total={len(y)}")

    X_sub, y_sub = stratified_subsample(X, y, MAX_SAMPLES, rng)
    n_allo = y_sub.sum()
    n_non  = (y_sub == 0).sum()
    print(f"  After subsample: allo={n_allo}, non-allo={n_non}, total={len(y_sub)}")

    print("  Running PCA + t-SNE …")
    emb = embed_tsne(X_sub, rng)

    non_mask  = y_sub == 0
    allo_mask = y_sub == 1

    ax.scatter(emb[non_mask,  0], emb[non_mask,  1],
               c="#3498db", alpha=0.35, s=6,
               label=f"Non-allosteric (n={n_non})", rasterized=True)
    ax.scatter(emb[allo_mask, 0], emb[allo_mask, 1],
               c="#e74c3c", alpha=0.65, s=6,
               label=f"Allosteric (n={n_allo})", rasterized=True)

    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_xlabel("t-SNE 1", fontsize=11)
    ax.set_ylabel("t-SNE 2", fontsize=11)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    ax.legend(fontsize=10, markerscale=3, framealpha=0.85)

fig.suptitle(
    "MACCS Fingerprint Embedding (bitwise OR): Allosteric vs Non-Allosteric Residues\n"
    "(PCA → t-SNE, fold 0)",
    fontsize=13, fontweight="bold",
)
plt.tight_layout()
out_path = os.path.join(OUT_DIR, "maccs_embedding.png")
plt.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"\nSaved: {out_path}")
