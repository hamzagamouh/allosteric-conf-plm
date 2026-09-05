import random
import numpy as np
from sklearn.metrics import (
    matthews_corrcoef, f1_score, accuracy_score,
    precision_score, recall_score, roc_auc_score,
)


def evaluate_metrics(clf, X, y):
    """Return a dict with MCC, F1, accuracy, precision, recall, and ROC-AUC."""
    y_true = y.ravel()
    y_pred = clf.predict(X)
    y_prob = clf.predict_proba(X)[:, 1]
    return {
        "mcc":       float(matthews_corrcoef(y_true, y_pred)),
        "f1":        float(f1_score(y_true, y_pred, zero_division=0)),
        "accuracy":  float(accuracy_score(y_true, y_pred)),
        "precision": float(precision_score(y_true, y_pred, zero_division=0)),
        "recall":    float(recall_score(y_true, y_pred, zero_division=0)),
        "roc_auc":   float(roc_auc_score(y_true, y_prob)),
    }


def evaluate_mcc(clf, X, y):
    return matthews_corrcoef(y.ravel(), clf.predict(X))


def drop_nan_rows(X, ref):
    return X[~np.isnan(ref).any(axis=1)]


def balance_subsample(allo_feats, non_allo_feats, seed=4652):
    """Downsample non-allosteric class to match allosteric class size."""
    random.seed(seed)
    selected = random.sample(range(len(non_allo_feats)), len(allo_feats))
    X = np.concatenate([allo_feats, non_allo_feats[selected]], axis=0)
    y = np.array([1] * len(allo_feats) + [0] * len(allo_feats)).reshape(-1, 1)
    return X, y


def build_val_arrays(allo_feats, non_allo_feats):
    X = np.concatenate([allo_feats, non_allo_feats], axis=0)
    y = np.array([1] * len(allo_feats) + [0] * len(non_allo_feats)).reshape(-1, 1)
    return X, y
