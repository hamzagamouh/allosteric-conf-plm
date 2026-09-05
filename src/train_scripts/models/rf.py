from sklearn.ensemble import RandomForestClassifier
from train_scripts.models.metrics import evaluate_metrics


def train_rf(train_feats, y_train, val_feats, y_val,
             n_estimators=10, n_jobs=-1, random_state=454356):
    """Train a Random Forest and return (train_metrics, val_metrics, feature_importances)."""
    rf = RandomForestClassifier(
        n_estimators=n_estimators, random_state=random_state, n_jobs=n_jobs
    )
    rf.fit(train_feats, y_train.ravel())
    train_metrics = evaluate_metrics(rf, train_feats, y_train)
    val_metrics = evaluate_metrics(rf, val_feats, y_val)
    return train_metrics, val_metrics, rf.feature_importances_
