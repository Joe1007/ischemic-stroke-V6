#!/usr/bin/env python3
"""Shared modelling helpers used by V6 (clinical loading, preprocessing, RFE, LightGBM, metrics).

Merged 2026-09-30 from the earlier iterations so V6 runs on its own. Only the definitions V6
actually calls are kept; their logic is unchanged. Origin of each block is marked below.
"""
from __future__ import annotations
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path

import lightgbm as lgb
import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    brier_score_loss,
    confusion_matrix,
    roc_auc_score,
)
from sklearn.neighbors import NearestNeighbors



# ---- from ischemic stroke analyseV2/leakage_free_experiment.py --------

ROOT = Path(__file__).resolve().parents[1]  # project root

DATA_DIR = ROOT / "data"

@dataclass(frozen=True)
class TaskSpec:
    target_column: str
    positive_rule: str
    final_features: int

TASKS = {
    "mRS": TaskSpec("mRS_discharge", ">2", 20),
    "NIHSS": TaskSpec("NIHSS_discharge", ">4", 35),
    "TE": TaskSpec("T/E", "E", 20),
    "LVO": TaskSpec("LVO-index", "1", 25),
}

def load_clinical() -> pd.DataFrame:
    """Reproduce source-notebook cleaning, but defer imputation to each fold."""
    clinical = pd.read_excel(DATA_DIR / "20231124is0155-500NIHSSfinal.xlsx", sheet_name=0)
    clinical = clinical.rename(columns={"Unnamed: 0": "ID"}).iloc[:, :-33].copy()
    delete_positions = [1, 2, 3, 11, 12, 13, 14, 15, 16, 17]
    clinical = clinical.drop(clinical.columns[delete_positions], axis=1)
    clinical = clinical.rename(
        columns={"LVO-index>50%": "LVO-index", "Tx (1tPA, 2IA)": "Tx"}
    )

    # Corrections present in the original notebook.
    clinical.iloc[[5, 155, 171], 35] = np.nan
    clinical.iloc[86, 19] = 1
    clinical.iloc[89, 23] = 1
    clinical.iloc[:, 12:] = clinical.iloc[:, 12:].apply(pd.to_numeric)
    clinical = clinical.drop(index=157).reset_index(drop=True)
    clinical["ID"] = clinical["ID"].astype(str)
    return clinical

def make_target(values: pd.Series, task: str) -> np.ndarray:
    if task == "mRS":
        return (pd.to_numeric(values) > 2).astype(np.int8).to_numpy()
    if task == "NIHSS":
        return (pd.to_numeric(values) > 4).astype(np.int8).to_numpy()
    if task == "TE":
        return (values.astype(str) == "E").astype(np.int8).to_numpy()
    if task == "LVO":
        return pd.to_numeric(values).astype(np.int8).to_numpy()
    raise KeyError(task)

def lightgbm_model(seed: int, n_jobs: int) -> lgb.LGBMClassifier:
    return lgb.LGBMClassifier(
        learning_rate=0.04,
        max_depth=3,
        num_leaves=6,
        reg_alpha=0.1,
        reg_lambda=1.0,
        colsample_bytree=0.5,
        random_state=seed,
        n_jobs=n_jobs,
        verbosity=-1,
        deterministic=True,
        force_col_wise=True,
    )

def fit_train_preprocessor(
    x_train: np.ndarray,
    x_test: np.ndarray,
    names: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    """Fit feature validity, medians and variance filtering on training rows only."""
    finite_counts = np.isfinite(x_train).sum(axis=0)
    keep = finite_counts > 0
    x_train = x_train[:, keep]
    x_test = x_test[:, keep]
    names = names[keep]

    medians = np.nanmedian(np.where(np.isfinite(x_train), x_train, np.nan), axis=0)
    train_bad = ~np.isfinite(x_train)
    test_bad = ~np.isfinite(x_test)
    if train_bad.any():
        x_train = x_train.copy()
        x_train[train_bad] = np.take(medians, np.where(train_bad)[1])
    if test_bad.any():
        x_test = x_test.copy()
        x_test[test_bad] = np.take(medians, np.where(test_bad)[1])

    varying = np.ptp(x_train, axis=0) > 0
    audit = {
        "raw_features": int(keep.size),
        "all_missing_train_removed": int((~keep).sum()),
        "zero_variance_train_removed": int((~varying).sum()),
        "eligible_features": int(varying.sum()),
        "test_values_imputed": int(test_bad[:, varying].sum()),
    }
    return x_train[:, varying], x_test[:, varying], names[varying], audit

def smote_training_only(
    x: np.ndarray, y: np.ndarray, seed: int, k_neighbors: int = 3
) -> tuple[np.ndarray, np.ndarray]:
    """Deterministic binary SMOTE equivalent for the already-selected feature set."""
    counts = np.bincount(y, minlength=2)
    minority = int(np.argmin(counts))
    majority_count = int(counts.max())
    minority_rows = x[y == minority]
    needed = majority_count - minority_rows.shape[0]
    if needed <= 0:
        return x, y
    k = min(k_neighbors, minority_rows.shape[0] - 1)
    if k < 1:
        raise ValueError("SMOTE requires at least two minority training samples")

    neighbors = NearestNeighbors(n_neighbors=k + 1).fit(minority_rows)
    neighbor_index = neighbors.kneighbors(
        minority_rows, n_neighbors=k + 1, return_distance=False
    )[:, 1:]
    rng = np.random.default_rng(seed)
    bases = rng.integers(0, minority_rows.shape[0], size=needed)
    neighbor_choices = rng.integers(0, k, size=needed)
    selected_neighbors = neighbor_index[bases, neighbor_choices]
    gaps = rng.random((needed, 1), dtype=np.float32)
    synthetic = minority_rows[bases] + gaps * (
        minority_rows[selected_neighbors] - minority_rows[bases]
    )
    return (
        np.vstack([x, synthetic]).astype(np.float32, copy=False),
        np.concatenate([y, np.full(needed, minority, dtype=y.dtype)]),
    )

def classification_metrics(y: np.ndarray, probability: np.ndarray) -> dict[str, float]:
    prediction = (probability >= 0.5).astype(np.int8)
    tn, fp, fn, tp = confusion_matrix(y, prediction, labels=[0, 1]).ravel()
    return {
        "auroc": float(roc_auc_score(y, probability)),
        "auprc": float(average_precision_score(y, probability)),
        "accuracy": float(accuracy_score(y, prediction)),
        "sensitivity": float(tp / (tp + fn)) if tp + fn else math.nan,
        "specificity": float(tn / (tn + fp)) if tn + fp else math.nan,
    }


# ---- from ischemic stroke analyseV2/atlas_group_experiment.py ---------

def clinical_feature_names(mode: str, task: str) -> list[str]:
    paper = [
        "age", "sex", "onset-to-MRI", "2D-EF", "AF", "CHF", "CAD", "HTN",
        "Lipid", "DM", "cancer", "RecurSt", "Tx", "pneumonia", "UTI",
        "smoking", "alcohol", "iSBP", "iDBP", "iHR", "iRR", "iBT", "iFS",
        "HbA1c", "Cr", "BUN", "Tchol", "LDL", "HDL", "TG", "WBC", "Lym",
        "Seg", "Hb", "Hct", "wSBP", "wDBP", "wHR", "BH", "BW", "BMI",
    ]
    if mode == "none":
        return []
    # The source notebook deliberately excludes AF from the T/E model because it
    # is close to being part of the etiologic label definition.
    if task == "TE":
        paper = [column for column in paper if column != "AF"]
    if mode == "paper":
        return paper
    if mode == "preoutcome":
        # Admission severity is known before discharge and is therefore legitimate
        # for a discharge-outcome model, unlike discharge scores themselves.
        return ["mRS_admission", "NIHSS_admission", *paper]
    raise ValueError(mode)


# ---- from ischemic stroke analyseV4/run_v4.py -------------------------

CORE_C = ('age', 'sex', 'AF', 'CHF', 'CAD', 'HTN', 'Lipid', 'DM',
          'cancer', 'RecurSt', 'smoking', 'alcohol')

def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

def write_json(path, value):
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding='utf-8')
    temp.replace(path)

def valid_ids(ids):
    if len(ids) != len(set(ids)):
        raise ValueError('Duplicate patient IDs: cannot guarantee patient-level separation')

def standardize(train, test):
    mu = train.mean(axis=0, dtype=np.float64)
    scale = train.std(axis=0, dtype=np.float64)
    scale[scale < 1e-8] = 1.
    return ((train-mu)/scale).astype('float32'), ((test-mu)/scale).astype('float32')

def staged_rfe(train, y, names, keep, seed, jobs):
    """V1 elimination targets and steps, fixed model seed, stable tie breaking."""
    active = np.arange(train.shape[1])
    trace = []
    for target, step in ((10000, 40000), (950, 2000), (50, 150), (keep, 1)):
        target = max(target, keep)
        while len(active) > target:
            model = lightgbm_model(seed, jobs)
            model.fit(train[:, active], y)
            rank = np.argsort(-model.feature_importances_, kind='stable')
            active = active[rank[:max(target, len(active)-step)]]
            trace.append(len(active))
    return active, trace

def fit_predict(train, test, y, seed, jobs):
    tx, ty = smote_training_only(train, y, seed)
    model = lightgbm_model(seed, jobs)
    model.fit(tx, ty)
    return model.predict_proba(test)[:, 1]

def clinical_matrix(frame, ids, mode, task):
    cols = list(CORE_C) if mode == 'core' else clinical_feature_names('paper', task)
    if task == 'TE':
        cols = [c for c in cols if c != 'AF']
    return (frame.set_index('ID').loc[ids, cols].apply(pd.to_numeric).to_numpy(dtype='float32'),
            np.array(['clinical|' + c for c in cols]))

def metrics(y, p):
    return {**classification_metrics(y, p), 'brier': float(brier_score_loss(y, p))}

def uncertainty(y, predictions, reference, iterations, seed):
    """Conditional, paired patient bootstrap; not a training-pipeline CI."""
    rng = np.random.default_rng(seed)
    groups = [np.flatnonzero(y == i) for i in (0, 1)]
    aucs = {m: [] for m in predictions}
    for _ in range(iterations):
        ix = np.concatenate([rng.choice(g, len(g), replace=True) for g in groups])
        for m, p in predictions.items():
            aucs[m].append(roc_auc_score(y[ix], p[ix]))
    rows = []
    for m, p in predictions.items():
        lo, hi = np.quantile(aucs[m], [.025, .975])
        delta = np.asarray(aucs[m])-np.asarray(aucs[reference])
        dl, dh = np.quantile(delta, [.025, .975])
        rows.append({'method': m, **metrics(y, p), 'auc_low': lo, 'auc_high': hi,
                     'reference': reference,
                     'delta_auc': roc_auc_score(y,p)-roc_auc_score(y,predictions[reference]),
                     'delta_low': dl, 'delta_high': dh})
    return rows


# ---- from ischemic stroke analyseV5/model_v5.py -----------------------

def signature(train,test,names,y,keep,seed,jobs):
    tx,vx,nn,prep=fit_train_preprocessor(train,test,names)
    tx,vx=standardize(tx,vx)
    selected,trace=staged_rfe(tx,y,nn,keep,seed,jobs)
    return tx[:,selected],vx[:,selected],{'preprocessing':prep,'selected_features':nn[selected].tolist(),'rfe_remaining':trace}
