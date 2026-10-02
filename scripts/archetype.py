"""Label-free conditioning: assign each cell to a transcriptional archetype.

Cell-type-conditional calibration needs the query cell's identity to exist in the
reference. Annotation vocabularies differ between studies, and disease states carry
labels no healthy atlas contains, so label matching fails exactly where it is needed.
Assigning each query cell to its nearest healthy reference centroid removes the
dependence on a shared label vocabulary.
"""
import os, glob
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, scanpy as sc


def _lognorm(A):
    B = A.copy()
    sc.pp.normalize_total(B, target_sum=1e4)
    sc.pp.log1p(B)
    return B


def symbols(A):
    return (A.var["feature_name"].astype(str).values if "feature_name" in A.var
            else np.asarray(A.var_names, dtype=str))


def build_centroids(paths, min_n=50):
    """Per-cell-type mean log-normalised profile, pooled over reference datasets."""
    univ = None
    for p in paths:
        s = set(symbols(ad.read_h5ad(p, backed="r")))
        univ = s if univ is None else (univ & s)
    univ = np.array(sorted(univ))
    idx = {g: i for i, g in enumerate(univ)}
    rows, names, counts = [], [], []
    for p in paths:
        A = _lognorm(ad.read_h5ad(p))
        sy = symbols(A)
        keep = np.array([g in idx for g in sy])
        cols = np.array([idx[g] for g in sy[keep]])
        X = A.X[:, keep]
        X = X.tocsc() if sp.issparse(X) else X
        ds = os.path.splitext(os.path.basename(p))[0]
        ct = A.obs["cell_type"].astype(str).values
        for t in pd.unique(ct):
            m = ct == t
            if m.sum() < min_n:
                continue
            v = np.zeros(len(univ))
            sub = X[m]
            mu = np.asarray(sub.mean(0)).ravel() if sp.issparse(sub) else np.asarray(sub).mean(0)
            v[cols] = mu
            rows.append(v); names.append(f"{ds}|{t}"); counts.append(int(m.sum()))
    C = pd.DataFrame(np.vstack(rows), index=names, columns=univ)
    return C, pd.Series(counts, index=names, name="n_cells")


def select_marker_genes(C, n_top=2000):
    """Genes that differ most across archetypes carry the identity signal."""
    v = C.values.std(0)
    return C.columns.values[np.argsort(v)[::-1][:n_top]]


def assign(query_path, C, genes):
    """Best-correlating archetype per query cell, with the margin to the runner-up."""
    A = _lognorm(ad.read_h5ad(query_path))
    sy = symbols(A)
    pos = pd.Series(np.arange(len(sy)), index=sy)
    pos = pos[~pos.index.duplicated()]
    shared = [g for g in genes if g in pos.index]
    Q = A.X[:, pos[shared].values]
    Q = np.asarray(Q.todense()) if sp.issparse(Q) else np.asarray(Q)
    M = C.loc[:, shared].values
    Qz = (Q - Q.mean(1, keepdims=True)) / (Q.std(1, keepdims=True) + 1e-9)
    Mz = (M - M.mean(1, keepdims=True)) / (M.std(1, keepdims=True) + 1e-9)
    R = Qz @ Mz.T / len(shared)
    o = np.argsort(R, 1)
    best = o[:, -1]
    return pd.DataFrame({"cell": np.asarray(A.obs_names, dtype=str),
                         "archetype": C.index.values[best],
                         "archetype_corr": R[np.arange(len(best)), best],
                         "archetype_margin": R[np.arange(len(best)), best] - R[np.arange(len(best)), o[:, -2]],
                         "n_shared_genes": len(shared)})
