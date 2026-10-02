
import h5py, numpy as np, scipy.sparse as sp, anndata as ad, pandas as pd, os

def _index_key(grp):
    """h5ad declares its index column in the group's _index attribute; the literal
    name varies across CELLxGENE releases, so resolve it rather than assuming."""
    ik = grp.attrs.get("_index", "_index")
    if isinstance(ik, bytes):
        ik = ik.decode()
    if ik not in grp:
        ik = next((k for k in ("_index", "index", "obs_names", "var_names") if k in grp), None)
    return ik

def _read_obs_col(f, name):
    o = f["obs"][name]
    if isinstance(o, h5py.Group):                      # categorical
        cats = o["categories"][:]
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in cats], dtype=object)
        codes = o["codes"][:]
        out = np.where(codes >= 0, cats[np.clip(codes, 0, None)], None)
        return pd.Series(out)
    v = o[:]
    if v.dtype.kind == "O":
        v = np.array([x.decode() if isinstance(x, bytes) else x for x in v], dtype=object)
    return pd.Series(v)

def read_obs(path, cols=None):
    with h5py.File(path, "r") as f:
        keys = list(f["obs"].keys())
        if cols is None:
            cols = [k for k in keys if not k.startswith("_")]
        cols = [c for c in cols if c in keys]
        df = pd.DataFrame({c: _read_obs_col(f, c) for c in cols})
        idx = f["obs"][_index_key(f["obs"])][:]
        df.index = [x.decode() if isinstance(x, bytes) else x for x in idx]
    return df

def _var_frame(f, root):
    g = f[root]["var"] if root == "raw" else f["var"]
    def col(n):
        o = g[n]
        if isinstance(o, h5py.Group):
            cats = np.array([c.decode() if isinstance(c, bytes) else c for c in o["categories"][:]], dtype=object)
            return cats[np.clip(o["codes"][:], 0, None)]
        v = o[:]
        return np.array([x.decode() if isinstance(x, bytes) else x for x in v], dtype=object) if v.dtype.kind == "O" else v
    idx = col(_index_key(g))
    d = {"var_index": idx}
    for n in ("ensembl_id", "feature_name", "feature_biotype"):
        if n in g: d[n] = col(n)
    return pd.DataFrame(d)

def read_rows_csr(path, group, rows):
    """Read selected rows of an h5ad CSR matrix without loading the whole thing."""
    rows = np.asarray(sorted(rows))
    with h5py.File(path, "r") as f:
        g = f[group]
        indptr = g["indptr"][:]
        ncol = int(g.attrs.get("shape", g.attrs.get("h5sparse_shape"))[1])
        dat, ind, newptr = [], [], [0]
        # group contiguous runs to cut the number of HDF5 reads
        runs, s = [], 0
        while s < len(rows):
            e = s
            while e + 1 < len(rows) and rows[e + 1] == rows[e] + 1: e += 1
            runs.append((rows[s], rows[e])); s = e + 1
        for a, b in runs:
            lo, hi = indptr[a], indptr[b + 1]
            blk_d = g["data"][lo:hi]; blk_i = g["indices"][lo:hi]
            off = indptr[a:b + 2] - lo
            for j in range(b - a + 1):
                dat.append(blk_d[off[j]:off[j + 1]]); ind.append(blk_i[off[j]:off[j + 1]])
                newptr.append(newptr[-1] + off[j + 1] - off[j])
    X = sp.csr_matrix((np.concatenate(dat), np.concatenate(ind), np.array(newptr)),
                      shape=(len(rows), ncol))
    return X, rows

def build_panel(path, out, name, assay_keep=None, per_type=300, cap=6000,
                min_counts=500, min_genes=250, max_pct_mt=20.0, min_type_n=40,
                disease_col="disease", seed=0):
    """Subsample one CELLxGENE h5ad into a compact raw-count panel dataset."""
    obs = read_obs(path)
    n0 = len(obs)
    keep = pd.Series(True, index=obs.index)
    if assay_keep is not None and "assay" in obs:
        keep &= obs["assay"].isin(assay_keep if isinstance(assay_keep, (list, tuple)) else [assay_keep])
    if "total_counts" in obs:      keep &= pd.to_numeric(obs["total_counts"], errors="coerce").fillna(0) >= min_counts
    if "n_genes_by_counts" in obs: keep &= pd.to_numeric(obs["n_genes_by_counts"], errors="coerce").fillna(0) >= min_genes
    if "pct_counts_mt" in obs:     keep &= pd.to_numeric(obs["pct_counts_mt"], errors="coerce").fillna(0) <= max_pct_mt
    obs_f = obs[keep]
    # stratified subsample: cell_type x donor x disease
    rng = np.random.default_rng(seed)
    strat = ["cell_type"]
    if disease_col in obs_f: strat.append(disease_col)
    counts = obs_f.groupby("cell_type", observed=True).size()
    ok_types = counts[counts >= min_type_n].index
    obs_f = obs_f[obs_f["cell_type"].isin(ok_types)]
    picks = []
    for key, grp in obs_f.groupby(strat, observed=True):
        take = min(per_type, len(grp))
        picks.append(grp.sample(take, random_state=int(rng.integers(1e9))))
    sel = pd.concat(picks) if picks else obs_f.iloc[:0]
    if len(sel) > cap:
        sel = sel.groupby(strat, observed=True, group_keys=False).apply(
            lambda g: g.sample(max(1, int(round(len(g) * cap / len(sel)))), random_state=0))
    pos = {c: i for i, c in enumerate(obs.index)}
    rows = [pos[c] for c in sel.index]
    with h5py.File(path, "r") as f:
        grp = "raw/X" if "raw" in f and "X" in f["raw"] else "X"
        vf = _var_frame(f, "raw" if grp.startswith("raw") else "X")
    X, rows_sorted = read_rows_csr(path, grp, rows)
    order = [obs.index[i] for i in rows_sorted]
    A = ad.AnnData(X=X, obs=sel.loc[order].copy(), var=vf.set_index("var_index"))
    A.obs["dataset"] = name
    A.uns["source_file"] = os.path.basename(path)
    A.uns["n_cells_original"] = int(n0)
    A.write_h5ad(out, compression="gzip")
    return A
