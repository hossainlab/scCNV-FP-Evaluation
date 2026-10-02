"""In-silico aneuploid positive controls.

Held-out healthy cells are given arm-level copy-number events by binomial thinning
(losses) or augmentation (gains) of their own UMI counts, then re-thinned to the
original library size. This preserves each cell's type identity, depth and sparsity
while imposing a known karyotype, giving ground-truth sensitivity and a genuine
alternative distribution against which FDR control can be checked.
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp

# acrocentric short arms carry no protein-coding genes worth modelling
EXCLUDE_ARMS = {("13", "p"), ("14", "p"), ("15", "p"), ("21", "p"), ("22", "p")}
AUTOSOMES = [str(i) for i in range(1, 23)]


def arm_list(arms_csv="chrom_arms.csv"):
    a = pd.read_csv(arms_csv, dtype={"chrom": str})
    out = []
    for _, r in a.iterrows():
        if r["chrom"] not in AUTOSOMES:
            continue
        for arm, s, e in (("p", r["p_start"], r["p_end"]), ("q", r["q_start"], r["q_end"])):
            if (r["chrom"], arm) in EXCLUDE_ARMS or not np.isfinite(s):
                continue
            out.append({"chrom": r["chrom"], "arm": arm, "start": int(s), "end": int(e)})
    return pd.DataFrame(out)


def make_karyotype(rng, arms, n_events, gain_frac=0.5):
    """Draw `n_events` distinct arm-level events; copy number 1 (loss) or 3 (gain), rarely 4."""
    idx = rng.choice(len(arms), size=n_events, replace=False)
    ev = []
    for i in idx:
        a = arms.iloc[i]
        if rng.random() < gain_frac:
            cn = 4 if rng.random() < 0.15 else 3
        else:
            cn = 1
        ev.append({"chrom": a["chrom"], "arm": a["arm"], "start": int(a["start"]),
                   "end": int(a["end"]), "cn": cn, "ratio": cn / 2.0})
    return ev


def gene_ratio_vector(karyotype, gene_df):
    """Map an arm-level karyotype onto a per-gene expression-dosage ratio."""
    r = np.ones(len(gene_df), dtype=np.float64)
    ch = gene_df["chrom"].astype(str).values
    mid = gene_df["mid"].values
    for ev in karyotype:
        m = (ch == ev["chrom"]) & (mid >= ev["start"]) & (mid <= ev["end"])
        r[m] = ev["ratio"]
    return r


def _thin_row(vals, p, rng):
    """Binomial thin/augment a count vector toward ratio p (elementwise, p>=0)."""
    vals = vals.astype(np.int64)
    out = vals.copy()
    lo = p < 1.0
    if lo.any():
        out[lo] = rng.binomial(vals[lo], np.clip(p[lo], 0.0, 1.0))
    hi = p > 1.0
    if hi.any():
        extra_p = np.clip(p[hi] - 1.0, 0.0, 1.0)
        whole = np.floor(p[hi] - 1.0).astype(np.int64)
        out[hi] = vals[hi] + whole * vals[hi] + rng.binomial(vals[hi], extra_p)
    return out


def spike_matrix(X, gene_r, rng, renorm_depth=True):
    """Apply a per-gene dosage ratio to every row of a CSR count matrix."""
    X = sp.csr_matrix(X)
    X.eliminate_zeros()
    dat = X.data.astype(np.int64).copy()
    ind = X.indices
    ptr = X.indptr
    newdat = np.empty_like(dat)
    for i in range(X.shape[0]):
        s, e = ptr[i], ptr[i + 1]
        if s == e:
            continue
        v = dat[s:e]
        p = gene_r[ind[s:e]]
        w = _thin_row(v, p, rng)
        if renorm_depth:
            tot0, tot1 = v.sum(), w.sum()
            if tot1 > 0 and tot0 > 0:
                f = tot0 / tot1
                if f < 1.0:
                    w = rng.binomial(w, f)
                elif f > 1.0:
                    w = w + rng.binomial(w, min(f - 1.0, 1.0))
        newdat[s:e] = w
    out = sp.csr_matrix((newdat.astype(np.float32), ind.copy(), ptr.copy()), shape=X.shape)
    out.eliminate_zeros()
    return out


def build_spikein(adata, gene_df, rng, burdens=(1, 2, 3, 5, 8), n_per_burden=150,
                  aneuploid_frac=None):
    """Return (X_new, obs_new) where a subset of cells carry known arm-level events.

    `burdens` is the number of arm events per simulated clone; one clone per burden
    level. Cells are drawn without replacement from `adata`. Remaining cells are kept
    as diploid companions so each simulated sample is a realistic mixture.
    """
    arms = arm_list()
    n = adata.n_obs
    order = rng.permutation(n)
    need = len(burdens) * n_per_burden
    assert n >= need, f"need {need} cells, have {n}"
    blocks, obs_rows, ptr = [], [], 0
    truth = []
    for b in burdens:
        sel = order[ptr:ptr + n_per_burden]; ptr += n_per_burden
        kary = make_karyotype(rng, arms, b)
        gr = gene_ratio_vector(kary, gene_df)
        Xs = spike_matrix(adata.X[sel], gr, rng)
        blocks.append(Xs)
        o = adata.obs.iloc[sel].copy()
        o["spike_burden"] = b
        o["spike_truth"] = "aneuploid"
        o["spike_clone"] = f"clone_b{b}"
        o["spike_karyotype"] = ";".join(f"{e['chrom']}{e['arm']}:CN{e['cn']}" for e in kary)
        obs_rows.append(o)
        truth.append({"clone": f"clone_b{b}", "burden": b, "n_cells": int(n_per_burden),
                      "events": kary,
                      "frac_genome": float(np.mean(gr != 1.0))})
    rest = order[ptr:]
    if len(rest):
        blocks.append(sp.csr_matrix(adata.X[rest]))
        o = adata.obs.iloc[rest].copy()
        o["spike_burden"] = 0
        o["spike_truth"] = "diploid"
        o["spike_clone"] = "diploid"
        o["spike_karyotype"] = ""
        obs_rows.append(o)
    X = sp.vstack(blocks, format="csr")
    obs = pd.concat(obs_rows)
    obs.index = [f"spk_{i}_{c}" for i, c in enumerate(obs.index)]
    return X, obs, truth
