"""Cell-type-conditional empirical-null calibration of malignant-cell scores.

A caller emits a per-cell aneuploidy score S. Practitioners threshold S globally, or
take the caller's internal two-class split. Both ignore that the null distribution of S
depends strongly on cell type: cells whose transcriptome is dominated by physically
clustered co-regulated genes (plasma cells, keratinocytes, erythroid precursors) carry
a high S with no copy-number change at all.

The fix conditions the null on cell identity and on sequencing depth, converts the score
to an empirical right-tail p-value against cancer-free reference cells of the same
identity, and controls FDR across cells within a sample.
"""
import numpy as np
import pandas as pd
from scipy import stats


# ---------------------------------------------------------------- depth adjustment
def fit_depth_model(ref, score_col="score", depth_col="log_depth", type_col="cell_type",
                    min_n=40):
    """Per-cell-type linear dependence of score on log depth, fitted on reference cells."""
    models = {}
    for t, g in ref.groupby(type_col, observed=True):
        if len(g) < min_n:
            continue
        x, y = g[depth_col].values, g[score_col].values
        if np.nanstd(x) < 1e-9:
            models[t] = (0.0, float(np.nanmean(y)))
            continue
        sl, ic = np.polyfit(x, y, 1)
        models[t] = (float(sl), float(ic))
    xall, yall = ref[depth_col].values, ref[score_col].values
    models["__pooled__"] = tuple(float(v) for v in np.polyfit(xall, yall, 1)) if np.nanstd(xall) > 1e-9 \
        else (0.0, float(np.nanmean(yall)))
    return models


def apply_depth_model(df, models, score_col="score", depth_col="log_depth",
                      type_col="cell_type"):
    """Residualise score on depth using the cell-type model (pooled fallback)."""
    out = np.empty(len(df))
    tv = df[type_col].astype(str).values
    dv = df[depth_col].values
    sv = df[score_col].values
    for i in range(len(df)):
        sl, ic = models.get(tv[i], models["__pooled__"])
        out[i] = sv[i] - (sl * dv[i] + ic)
    return out


# ---------------------------------------------------------------- empirical null
class TypeConditionalNull:
    """Right-tail empirical null of a malignancy score, conditional on cell type.

    Fitted on cancer-free reference cells only. `min_n` guards against thin strata:
    a cell type with too few reference cells falls back to the pooled null, which is
    conservative for the clustered-locus types (they inflate the pooled tail).
    """

    def __init__(self, min_n=100, residualise_depth=True):
        self.min_n = min_n
        self.residualise_depth = residualise_depth
        self.null_ = {}
        self.depth_models_ = None
        self.pooled_ = None

    def fit(self, ref, score_col="score", type_col="cell_type", depth_col="log_depth"):
        ref = ref.dropna(subset=[score_col]).copy()
        if self.residualise_depth:
            self.depth_models_ = fit_depth_model(ref, score_col, depth_col, type_col)
            ref["_adj"] = apply_depth_model(ref, self.depth_models_, score_col, depth_col, type_col)
        else:
            ref["_adj"] = ref[score_col].values
        self.pooled_ = np.sort(ref["_adj"].values)
        for t, g in ref.groupby(type_col, observed=True):
            if len(g) >= self.min_n:
                self.null_[str(t)] = np.sort(g["_adj"].values)
        self.n_types_ = len(self.null_)
        self.n_ref_ = len(ref)
        return self

    def _adjust(self, df, score_col, type_col, depth_col):
        if self.residualise_depth:
            return apply_depth_model(df, self.depth_models_, score_col, depth_col, type_col)
        return df[score_col].values

    def pvalues(self, query, score_col="score", type_col="cell_type", depth_col="log_depth",
                tail="gpd"):
        """Right-tail p-value per cell against its stratum null.

        `tail="gpd"` extrapolates the extreme tail beyond the empirical resolution
        floor of 1/n_ref; `tail="empirical"` uses the (r+1)/(n+1) estimate only.
        """
        adj = self._adjust(query, score_col, type_col, depth_col)
        tv = query[type_col].astype(str).values
        p = np.empty(len(query))
        src = np.empty(len(query), dtype=object)
        keys = np.array([t if t in self.null_ else "__pooled__" for t in tv])
        self.tail_fits_ = {}
        for k in pd.unique(keys):
            m = keys == k
            ref = self.pooled_ if k == "__pooled__" else self.null_[k]
            src[m] = "pooled" if k == "__pooled__" else "type"
            if tail == "gpd":
                pk, info = gpd_pvalues(adj[m], ref)
                self.tail_fits_[k] = info
            else:
                n = len(ref)
                r = n - np.searchsorted(ref, adj[m], side="left")
                pk = (r + 1.0) / (n + 1.0)
            p[m] = pk
        return p, adj, src


# ---------------------------------------------------------------- FDR
def bh_qvalues(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def calibrated_calls(query, null, q=0.05, score_col="score", type_col="cell_type",
                     depth_col="log_depth", by=None):
    """Attach empirical p, BH q and a calibrated malignant call to `query`.

    `by` groups the FDR correction (normally the sample): FDR is a per-sample statement,
    so pooling cells from several samples into one BH run would be wrong.
    """
    out = query.copy()
    p, adj, src = null.pvalues(query, score_col, type_col, depth_col)
    out["score_adj"] = adj
    out["p_emp"] = p
    out["null_source"] = src
    if by is None:
        out["q_emp"] = bh_qvalues(p)
    else:
        out["q_emp"] = np.nan
        for _, idx in out.groupby(by, observed=True).groups.items():
            out.loc[idx, "q_emp"] = bh_qvalues(out.loc[idx, "p_emp"].values)
    out["call_calibrated"] = (out["q_emp"] <= q).map({True: "aneuploid", False: "diploid"})
    return out


# ---------------------------------------------------------------- label-free archetypes
def reference_centroids(ref_X, ref_types, genes, min_n=50):
    """Mean log-normalised profile per reference cell type."""
    df = pd.DataFrame(ref_X, columns=genes)
    df["_t"] = np.asarray(ref_types).astype(str)
    keep = df["_t"].value_counts()
    keep = keep[keep >= min_n].index
    cen = df[df["_t"].isin(keep)].groupby("_t").mean()
    return cen


def assign_archetype(query_X, centroids, genes, min_corr=0.0):
    """Assign each query cell to its best-correlating reference cell type.

    Used for the label-free variant of the calibration: a malignant cell still resembles
    some normal archetype, so the question becomes whether its score is extreme relative
    to cancer-free cells of that same archetype.
    """
    shared = [g for g in centroids.columns if g in set(genes)]
    gi = {g: i for i, g in enumerate(genes)}
    cols = [gi[g] for g in shared]
    Q = np.asarray(query_X)[:, cols]
    C = centroids[shared].values
    Qz = (Q - Q.mean(1, keepdims=True)) / (Q.std(1, keepdims=True) + 1e-9)
    Cz = (C - C.mean(1, keepdims=True)) / (C.std(1, keepdims=True) + 1e-9)
    R = Qz @ Cz.T / Q.shape[1]
    best = R.argmax(1)
    return pd.DataFrame({"archetype": centroids.index.values[best],
                         "archetype_corr": R[np.arange(len(best)), best],
                         "archetype_margin": np.sort(R, 1)[:, -1] - np.sort(R, 1)[:, -2]})


# ---------------------------------------------------------------- evaluation helpers
def fpr_table(df, call_col, group_cols, truth_col=None):
    """False-positive rate per group. On cancer-free panels every positive is false."""
    def _f(g):
        pos = (g[call_col] == "aneuploid")
        if truth_col is not None:
            neg = g[truth_col] == "diploid"
            return pd.Series({"n": int(neg.sum()),
                              "n_fp": int((pos & neg).sum()),
                              "fpr": float((pos & neg).sum() / max(neg.sum(), 1))})
        return pd.Series({"n": len(g), "n_fp": int(pos.sum()),
                          "fpr": float(pos.mean())})
    return df.groupby(group_cols, observed=True).apply(_f, include_groups=False).reset_index()


def realized_vs_nominal(df, q_grid, truth_col="spike_truth", by=None):
    """Realized FDR and sensitivity across a grid of nominal q, for calibration curves."""
    rows = []
    for q in q_grid:
        call = df["q_emp"] <= q
        tp = int((call & (df[truth_col] == "aneuploid")).sum())
        fp = int((call & (df[truth_col] == "diploid")).sum())
        pos = int((df[truth_col] == "aneuploid").sum())
        rows.append({"q_nominal": q, "n_called": int(call.sum()), "tp": tp, "fp": fp,
                     "fdr_realized": fp / max(tp + fp, 1),
                     "sensitivity": tp / max(pos, 1)})
    return pd.DataFrame(rows)


def wilson_ci(k, n, alpha=0.05):
    """Wilson score interval — FP counts are often 0, where Wald intervals are useless."""
    if n == 0:
        return (np.nan, np.nan)
    z = stats.norm.ppf(1 - alpha / 2)
    ph = k / n
    d = 1 + z**2 / n
    c = (ph + z**2 / (2 * n)) / d
    h = z * np.sqrt(ph * (1 - ph) / n + z**2 / (4 * n**2)) / d
    return (max(0.0, c - h), min(1.0, c + h))


# --- extreme-tail p-values -------------------------------------------------
# Empirical p-values cannot go below 1/n_null, so with a few thousand reference
# cells no cell can survive Benjamini-Hochberg across a whole dataset. Fitting a
# generalized Pareto distribution to the upper tail of the null extrapolates
# beyond the sampling resolution (Knijnenburg et al. 2009, Bioinformatics 25:i161).
from scipy import stats as _st

def gpd_pvalues(x_query, x_null, tail_frac=0.10, min_exc=50):
    x_query = np.asarray(x_query, float); x_null = np.asarray(x_null, float)
    n = len(x_null)
    emp = (1.0 + (x_null[None, :] >= x_query[:, None]).sum(1)) / (n + 1.0)
    n_exc = int(round(tail_frac * n))
    if n_exc < min_exc:
        return emp, {"fit": "none", "reason": "too few null cells"}
    u = np.quantile(x_null, 1.0 - tail_frac)
    exc = x_null[x_null > u] - u
    if len(exc) < min_exc:
        return emp, {"fit": "none", "reason": "too few exceedances"}
    try:
        c, loc, scale = _st.genpareto.fit(exc, floc=0.0)
        ad = _st.anderson_ksamp([exc, _st.genpareto.rvs(c, loc, scale, size=min(2000, 10*len(exc)),
                                                        random_state=0)])
        ok = ad.pvalue > 0.05
    except Exception as e:
        return emp, {"fit": "none", "reason": str(e)[:60]}
    if not ok:
        return emp, {"fit": "rejected", "ad_p": float(ad.pvalue)}
    p = emp.copy()
    m = x_query > u
    p[m] = tail_frac * _st.genpareto.sf(x_query[m] - u, c, loc, scale)
    p = np.clip(p, 1e-300, 1.0)
    return p, {"fit": "gpd", "shape": float(c), "scale": float(scale), "u": float(u),
               "ad_p": float(ad.pvalue), "n_exc": int(len(exc))}
