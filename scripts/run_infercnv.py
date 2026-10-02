"""Run the inferCNV-style CNV score (infercnvpy) on one panel dataset.

Two reference modes are exercised because both are used in practice:
  reffree  - no normal reference supplied; the baseline is the mean of all cells
             (what a practitioner gets when the tissue is assumed to contain tumour)
  refimm   - immune cells supplied as the normal reference (standard in tumour studies)

Three malignancy decision rules are applied to the same score, because the tools
publish a score but leave the cut to the user:
  R1_clusterGMM   two-component Gaussian mixture on per-cluster mean score
  R2_cellQ90      top decile of cells by score
  R3_refSD        score > mean + 3 SD of the reference cells (reference mode only)
"""
import os, sys, json, time
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad, scanpy as sc
import infercnvpy as cnv
import infercnvpy.tl._infercnv as _ic


def _serial_map(fn, *iterables, **kw):
    """Serial stand-in for tqdm's process_map: the sandbox blocks the named pipes a
    ProcessPoolExecutor needs. Identical computation, one process."""
    kw.pop("max_workers", None); kw.pop("tqdm_class", None); kw.pop("chunksize", None)
    return [fn(*args) for args in zip(*iterables)]


_ic.process_map = _serial_map
from sklearn.mixture import GaussianMixture

import re as _re
IMMUNE_PAT = ("t cell", "b cell", "nk cell", "natural killer", "macrophage", "monocyte",
              "dendritic", "mast cell", "plasma cell", "neutrophil", "lymphocyte",
              "granulocyte", "myeloid", "leukocyte", "microglial", "kupffer",
              "langerhans", "thymocyte", "basophil", "eosinophil")
# word-boundary match: plain substring matching mislabels epithelial types
# ("goblet cell" contains "t cell", "club cell" contains "b cell")
_IMMUNE_RX = _re.compile(r"\b(" + "|".join(p.replace(" ", r"\s") for p in IMMUNE_PAT) + r")\b")


def is_immune(ct):
    return bool(_IMMUNE_RX.search(str(ct).lower()))



def attach_positions(A, gene_pos):
    """Map var to chromosome/start/end; keep only genes with a placed autosomal locus."""
    gp = gene_pos.dropna(subset=["symbol"]).drop_duplicates("symbol").set_index("symbol")
    sym = A.var["feature_name"].astype(str) if "feature_name" in A.var else pd.Series(A.var_names, index=A.var_names)
    hit = sym.isin(gp.index)
    A = A[:, hit.values].copy()
    s = sym[hit.values]
    A.var["chromosome"] = ("chr" + gp.loc[s, "chrom"].astype(str)).values
    A.var["start"] = gp.loc[s, "start"].values
    A.var["end"] = gp.loc[s, "end"].values
    A = A[:, A.var["chromosome"].isin(["chr" + str(i) for i in range(1, 23)])].copy()
    order = np.lexsort((A.var["start"].values,
                        [int(c[3:]) for c in A.var["chromosome"].values]))
    return A[:, order].copy()


def normalize(A):
    A.layers["counts"] = A.X.copy()
    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)
    return A


def cell_scores(A, key="X_cnv"):
    M = A.obsm[key]
    M = M.toarray() if sp.issparse(M) else np.asarray(M)
    C = A.layers["counts"]
    tot = np.asarray(C.sum(1)).ravel() if sp.issparse(C) else np.asarray(C).sum(1)
    ng = np.asarray((C > 0).sum(1)).ravel() if sp.issparse(C) else (np.asarray(C) > 0).sum(1)
    return pd.DataFrame({"cnv_mean_abs": np.abs(M).mean(1),
                         "cnv_mean_sq": (M ** 2).mean(1),
                         "n_counts": tot,
                         "n_genes_det": ng,
                         "log_depth": np.log10(tot + 1)}, index=A.obs_names)


def cluster_labels(A, resolution=1.0, seed=0):
    sc.pp.neighbors(A, use_rep="X_cnv", key_added="cnv_neighbors", random_state=seed)
    sc.tl.leiden(A, resolution=resolution, key_added="cnv_leiden",
                 neighbors_key="cnv_neighbors", random_state=seed, flavor="igraph", n_iterations=2)
    return A.obs["cnv_leiden"].astype(str)


def rule_cluster_gmm(score, clusters, seed=0):
    """Two-component GMM on per-cluster mean score; higher-mean component = malignant."""
    cm = score.groupby(clusters).mean()
    if cm.size < 2:
        return pd.Series(False, index=score.index), {}
    g = GaussianMixture(2, random_state=seed, n_init=5).fit(cm.values.reshape(-1, 1))
    hi = int(np.argmax(g.means_.ravel()))
    mal_clusters = set(cm.index[g.predict(cm.values.reshape(-1, 1)) == hi])
    return clusters.isin(mal_clusters), {"malignant_clusters": sorted(mal_clusters),
                                         "cluster_means": {str(k): float(v) for k, v in cm.items()}}


def rule_twogroup(M, score, seed=0):
    """Forced two-group split: Ward linkage on the CNV profile cut into two groups,
    the higher-mean-score group called malignant. This is the shape of the decision
    CopyKAT and similar tools make, and it cannot return 'no tumour present'."""
    from sklearn.decomposition import PCA
    from scipy.cluster.hierarchy import linkage, fcluster
    M = M.toarray() if sp.issparse(M) else np.asarray(M)
    P = PCA(n_components=min(30, M.shape[1], M.shape[0] - 1), random_state=seed).fit_transform(M)
    lab = fcluster(linkage(P, method="ward"), 2, criterion="maxclust")
    means = {g: score.values[lab == g].mean() for g in np.unique(lab)}
    hi = max(means, key=means.get)
    return lab == hi, {str(k): float(v) for k, v in means.items()}


def rule_topcluster(score, clusters):
    """The single leiden cluster with the highest mean CNV score is called malignant --
    the heatmap-inspection heuristic, also a forced choice."""
    cm = score.groupby(clusters).mean()
    top = cm.idxmax()
    return (clusters == top).values, {"top_cluster": str(top), "top_mean": float(cm.max())}


def rule_tirosh(M, score, ref_mask, t_cor=0.4, top_frac=0.05):
    """Tirosh/Patel two-criterion rule: a cell is malignant if its CNV signal exceeds
    mean+2SD of the reference set AND its CNV profile correlates (> t_cor) with the
    mean profile of the top-signal cells. Where no reference is supplied, the
    bottom quintile of cells by signal serves as the internal reference, which is
    what is done in practice when the tissue has no annotated normal compartment."""
    M = M.toarray() if sp.issparse(M) else np.asarray(M)
    if ref_mask is None or ref_mask.sum() < 50:
        ref_mask = (score <= np.quantile(score, 0.20)).values
    t_sig = score[ref_mask].mean() + 2 * score[ref_mask].std()
    k = max(10, int(round(top_frac * M.shape[0])))
    prof = M[np.argsort(-score.values)[:k]].mean(0)
    Mc = M - M.mean(1, keepdims=True)
    pc = prof - prof.mean()
    denom = np.sqrt((Mc ** 2).sum(1) * (pc ** 2).sum())
    corr = np.where(denom > 0, (Mc @ pc) / np.where(denom > 0, denom, 1), 0.0)
    return (score.values > t_sig) & (corr > t_cor), corr, float(t_sig)


def window_coords(A):
    """Map X_cnv columns back to genomic coordinates. infercnvpy stores the first
    column index of each chromosome in .uns['cnv']['chr_pos']; within a chromosome the
    columns are evenly spaced windows over that chromosome's genes in order."""
    chr_pos = A.uns["cnv"]["chr_pos"]
    ncol = A.obsm["X_cnv"].shape[1]
    order = sorted(chr_pos.items(), key=lambda kv: kv[1])
    chrom_out = np.empty(ncol, dtype=object)
    pos_out = np.zeros(ncol, dtype=np.int64)
    for i, (c, s) in enumerate(order):
        e = order[i + 1][1] if i + 1 < len(order) else ncol
        g = A.var.loc[A.var["chromosome"] == c, "start"].values
        n = e - s
        if n <= 0 or len(g) == 0:
            continue
        idx = np.clip(((np.arange(n) + 0.5) * len(g) / n).astype(int), 0, len(g) - 1)
        chrom_out[s:e] = c
        pos_out[s:e] = g[idx]
    return chrom_out, pos_out


def save_profiles(A, frame, name, mode, out_dir):
    """Mean CNV profile of called vs uncalled cells, for the locus attribution."""
    M = A.obsm["X_cnv"]
    M = M.toarray() if sp.issparse(M) else np.asarray(M)
    chrom, pos = window_coords(A)
    d = {"chrom": np.asarray(chrom, dtype=str), "pos": pos, "mean_all": M.mean(0)}
    for rule in ("R1_twogroup", "R4_tirosh"):
        m = frame[rule].values.astype(bool)
        d[f"mean_called_{rule}"] = M[m].mean(0) if m.sum() else np.zeros(M.shape[1])
        d[f"mean_uncalled_{rule}"] = M[~m].mean(0) if (~m).sum() else np.zeros(M.shape[1])
        d[f"n_called_{rule}"] = np.array([m.sum()])
    np.savez_compressed(f"{out_dir}/{name}_{mode}_profile.npz", **d)


def run_one(path, gene_pos, out_dir, seed=0):
    name = os.path.splitext(os.path.basename(path))[0]
    A0 = ad.read_h5ad(path)
    A0 = attach_positions(A0, gene_pos)
    A0 = normalize(A0)
    ct = A0.obs["cell_type"].astype(str)
    A0.obs["is_immune"] = [is_immune(c) for c in ct]

    out = {"dataset": name, "n_cells": int(A0.n_obs), "n_genes_placed": int(A0.n_vars),
           "n_immune": int(A0.obs["is_immune"].sum()), "modes": {}}
    frames = []
    for mode in ("reffree", "refimm"):
        if mode == "refimm" and A0.obs["is_immune"].sum() < 50:
            out["modes"][mode] = {"skipped": "fewer than 50 immune cells"}
            continue
        A = A0.copy()
        t0 = time.time()
        if mode == "reffree":
            cnv.tl.infercnv(A, window_size=100, step=10, lfc_clip=3, dynamic_threshold=1.5,
                            n_jobs=1)
        else:
            cnv.tl.infercnv(A, reference_key="is_immune", reference_cat=[True],
                            window_size=100, step=10, lfc_clip=3, dynamic_threshold=1.5,
                            n_jobs=1)
        sc_ = cell_scores(A)
        clusters = cluster_labels(A, seed=seed)
        s = sc_["cnv_mean_sq"]

        r1, gmeans = rule_twogroup(A.obsm["X_cnv"], s, seed=seed)
        meta = {"twogroup_means": gmeans}
        r1b, tmeta = rule_topcluster(s, clusters)
        meta.update(tmeta)
        r2 = s >= s.quantile(0.90)
        if mode == "refimm":
            ref = s[A.obs["is_immune"].values]
            thr = ref.mean() + 3 * ref.std()
            r3 = s > thr
            meta["refSD_threshold"] = float(thr)
            ref_mask = A.obs["is_immune"].values
        else:
            r3 = pd.Series(np.nan, index=s.index)
            ref_mask = None
        r4, corr, t_sig = rule_tirosh(A.obsm["X_cnv"], s, ref_mask)
        meta["tirosh_signal_threshold"] = t_sig

        f = pd.DataFrame({"cell": s.index, "dataset": name, "mode": mode,
                          "cell_type": ct.values, "is_immune": A0.obs["is_immune"].values,
                          "donor_id": A0.obs["donor_id"].astype(str).values if "donor_id" in A0.obs else "NA",
                          "cnv_mean_abs": sc_["cnv_mean_abs"].values, "cnv_mean_sq": s.values,
                          "cnv_cluster": clusters.values,
                          "cnv_corr_top": corr,
                          "R1_twogroup": r1, "R1b_topcluster": r1b, "R2_cellQ90": r2.values,
                          "R3_refSD": r3.values, "R4_tirosh": r4})
        frames.append(f)
        os.makedirs(out_dir, exist_ok=True)
        save_profiles(A, f, name, mode, out_dir)
        meta["runtime_s"] = round(time.time() - t0, 1)
        meta["fpr"] = {r: float(np.nanmean(f[r].astype(float))) for r in
                       ("R1_twogroup", "R1b_topcluster", "R2_cellQ90", "R3_refSD", "R4_tirosh")}
        out["modes"][mode] = meta
        del A

    os.makedirs(out_dir, exist_ok=True)
    pd.concat(frames, ignore_index=True).to_parquet(f"{out_dir}/{name}_infercnv.parquet", index=False)
    json.dump(out, open(f"{out_dir}/{name}_infercnv.json", "w"), indent=1)
    return out


if __name__ == "__main__":
    gene_pos = pd.read_csv("gene_positions.csv")
    targets = sys.argv[1:] or sorted(
        os.path.join("panel", f) for f in os.listdir("panel") if f.endswith(".h5ad"))
    for p in targets:
        nm = os.path.splitext(os.path.basename(p))[0]
        if os.path.exists(f"calls/{nm}_infercnv.json"):
            print(f"{nm} SKIP", flush=True); continue
        try:
            o = run_one(p, gene_pos, "calls")
            print(f"{nm} OK cells={o['n_cells']} " +
                  " ".join(f"{m}:{v.get('runtime_s','-')}s" for m, v in o["modes"].items()), flush=True)
        except Exception as e:
            import traceback
            print(f"{nm} FAIL {type(e).__name__}: {e}", flush=True)
            traceback.print_exc()
