"""Attribute false malignant calls to cell types and to genomic loci.

Two questions:
  (1) which cell types are called malignant in tumour-free tissue, relative to their
      abundance in the panel;
  (2) where in the genome the spurious CNV signal sits, and whether it concentrates at
      clustered loci of co-regulated gene families (immunoglobulin, HLA, histone,
      keratin, epidermal differentiation complex, ...) rather than spreading evenly.
"""
import os, glob, json
import numpy as np, pandas as pd

BIN = 10_000_000  # genomic bin size for cross-dataset profile aggregation


# ---------------------------------------------------------------- cell types
def celltype_enrichment(D, rule, mode):
    """Log2 enrichment of each cell type among called cells vs its panel share."""
    d = D[D["mode"] == mode]
    n = len(d)
    called = d[d[rule].astype(bool)]
    if not len(called):
        return pd.DataFrame()
    base = d["cell_type"].value_counts() / n
    hit = called["cell_type"].value_counts() / len(called)
    tab = pd.DataFrame({"panel_frac": base, "called_frac": hit}).fillna(0.0)
    tab["n_cells"] = d["cell_type"].value_counts()
    tab["n_called"] = called["cell_type"].value_counts().reindex(tab.index).fillna(0).astype(int)
    tab["rate"] = tab["n_called"] / tab["n_cells"]
    eps = 0.5 / n
    tab["log2_enrich"] = np.log2((tab["called_frac"] + eps) / (tab["panel_frac"] + eps))
    tab["rule"], tab["mode"] = rule, mode
    return tab.reset_index(names="cell_type")


# ---------------------------------------------------------------- genomic loci
def binned_profiles(profile_dir="calls", rule="R4_tirosh", mode="reffree", groups=None):
    """Per-dataset called-minus-uncalled CNV difference, aggregated into fixed bins."""
    rows = []
    for f in sorted(glob.glob(os.path.join(profile_dir, f"*_{mode}_profile.npz"))):
        name = os.path.basename(f).replace(f"_{mode}_profile.npz", "")
        if groups is not None and name not in groups:
            continue
        z = np.load(f, allow_pickle=False)
        if f"mean_called_{rule}" not in z or int(z[f"n_called_{rule}"][0]) < 20:
            continue
        diff = z[f"mean_called_{rule}"] - z[f"mean_uncalled_{rule}"]
        chrom = np.array([c.replace("chr", "") for c in z["chrom"].astype(str)])
        b = (z["pos"] // BIN).astype(int)
        df = pd.DataFrame({"chrom": chrom, "bin": b, "diff": diff, "dataset": name})
        df = df[df["chrom"] != ""]
        rows.append(df.groupby(["chrom", "bin"], as_index=False)["diff"].mean().assign(dataset=name))
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def family_bins(fam, bin_size=BIN):
    """Set of (chrom, bin) occupied by each clustered gene family."""
    out = {}
    for _, r in fam[fam["clustered"]].iterrows():
        c = str(r["cluster_chrom"])
        bs = int(r["cluster_start"]) // bin_size
        be = int(r["cluster_end"]) // bin_size
        out[r["family"]] = {(c, b) for b in range(bs, be + 1)}
    return out


def locus_enrichment(bp, fam, n_perm=2000, seed=0):
    """Is |called-minus-uncalled| signal larger in clustered-family bins than elsewhere?

    Permutation test on the dataset-averaged absolute difference, circularly shifting
    bins within each chromosome so autocorrelation along the genome is preserved.
    """
    m = bp.groupby(["chrom", "bin"], as_index=False)["diff"].mean()
    m["abs"] = m["diff"].abs()
    fb = family_bins(fam)
    allfb = set().union(*fb.values()) if fb else set()
    key = list(zip(m["chrom"], m["bin"]))
    m["in_family"] = [k in allfb for k in key]
    obs = m.loc[m["in_family"], "abs"].mean() - m.loc[~m["in_family"], "abs"].mean()

    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    by_chr = {c: g.sort_values("bin").reset_index(drop=True) for c, g in m.groupby("chrom")}
    for i in range(n_perm):
        parts = []
        for c, g in by_chr.items():
            s = rng.integers(0, len(g))
            gg = g.copy()
            gg["abs"] = np.roll(g["abs"].values, s)
            parts.append(gg)
        p = pd.concat(parts, ignore_index=True)
        null[i] = p.loc[p["in_family"], "abs"].mean() - p.loc[~p["in_family"], "abs"].mean()
    pval = (1 + int((null >= obs).sum())) / (1 + n_perm)

    per_fam = []
    for name, bins in fb.items():
        sel = m[[k in bins for k in key]]
        if len(sel):
            per_fam.append({"family": name, "n_bins": len(sel), "mean_abs_diff": sel["abs"].mean()})
    per_fam = pd.DataFrame(per_fam).sort_values("mean_abs_diff", ascending=False)
    bg = m.loc[~m["in_family"], "abs"].mean()
    per_fam["fold_vs_background"] = per_fam["mean_abs_diff"] / bg
    return {"obs_diff": float(obs), "p_perm": float(pval), "background_mean_abs": float(bg),
            "null_mean": float(null.mean()), "null_sd": float(null.std())}, m, per_fam
