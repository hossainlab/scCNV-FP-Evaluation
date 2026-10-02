"""Build in-silico aneuploid positive controls from held-out healthy tissues."""
import os, json, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
from spikein import build_spikein

BACKGROUNDS = ["TS_Lung", "TS_Liver", "TS_Spleen"]
BURDENS = (1, 2, 3, 5, 8, 12, 18)
N_PER = 120


def gene_frame(A, gene_pos):
    gp = gene_pos.dropna(subset=["symbol"]).drop_duplicates("symbol").set_index("symbol")
    sym = A.var["feature_name"].astype(str).values if "feature_name" in A.var else np.array(A.var_names)
    chrom = np.array(["NA"] * len(sym), dtype=object)
    mid = np.zeros(len(sym), dtype=np.int64)
    hit = pd.Index(sym).isin(gp.index)
    chrom[hit] = gp.loc[sym[hit], "chrom"].astype(str).values
    mid[hit] = gp.loc[sym[hit], "mid"].values
    return pd.DataFrame({"symbol": sym, "chrom": chrom, "mid": mid}), int(hit.sum())


if __name__ == "__main__":
    gene_pos = pd.read_csv("gene_positions.csv")
    os.makedirs("spike", exist_ok=True)
    summary = {}
    for i, bg in enumerate(BACKGROUNDS):
        rng = np.random.default_rng(100 + i)
        A = ad.read_h5ad(f"panel/{bg}.h5ad")
        gdf, n_placed = gene_frame(A, gene_pos)
        X, obs, truth = build_spikein(A, gdf, rng, burdens=BURDENS, n_per_burden=N_PER)
        B = ad.AnnData(X=X, obs=obs, var=A.var.copy())
        name = f"SPK_{bg}"
        B.write_h5ad(f"spike/{name}.h5ad", compression="gzip")
        json.dump(truth, open(f"spike/{name}_truth.json", "w"), indent=1)
        summary[name] = {"background": bg, "n_cells": int(B.n_obs), "n_genes_placed": n_placed,
                         "n_aneuploid": int((obs["spike_truth"] == "aneuploid").sum()),
                         "n_diploid": int((obs["spike_truth"] == "diploid").sum()),
                         "frac_genome_by_burden": {str(t["burden"]): round(t["frac_genome"], 4) for t in truth}}
        print(name, summary[name]["n_cells"], "cells,",
              summary[name]["n_aneuploid"], "aneuploid,",
              summary[name]["frac_genome_by_burden"], flush=True)
    json.dump(summary, open("spike/spike_summary.json", "w"), indent=1)
