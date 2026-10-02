"""Export each panel dataset as a Matrix-Market triplet that R can read without
Bioconductor: genes x cells raw counts, gene symbols, cell barcodes, and metadata."""
import os, sys, numpy as np, pandas as pd, scipy.sparse as sp, scipy.io as sio
import anndata as ad

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



def export(path, out_dir, min_cells=3):
    name = os.path.splitext(os.path.basename(path))[0]
    d = os.path.join(out_dir, name)
    os.makedirs(d, exist_ok=True)
    A = ad.read_h5ad(path)
    sym = A.var["feature_name"].astype(str).values if "feature_name" in A.var else np.array(A.var_names)
    X = A.X.tocsc() if sp.issparse(A.X) else sp.csc_matrix(A.X)
    keep = np.asarray((X > 0).sum(0)).ravel() >= min_cells
    # collapse duplicate symbols to the more highly detected copy
    order = np.argsort(-np.asarray(X.sum(0)).ravel())
    seen, uniq = set(), np.zeros(X.shape[1], bool)
    for i in order:
        if keep[i] and sym[i] not in seen and sym[i] not in ("nan", "None", ""):
            seen.add(sym[i]); uniq[i] = True
    X = X[:, uniq]; sym = sym[uniq]
    sio.mmwrite(os.path.join(d, "counts.mtx"), sp.csr_matrix(X.T).astype(np.int32), field="integer")
    pd.Series(sym).to_csv(os.path.join(d, "genes.txt"), index=False, header=False)
    pd.Series(A.obs_names).to_csv(os.path.join(d, "cells.txt"), index=False, header=False)
    pd.DataFrame({"cell": A.obs_names,
                  "cell_type": A.obs["cell_type"].astype(str).values,
                  "donor_id": A.obs["donor_id"].astype(str).values if "donor_id" in A.obs else "NA",
                  "is_immune": [is_immune(c) for c in A.obs["cell_type"].astype(str)]}
                 ).to_csv(os.path.join(d, "meta.csv"), index=False)
    return name, X.shape[0], X.shape[1]


if __name__ == "__main__":
    out = "rmat"
    targets = sys.argv[1:] or sorted(os.path.join("panel", f) for f in os.listdir("panel") if f.endswith(".h5ad"))
    for p in targets:
        print("%s cells=%d genes=%d" % export(p, out), flush=True)
