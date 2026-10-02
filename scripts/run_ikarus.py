"""Run ikarus (supervised tumour/normal gene-signature classifier) on the panel.

ikarus takes no normal reference and emits a hard Tumor/Normal label per cell, so
there is a single mode and no user-chosen threshold. Its published signature set
and pretrained core model are used unchanged.
"""
import os, sys, json, time
sys.path.insert(0, os.path.abspath("pyshim"))
import numpy as np, pandas as pd, anndata as ad, scanpy as sc
from ikarus import classifier, data

SIG = os.path.abspath("ikarus_assets/signatures.gmt")
MODEL = os.path.abspath("ikarus_assets/core_model.joblib")

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



def run_one(path, out_dir="calls", seed=0):
    name = os.path.splitext(os.path.basename(path))[0]
    A = ad.read_h5ad(path)
    if "feature_name" in A.var:
        A.var["ensembl"] = A.var_names
        A.var_names = A.var["feature_name"].astype(str).values
        A.var_names_make_unique()
    A.var["gene_symbol"] = A.var_names.astype(str)
    A = data.preprocess_adata(A)
    t0 = time.time()
    model = classifier.Ikarus(signatures_gmt=SIG, out_dir="ikarus_out", n_neighbors=100,
                              adapt_signatures=True)
    model.load_core_model(MODEL)
    pred = model.predict(A, name, save=False)
    el = round(time.time() - t0, 1)

    ct = A.obs["cell_type"].astype(str)
    out = pd.DataFrame({"cell": A.obs_names, "dataset": name, "mode": "supervised",
                        "cell_type": ct.values,
                        "is_immune": [is_immune(c) for c in ct],
                        "donor_id": A.obs["donor_id"].astype(str).values if "donor_id" in A.obs else "NA",
                        "ikarus_pred": np.asarray(pred).ravel()})
    for c in ("core_pred", "Tumor", "Normal"):
        if hasattr(model, "results") and c in model.results.columns:
            out["ikarus_" + c] = model.results[c].values
    os.makedirs(out_dir, exist_ok=True)
    out.to_csv(f"{out_dir}/{name}_ikarus.csv", index=False)
    meta = {"dataset": name, "n_cells": int(A.n_obs), "runtime_s": el,
            "counts": {str(k): int(v) for k, v in pd.Series(out["ikarus_pred"]).value_counts().items()},
            "fpr_tumor": float((out["ikarus_pred"] == "Tumor").mean())}
    json.dump(meta, open(f"{out_dir}/{name}_ikarus.json", "w"), indent=1)
    return meta


if __name__ == "__main__":
    targets = sys.argv[1:] or sorted(
        os.path.join("panel", f) for f in os.listdir("panel") if f.endswith(".h5ad"))
    for p in targets:
        nm = os.path.splitext(os.path.basename(p))[0]
        if os.path.exists(f"calls/{nm}_ikarus.json"):
            print(nm, "SKIP", flush=True); continue
        try:
            m = run_one(p)
            print(f"{nm} OK cells={m['n_cells']} tumour={m['fpr_tumor']:.3f} {m['runtime_s']}s", flush=True)
        except Exception as e:
            import traceback; print(f"{nm} FAIL {type(e).__name__}: {e}", flush=True); traceback.print_exc()
