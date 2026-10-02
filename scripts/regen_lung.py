"""Regenerate the per-cell inferCNV scores for TS_Lung (Fig 1c input).

The per-cell score table was lost with the workspace; every other figure input
survives in the supplement archive. The panel build is seeded (seed=0) and the
scoring is deterministic, so this should reproduce the published TS_Lung rate
(R3_refSD = 6.26% under the immune-reference mode) exactly. That equality is
checked at the end and is the gate on using the result.
"""
import os, sys, json, time, urllib.request, urllib.parse
import numpy as np, pandas as pd

sys.path.insert(0, ".")
URL = "https://datasets.cellxgene.cziscience.com/fddf43a3-b7c8-49e9-ac0a-5e24b390849f.h5ad"
RAW = "raw/TS_Lung.h5ad"

# ---- 1. gene coordinate table -------------------------------------------------
if not os.path.exists("gene_positions.csv"):
    if not os.path.exists("gene_pos_biomart.tsv"):
        xml = ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
               '<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1"'
               ' count="" datasetConfigVersion="0.6">'
               '<Dataset name="hsapiens_gene_ensembl" interface="default">'
               '<Attribute name="ensembl_gene_id"/><Attribute name="external_gene_name"/>'
               '<Attribute name="chromosome_name"/><Attribute name="start_position"/>'
               '<Attribute name="end_position"/><Attribute name="gene_biotype"/>'
               '</Dataset></Query>')
        u = "https://www.ensembl.org/biomart/martservice?" + urllib.parse.urlencode({"query": xml})
        with urllib.request.urlopen(u, timeout=900) as r:
            txt = r.read().decode()
        open("gene_pos_biomart.tsv", "w", newline="").write(txt)
        print("biomart rows", txt.count("\n") - 1, flush=True)
    import loci
    loci.build()
print("gene table", len(pd.read_csv("gene_positions.csv")), flush=True)

# ---- 2. raw dataset -----------------------------------------------------------
os.makedirs("raw", exist_ok=True)
EXPECTED = 3196847981          # size returned by the store on the completed transfer
for attempt in range(8 if not (os.path.exists(RAW) and os.path.getsize(RAW) == EXPECTED) else 0):
    have = os.path.getsize(RAW) if os.path.exists(RAW) else 0
    try:
        req = urllib.request.Request(URL, headers={"Range": f"bytes={have}-"} if have else {})
        with urllib.request.urlopen(req, timeout=300) as r:
            if have and r.status != 206:      # server ignored the range: start over
                have = 0
            with open(RAW, "ab" if have else "wb") as fh:
                while True:
                    blk = r.read(1 << 22)
                    if not blk:
                        break
                    fh.write(blk)
        break
    except Exception as e:
        print("retry", attempt, type(e).__name__, str(e)[:90], flush=True)
        time.sleep(5)
print("raw bytes", os.path.getsize(RAW), flush=True)

# ---- 3. panel subsample (same parameters as assemble.py) ----------------------
from prep import build_panel
os.makedirs("panel", exist_ok=True)
A = build_panel(RAW, "panel/TS_Lung.h5ad", "TS_Lung", assay_keep=["10x 3' v3"], per_type=200, cap=3000,
                min_counts=500, min_genes=250, max_pct_mt=50.0)
print("panel cells", A.n_obs, "genes", A.n_vars, flush=True)
del A

# ---- 4. CNV scoring -----------------------------------------------------------
import run_infercnv as R
o = R.run_one("panel/TS_Lung.h5ad", pd.read_csv("gene_positions.csv"), "calls")
print("modes", json.dumps({m: v.get("fpr") for m, v in o["modes"].items()}), flush=True)

# ---- 5. reproduction check ----------------------------------------------------
d = pd.read_parquet("calls/TS_Lung_infercnv.parquet")
if "dataset" not in d.columns:
    d["dataset"] = "TS_Lung"
d.to_parquet("calls/all_infercnv_cells.parquet", index=False)
sub = d[(d.dataset == "TS_Lung") & (d["mode"] == "refimm")]
got = 100 * sub.R3_refSD.astype(float).mean()
T1 = pd.read_csv("results/T1_infercnv_fpr.csv")
want = float(T1[(T1.dataset == "TS_Lung") & (T1["mode"] == "refimm")].R3_refSD.iloc[0])
print(f"CHECK R3_refSD regenerated={got:.2f}% published={want:.2f}% delta={got-want:+.2f}", flush=True)
print("cells", len(sub), flush=True)
