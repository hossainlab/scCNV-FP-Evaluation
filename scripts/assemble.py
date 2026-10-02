"""Assemble the null-control panel: download each CELLxGENE h5ad, subsample to a
compact raw-count dataset, then delete the raw file to stay inside the disk budget."""
import json, os, sys, time, traceback
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dl import fetch
from prep import read_obs, build_panel

# The atlases are already author-QC'd. A tight mitochondrial cap would remove
# precisely the mitochondria-rich epithelial types most at risk of a false
# aneuploidy call, biasing the audit toward a null result, so the cap is lenient.
MT_CAP = 50.0
panel = json.load(open("config/panel_def.json"))
os.makedirs("raw", exist_ok=True)
os.makedirs("panel", exist_ok=True)
LOG = "assemble_log.txt"

def log(msg):
    with open(LOG, "a") as f:
        f.write(f"[{time.strftime('%H:%M:%S')}] {msg}\n")

summary = {}
for p in panel:
    name, out = p["name"], f"panel/{p['name']}.h5ad"
    if os.path.exists(out):
        log(f"{name} SKIP (exists)"); continue
    rawf = f"raw/{name}.h5ad"
    try:
        t0 = time.time()
        fetch(p["url"], rawf)
        dl = time.time() - t0
        # pick dominant assay when not pre-specified, so one protocol per dataset
        assay = p["assay"]
        if assay is None:
            o = read_obs(rawf, ["assay"])
            assay = [o["assay"].value_counts().index[0]]
        A = build_panel(rawf, out, name, assay_keep=assay, per_type=200, cap=3000,
                        min_counts=500, min_genes=250, max_pct_mt=MT_CAP)
        ct = A.obs["cell_type"].value_counts()
        summary[name] = {"group": p["group"], "assay": assay,
                         "n_cells": int(A.n_obs), "n_genes": int(A.n_vars),
                         "n_cell_types": int(ct.size),
                         "n_donors": int(A.obs["donor_id"].nunique()) if "donor_id" in A.obs else None,
                         "diseases": sorted(map(str, A.obs["disease"].unique())) if "disease" in A.obs else None,
                         "cell_types": {str(k): int(v) for k, v in ct.items()},
                         "dl_s": round(dl, 1)}
        log(f"{name} OK cells={A.n_obs} types={ct.size} dl={dl:.0f}s")
        json.dump(summary, open("results/panel_summary.json", "w"), indent=1)
    except Exception as e:
        log(f"{name} FAIL {type(e).__name__}: {e}")
        log(traceback.format_exc()[:1500])
    finally:
        if os.path.exists(rawf):
            os.remove(rawf)
log("DONE")
