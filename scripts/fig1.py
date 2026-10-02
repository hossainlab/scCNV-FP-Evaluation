"""Figure 1 — malignant-cell callers assign tumour labels to tumour-free tissue."""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from figstyle import apply_figure_style, panel_letter, META_GREY
apply_figure_style(sizes=(8, 7, 6))

T1 = pd.read_csv("results/T1_infercnv_fpr.csv")
T2 = pd.read_csv("results/T2_ikarus_fpr.csv")
T2b = pd.read_csv("results/T2b_copykat_fpr.csv")
T0 = pd.read_csv("results/T0_datasets.csv")
# panel (c) uses only TS_Lung; fall back to the shipped per-dataset table
import os
D = pd.read_parquet("calls/all_infercnv_cells.parquet" if os.path.exists("calls/all_infercnv_cells.parquet")
                    else "calls/TS_Lung_infercnv.parquet")

PRETTY = {"TS_Bone_Marrow": "bone marrow", "TS_Kidney": "kidney",
          "TS_Large_Intestine": "large intestine", "TS_Liver": "liver", "TS_Lung": "lung",
          "TS_Pancreas": "pancreas", "TS_Skin": "skin", "TS_Small_Intestine": "small intestine",
          "TS_Spleen": "spleen", "TS_Stomach": "stomach",
          "UC_colon": "colon (colitis)", "PSO_skin": "skin (psoriasis)",
          "CD_epithelial": "colon (Crohn's, epithelium)", "CD_immune": "colon (Crohn's, immune)"}
C_HEALTHY, C_INFLAMED = "#3C6E9F", "#C1662F"
C_ACC = "#8A4F9E"
C_CK = "#5C8A4A"

fig = plt.figure(figsize=(7.2, 5.6))
gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0], width_ratios=[1.15, 1.0],
                      hspace=0.72, wspace=0.34)

# ---- (a) per-dataset false-positive rate, best-practice rule and the classifier
ax = fig.add_subplot(gs[0, :])
t = T1[(T1["mode"] == "refimm")].set_index("dataset")["R4_tirosh"]
t_free = T1[(T1["mode"] == "reffree")].set_index("dataset")["R4_tirosh"]
t = t.combine_first(t_free)
k = T2.set_index("dataset")["tumour_pct"]
order = [d for d in PRETTY if d in t.index]
x = np.arange(len(order))
grp = {d: ("inflamed" if d.startswith(("UC_", "PSO_", "CD_")) else "healthy") for d in order}
cols = [C_INFLAMED if grp[d] == "inflamed" else C_HEALTHY for d in order]
ax.vlines(x, 0, [t[d] for d in order], color=cols, lw=1.0, alpha=.55)
ax.scatter(x, [t[d] for d in order], s=26, c=cols, zorder=3, label="inferCNV (copy-number method)")
ck = T2b.set_index("dataset")["aneuploid_pct"]
ax.scatter(x, [k[d] for d in order], s=26, facecolors="none", edgecolors=C_ACC, lw=1.1,
           zorder=3, label="ikarus (trained classifier)")
ax.scatter(x, [ck[d] for d in order], s=22, marker="s", facecolors="none", edgecolors=C_CK,
           lw=1.1, zorder=3, label="CopyKAT (copy-number method)")
ax.set_xticks(x); ax.set_xticklabels([PRETTY[d] for d in order], rotation=38, ha="right")
ax.set_ylabel("cells called malignant (%)")
ax.axhline(0, color=META_GREY, lw=.6)
ax.set_ylim(-4, 104); ax.margins(x=0.02)
n_heal = sum(1 for d in order if grp[d] == "healthy")
ax.axvspan(n_heal - .5, len(order) - .5, color=C_INFLAMED, alpha=.07, zorder=0, lw=0)
ax.text(len(order) - .6, 96, "inflamed, tumour-free", color=C_INFLAMED, fontsize=6, va="top", ha="right")
ax.text(-.15, 96, "healthy", color=C_HEALTHY, fontsize=6, va="top", ha="left")
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.88), handletextpad=.4, ncol=3, columnspacing=1.2)
panel_letter(ax, "a")

# ---- (b) pooled rate by decision rule
ax = fig.add_subplot(gs[1, 0])
RULES = ["R1_twogroup", "R1b_topcluster", "R3_refSD", "R4_tirosh"]
RLAB = {"R1_twogroup": "forced 2-group split", "R1b_topcluster": "highest-scoring cluster",
        "R3_refSD": "cutoff from normal cells", "R4_tirosh": "2-criterion rule (Tirosh)"}
# cell-weighted pooling over the panel, from the shipped per-dataset tables
M = T1.merge(T0[["dataset", "n_cells"]], on="dataset", how="left")
assert M.n_cells.notna().all()
vals = []
for r in RULES:
    for m, lab in [("reffree", "no normal reference"), ("refimm", "immune cells as reference")]:
        s = M[M["mode"] == m][[r, "n_cells"]].dropna()
        if len(s):
            vals.append({"rule": RLAB[r], "mode": lab,
                         "pct": float((s[r] * s.n_cells).sum() / s.n_cells.sum())})
V = pd.DataFrame(vals)
w, xs = 0.36, np.arange(len(RULES))
for i, (m, cc) in enumerate([("no normal reference", "#8FA9C4"), ("immune cells as reference", C_HEALTHY)]):
    sub = V[V["mode"] == m].set_index("rule").reindex([RLAB[r] for r in RULES])["pct"]
    ax.bar(xs + (i - .5) * w, sub.values, w, color=cc, label=m, zorder=2)
for i, r in enumerate(RULES):
    if M[M["mode"] == "reffree"][r].isna().all():
        ax.text(xs[i] - .5 * w, 1.5, "n.d.", ha="center", fontsize=6, color=META_GREY, rotation=90)
ax.set_xticks(xs); ax.set_xticklabels([RLAB[r] for r in RULES], rotation=30, ha="right")
ax.set_ylabel("cells wrongly called\nmalignant, all datasets (%)")
ax.legend(frameon=False, loc="upper right", handlelength=1.0)
panel_letter(ax, "b")

# ---- (c) the score being thresholded is unimodal
ax = fig.add_subplot(gs[1, 1])
d = D[(D.dataset == "TS_Lung") & (D["mode"] == "refimm")]
v = np.log10(d.cnv_mean_sq.values + 1e-8)
ax.hist(v, bins=60, color="#C8D3DE", edgecolor="none")
pos = d.R3_refSD.astype(bool).values
cut = np.log10(d.cnv_mean_sq.values[pos].min() + 1e-8) if pos.any() else np.nan
ax.axvline(cut, color=C_ACC, lw=1.0, ls=(0, (4, 2)))
ymax = ax.get_ylim()[1]
ax.annotate("cutoff from\nnormal cells", (cut - .04, ymax * .99), fontsize=6, color=C_ACC,
            ha="right", va="top")
ax.annotate("%.1f%% of cells\nfall above it" % (100 * pos.mean()), (cut + .04, ymax * .99),
            fontsize=6, color=META_GREY, ha="left", va="top")
ax.set_xlabel("copy-number score per cell (log$_{10}$)")
ax.set_ylabel("cells")
panel_letter(ax, "c")

# put each panel letter just outside that panel's own labels, with (a) and (b)
# — the two panels sharing the left margin — aligned on one vertical line
letters = {}
for ax_, ch in zip(fig.axes, "abc"):
    for txt in ax_.texts:
        if txt.get_text() == ch and txt.get_fontweight() in ("bold", 700):
            letters[ch] = (ax_, txt)
for _, txt in letters.values():
    txt.set_visible(False)
fig.canvas.draw()
rend, inv = fig.canvas.get_renderer(), fig.transFigure.inverted()
tb = {ch: a_.get_tightbbox(rend).transformed(inv) for ch, (a_, _) in letters.items()}
GAP = 0.014
left_col = min(tb["a"].x0, tb["b"].x0) - GAP
for ch, (ax_, txt) in letters.items():
    x = left_col if ch in ("a", "b") else tb[ch].x0 - GAP
    p = ax_.get_position()
    txt.set_x((x - p.x0) / p.width)
    txt.set_visible(True)

fig.savefig("figures/Fig1_false_positives.png", dpi=330, bbox_inches="tight")
fig.savefig("figures/Fig1_false_positives.pdf", bbox_inches="tight")
print("saved")
