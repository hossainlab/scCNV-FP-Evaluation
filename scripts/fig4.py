"""Figure 4 — a calibrated decision procedure restores nominal false-positive control."""
import numpy as np, pandas as pd, matplotlib.pyplot as plt

from figstyle import apply_figure_style, panel_letter, META_GREY
apply_figure_style(sizes=(8, 7, 6))
SW = pd.read_csv("results/T5c_alpha_sweep.csv")
CM = pd.read_csv("results/T10_condition_matched.csv")

C_RULE, C_CAL, C_ARCH = "#6E7B8B", "#3C6E9F", "#7BA05B"
PRETTY = {"TS_Bone_Marrow": "bone marrow", "TS_Kidney": "kidney",
          "TS_Large_Intestine": "large intestine", "TS_Liver": "liver", "TS_Lung": "lung",
          "TS_Pancreas": "pancreas", "TS_Skin": "skin", "TS_Small_Intestine": "small intestine",
          "TS_Spleen": "spleen", "TS_Stomach": "stomach", "UC_colon": "colon (colitis)",
          "PSO_skin": "skin (psoriasis)", "CD_epithelial": "colon (Crohn's, epi.)",
          "CD_immune": "colon (Crohn's, imm.)"}

fig = plt.figure(figsize=(7.2, 5.2))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.68, wspace=0.34)
axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, :]), fig.add_subplot(gs[0, 1])]

# ---- (a) realized vs nominal on held-out tissues
ax = axes[0]
for src, cc, lab, mk in [("labels", C_CAL, "compared within cell type", "o"),
                         ("archetype", C_ARCH, "compared within cell group\n(no labels needed)", "s")]:
    T = SW[(SW.source == src) & (SW.dataset != "PSO_skin")]
    g = T.groupby("alpha")["realized"].agg(["mean", "std"]).reset_index()
    ax.plot(100 * g.alpha, 100 * g["mean"], "-", marker=mk, ms=3.4, lw=1.2, color=cc, label=lab)
    ax.fill_between(100 * g.alpha, 100 * (g["mean"] - g["std"]), 100 * (g["mean"] + g["std"]),
                    color=cc, alpha=.16, lw=0)
ax.plot([0, 10], [0, 10], color=META_GREY, lw=.8, ls=(0, (3, 2)))
ax.text(9.6, 9.9, "asked = got", fontsize=6, color=META_GREY, ha="right", va="bottom", rotation=38)
ax.text(0.4, 12.6, "13 datasets (psoriasis set: see b, c)",
        fontsize=5.8, color=META_GREY, va="top")
ax.set_xlabel("error rate we asked for, $\\alpha$ (%)")
ax.set_ylabel("error rate we got (%)")
ax.legend(frameon=False, loc="upper left", handlelength=1.4, bbox_to_anchor=(-0.02, 0.82))
ax.set_xlim(0, 10.5); ax.set_ylim(0, 13)
panel_letter(ax, "a")

# ---- (b) per-tissue: published rule vs calibrated at alpha = 5%
ax = axes[1]
a5 = SW[(SW.source == "labels") & np.isclose(SW.alpha, 0.05)].set_index("dataset")["realized"] * 100
T1 = pd.read_csv("results/T1_infercnv_fpr.csv")
rule = (T1[T1["mode"] == "refimm"].set_index("dataset")["R4_tirosh"]
        .combine_first(T1[T1["mode"] == "reffree"].set_index("dataset")["R4_tirosh"]))
ds = [d for d in PRETTY if d in a5.index]
x = np.arange(len(ds))
ax.vlines(x, [a5[d] for d in ds], [rule[d] for d in ds], color=META_GREY, lw=.8, zorder=1)
ax.scatter(x, [rule[d] for d in ds], s=26, color=C_RULE, zorder=3, label="published cutoff rule (Tirosh)")
ax.scatter(x, [a5[d] for d in ds], s=26, color=C_CAL, zorder=3, label="calibrated, asked for 5%")
ax.axhline(5, color=C_CAL, lw=.8, ls=(0, (2, 2)), zorder=0)
ax.set_xticks(x); ax.set_xticklabels([PRETTY[d] for d in ds], rotation=35, ha="right")
ax.set_ylabel("wrongly called cells (%)")
ax.legend(frameon=False, loc="upper left", handletextpad=.3, ncol=2)
ax.set_ylim(-2, 34); ax.margins(x=.02)
ax.annotate("reference from other tissues\nfails when the protocol differs", (ds.index("PSO_skin"), 26.7),
            xytext=(-8, -12), textcoords="offset points", fontsize=5.8, color=C_CAL, ha="right", va="top")
panel_letter(ax, "b")

# ---- (c) inflamed tissue needs a condition-matched null
ax = axes[2]
g = CM.set_index(["dataset", "null_source"])["realized"].unstack() * 100
g = g.sort_values("healthy", ascending=False)
lab = [PRETTY.get(d, d) for d in g.index]
xs = np.arange(len(g)); w = .36
ax.bar(xs - w / 2, g["healthy"], w, color=C_RULE, label="reference: other tissues")
ax.bar(xs + w / 2, g["matched"], w, color=C_CAL, label="reference: other donors,\nsame dataset")
ax.axhline(5, color=META_GREY, lw=.8, ls=(0, (2, 2)))
ax.text(len(g) - .4, 5.8, "$\\alpha$ = 5%", fontsize=6, color=META_GREY, ha="right")
ax.set_ylim(0, 32)
ax.set_xticks(xs); ax.set_xticklabels(lab, rotation=32, ha="right", fontsize=5.8)
ax.set_ylabel("wrongly called cells (%)")
ax.legend(frameon=False, loc="upper right", handlelength=1.0, labelspacing=.6)
panel_letter(ax, "c")

# align the three panel letters on a common figure-x, using panel (a) as anchor
p0 = axes[0].get_position()
target_x = p0.x0 - 0.18 * p0.width
for ax_, ch in zip(axes[:2], "ab"):
    p = ax_.get_position()
    for tx in ax_.texts:
        if tx.get_text() == ch and tx.get_fontweight() in ("bold", 700):
            tx.set_x((target_x - p.x0) / p.width)

fig.savefig("figures/Fig4_calibration.png", dpi=330, bbox_inches="tight")
fig.savefig("figures/Fig4_calibration.pdf", bbox_inches="tight")
print("saved")
