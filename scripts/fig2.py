"""Figure 2 — in-silico aneuploid spike-ins: what the callers can actually detect."""
import numpy as np, pandas as pd, matplotlib.pyplot as plt

from figstyle import apply_figure_style, panel_letter, META_GREY
apply_figure_style(sizes=(8, 7, 6))
AUC = pd.read_csv("results/T3_spike_auc.csv")
OP = pd.read_csv("results/T4_matched_fpr.csv")
KS_ = pd.read_csv("results/T3b_ikarus_burden.csv")
CB = pd.read_csv("results/T3c_copykat_burden.csv")

TIS = {"SPK_TS_Lung": "lung", "SPK_TS_Liver": "liver", "SPK_TS_Spleen": "spleen"}
COL = {"SPK_TS_Lung": "#3C6E9F", "SPK_TS_Liver": "#7BA05B", "SPK_TS_Spleen": "#B0567A"}
C_SCORE, C_RULE = "#8A4F9E", "#6E7B8B"

fig = plt.figure(figsize=(7.2, 4.5))
gs = fig.add_gridspec(2, 2, height_ratios=[1, .78], hspace=.58, wspace=.30)
axa, axb, axc = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])

# ---- (a) discrimination rises with karyotype burden but stays modest
ax = axa
a = AUC[AUC["mode"] == "reffree"]
for ds, g in a.groupby("dataset"):
    g = g.sort_values("burden")
    ax.plot(g.burden, g.auc_meansq, "-o", ms=3.2, lw=1.3, color=COL[ds], label=TIS[ds])
ax.axhline(0.5, color=META_GREY, lw=.7, ls=(0, (3, 2)))
ax.text(20.7, 0.513, "chance", fontsize=6, color=META_GREY, ha="right", va="bottom")
ax.set_xlabel("arm changes added per cell")
ax.set_ylabel("altered cells separated\nfrom normal (AUC)")
ax.set_ylim(0.45, 0.95); ax.set_xlim(0, 21); ax.set_xticks([1, 5, 10, 15, 18])
ax.legend(frameon=False, loc="upper left", handlelength=1.1, labelspacing=.25,
          borderpad=0, handletextpad=.5)

# ---- (b) neither classifier-style caller responds to dosage
ax = axb
for ds, g in KS_.groupby("dataset"):
    g = g.sort_values("spike_burden")
    ax.plot(g.spike_burden, g.tumour_pct, "-o", ms=3.2, lw=1.3, color=COL[ds])
for ds, g in CB.groupby("dataset"):
    g = g.sort_values("spike_burden")
    ax.plot(g.spike_burden, g.aneu_pct, "--s", ms=3.0, lw=1.0, color=COL[ds], alpha=.85)
ax.set_xlabel("arm changes added per cell")
ax.set_ylabel("cells called malignant (%)")
ax.set_xlim(-1, 19.8); ax.set_ylim(-6, 108); ax.set_xticks([0, 5, 10, 15, 18])
ax.set_yticks([0, 25, 50, 75, 100])
ax.plot([], [], "-o", ms=3.2, lw=1.3, color=META_GREY, label="ikarus")
ax.plot([], [], "--s", ms=3.0, lw=1.0, color=META_GREY, label="CopyKAT")
ax.legend(frameon=False, loc="center right", bbox_to_anchor=(1.0, .52), handlelength=1.6,
          labelspacing=.25, borderpad=0, handletextpad=.5)
ax.text(19.3, 30, "colours: same three tissues as (a)", fontsize=6, color=META_GREY,
        ha="right", va="center")

# ---- (c) at matched false-positive rate the raw score beats every rule
ax = axc
lab = {"R1b_topcluster": "highest-scoring cluster", "R4_tirosh": "2-criterion rule (Tirosh)",
       "R1_twogroup": "forced 2-group split", "ikarus": "ikarus"}
m = OP.groupby("method")[["fpr", "sens_hi"]].mean()
pairs = [(r, f"score@FPR({r})") for r in ["R1b_topcluster", "R4_tirosh", "R1_twogroup", "ikarus"]]
XMAX = 86
for i, (r, sc) in enumerate(pairs):
    ax.plot([m.loc[r, "sens_hi"], m.loc[sc, "sens_hi"]], [i, i], color=META_GREY, lw=.9, zorder=1)
    ax.scatter(m.loc[r, "sens_hi"], i, s=34, color=C_RULE, zorder=3)
    ax.scatter(m.loc[sc, "sens_hi"], i, s=34, color=C_SCORE, zorder=3)
    ax.text(XMAX, i, f"{m.loc[r, 'fpr']:.0f}%", fontsize=6.5, color=META_GREY,
            ha="right", va="center")
ax.text(XMAX, len(pairs) - .45, "false calls\nmade to get there", fontsize=6, color=META_GREY,
        ha="right", va="bottom", linespacing=1.25)
ax.set_yticks(np.arange(len(pairs))); ax.set_yticklabels([lab[r] for r, _ in pairs])
ax.set_xlabel("altered cells found (%), among cells with 12–18 arm changes")
ax.scatter([], [], s=34, color=C_RULE, label="published cutoff rule")
ax.scatter([], [], s=34, color=C_SCORE, label="ranking by score, same false-call rate")
ax.legend(frameon=False, loc="lower right", bbox_to_anchor=(.86, -.04),
          scatterpoints=1, handletextpad=.4, labelspacing=.35, borderpad=0)
ax.set_xlim(0, XMAX + 1); ax.set_ylim(-0.6, len(pairs) - 0.25)
ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)

for ax_, ch in [(axa, "a"), (axb, "b"), (axc, "c")]:
    panel_letter(ax_, ch)

# letters: just outside each panel's own labels, (a) and (c) on one vertical line
letters = {ch: (ax_, [t for t in ax_.texts
                     if t.get_text() == ch and t.get_fontweight() in ("bold", 700)][0])
           for ax_, ch in [(axa, "a"), (axb, "b"), (axc, "c")]}
for _, t in letters.values():
    t.set_visible(False)
fig.canvas.draw()
rend, inv = fig.canvas.get_renderer(), fig.transFigure.inverted()
tb = {ch: ax_.get_tightbbox(rend).transformed(inv) for ch, (ax_, _) in letters.items()}
GAP = 0.014
left_col = min(tb["a"].x0, tb["c"].x0) - GAP
for ch, (ax_, t) in letters.items():
    x = left_col if ch in ("a", "c") else tb[ch].x0 - GAP
    p = ax_.get_position()
    t.set_x((x - p.x0) / p.width)
    t.set_visible(True)

fig.savefig("figures/Fig2_spikein.png", dpi=330, bbox_inches="tight")
fig.savefig("figures/Fig2_spikein.pdf", bbox_inches="tight")
print("saved")
