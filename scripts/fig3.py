"""Figure 3 — where the false signal comes from: gene density, not gene families."""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy import stats

from figstyle import apply_figure_style, panel_letter, META_GREY
apply_figure_style(sizes=(8, 7, 6))
M = pd.read_csv("results/T6_density_vs_bias.csv")
BP = pd.read_csv("results/T7_binned_profiles.csv")
XC = pd.read_csv("results/T8_cross_tissue_corr.csv", index_col=0)
FAM = pd.read_csv("results/T9_family_enrichment.csv")

C_MAIN, C_ACC, C_NEG = "#3C6E9F", "#8A4F9E", "#C1662F"

fig = plt.figure(figsize=(7.2, 4.9))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], width_ratios=[1.45, 1],
                      hspace=0.52, wspace=0.4)

# ---- (a) genome-wide mean difference profile, recurrent loci labelled
ax = fig.add_subplot(gs[0, :])
prof = BP.groupby(["chrom", "bin"], as_index=False).agg(diff=("diff", "mean"),
                                                        start=("start", "mean"))
prof["chrom"] = prof["chrom"].astype(str)
CHR = [str(i) for i in range(1, 23)]
prof = prof[prof.chrom.isin(CHR)].copy()
prof["ci"] = prof.chrom.map({c: i for i, c in enumerate(CHR)})
prof = prof.sort_values(["ci", "start"]).reset_index(drop=True)
prof["x"] = np.arange(len(prof))
for i, c in enumerate(CHR):
    sub = prof[prof.chrom == c]
    if i % 2 == 0:
        ax.axvspan(sub.x.min() - .5, sub.x.max() + .5, color="#000000", alpha=.035, lw=0)
ax.bar(prof.x, prof["diff"], width=1.0,
       color=np.where(prof["diff"] < 0, C_MAIN, C_NEG))
ax.axhline(0, color="black", lw=.6)
ax.set_xticks(prof.groupby("chrom", sort=False).x.mean().values)
tl = ax.set_xticklabels(CHR, fontsize=5.5)
# chr19-22 are short enough that their centred labels touch; drop every other
# one of the crowded tail onto a second row so each stays separately readable
for c, t in zip(CHR, tl):
    if c in ("19", "21"):
        t.set_y(t.get_position()[1] - 0.030)
ax.set_xlabel("chromosome")
ax.set_ylabel("apparent copy-number shift\n(called cells − uncalled cells)")
ax.set_xlim(-1, len(prof))
ylo, yhi = ax.get_ylim(); ax.set_ylim(ylo*1.18, yhi*1.05)
top = prof.reindex(prof["diff"].abs().sort_values(ascending=False).index).head(3)
for _, r in top.iterrows():
    ax.annotate("chr%s:%d–%d Mb" % (r.chrom, r["bin"]*10, r["bin"]*10+10),
                (r.x, r["diff"]), fontsize=5.5, ha="center",
                va="top" if r["diff"] < 0 else "bottom",
                xytext=(0, -3 if r["diff"] < 0 else 3), textcoords="offset points",
                color="#333333")
ax.text(0.004, 0.05, "apparent loss", transform=ax.transAxes, fontsize=6, color=C_MAIN)
ax.text(0.004, 0.90, "apparent gain", transform=ax.transAxes, fontsize=6, color=C_NEG)
panel_letter(ax, "a")

# ---- (b) gene density explains the direction of the bias
ax = fig.add_subplot(gs[1, 0])
ax.scatter(M.n_genes, M["mean_diff"], s=8, color=C_MAIN, alpha=.45, edgecolors="none")
sl, ic, r, p, _ = stats.linregress(M.n_genes, M["mean_diff"])
xx = np.linspace(M.n_genes.min(), M.n_genes.max(), 50)
ax.plot(xx, ic + sl * xx, color=C_ACC, lw=1.3)
ax.axhline(0, color=META_GREY, lw=.6)
ax.set_xlabel("genes per 10-Mb bin")
ax.set_ylabel("apparent copy-number shift")
ax.annotate(f"Spearman $\\rho$ = {stats.spearmanr(M.n_genes, M['mean_diff'])[0]:.2f}\n"
            f"$P$ = {stats.spearmanr(M.n_genes, M['mean_diff'])[1]:.1e}",
            (0.97, 0.95), xycoords="axes fraction", ha="right", va="top", fontsize=6.5)
panel_letter(ax, "b", dx=-0.30)

# ---- (c) the same bins misbehave in independent tissues
ax = fig.add_subplot(gs[1, 1])
PLAIN = {"CD_epithelial": "colon (Crohn's, epi.)", "CD_immune": "colon (Crohn's, imm.)",
         "PSO_skin": "skin (psoriasis)", "UC_colon": "colon (colitis)"}
lab = [PLAIN.get(c, c.replace("TS_", "").replace("_", " ").lower()) for c in XC.columns]
im = ax.imshow(XC.values, cmap="RdBu_r", vmin=-.6, vmax=.6)
ax.set_xticks(range(len(lab))); ax.set_xticklabels(lab, rotation=60, ha="right", fontsize=5.5)
ax.set_yticks(range(len(lab))); ax.set_yticklabels(lab, fontsize=5.5)
ax.text(0, 1.005, "typical tissue-to-tissue agreement: $\\rho$ = %.2f" % np.median(XC.values[np.triu_indices_from(XC.values, k=1)]), transform=ax.transAxes, fontsize=6, color=META_GREY, va="bottom")
cb = fig.colorbar(im, ax=ax, fraction=.046, pad=.04)
cb.set_label("agreement between two\ntissues (Spearman $\\rho$)", fontsize=6)
cb.ax.tick_params(labelsize=5.5)
panel_letter(ax, "c")

fig.savefig("figures/Fig3_mechanism.png", dpi=330, bbox_inches="tight")
fig.savefig("figures/Fig3_mechanism.pdf", bbox_inches="tight")
print("saved")
