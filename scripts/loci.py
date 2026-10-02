"""Build the gene-position table and the clustered co-regulated gene-family annotation
used to attribute spurious CNV segments to genomic loci."""
import re
import numpy as np
import pandas as pd

CANON = [str(i) for i in range(1, 23)] + ["X", "Y"]

def gene_table(path="gene_pos_biomart.tsv"):
    g = pd.read_csv(path, sep="\t")
    g.columns = ["ensembl", "symbol", "chrom", "start", "end", "biotype"]
    g = g[g["chrom"].astype(str).isin(CANON)].copy()
    g["chrom"] = g["chrom"].astype(str)
    g["start"] = g["start"].astype(int); g["end"] = g["end"].astype(int)
    g["mid"] = (g["start"] + g["end"]) // 2
    g["symbol"] = g["symbol"].fillna("")
    g = g.drop_duplicates("ensembl").reset_index(drop=True)
    return g

# --- a-priori families: named loci known to be co-regulated and physically clustered,
# --- plus dispersed programs (ribosomal, cell cycle, interferon) kept as contrast sets.
CURATED = {
    # immunoglobulin loci
    "IG_heavy":        (r"^IGH[VDJGAEM]", None),
    "IG_kappa":        (r"^IGK[VDJC]", None),
    "IG_lambda":       (r"^IGL[VDJC]", None),
    # T-cell receptor loci
    "TCR_alpha":       (r"^TR[AD][VDJC]", None),
    "TCR_beta":        (r"^TRB[VDJC]", None),
    "TCR_gamma":       (r"^TRG[VJC]", None),
    # MHC
    "HLA_classI_II":   (r"^HLA-", "6"),
    # histone cluster
    "histone_HIST1":   (r"^(H1-|H2A|H2B|H3C|H4C|H2AC|H2BC|HIST1)", "6"),
    # epidermal differentiation complex (1q21.3)
    "EDC_S100":        (r"^S100A", "1"),
    "EDC_SPRR_LCE":    (r"^(SPRR|LCE|IVL|FLG|LOR|CRCT1|TCHH)", "1"),
    # keratins
    "KRT_typeII_12q":  (r"^KRT", "12"),
    "KRT_typeI_17q":   (r"^KRT", "17"),
    # protocadherins
    "PCDH_5q31":       (r"^PCDH[ABG]", "5"),
    # haemoglobin loci
    "HBB_cluster":     (r"^HB[BDGEQ]", "11"),
    "HBA_cluster":     (r"^HB[AMZQ]", "16"),
    # secretory / digestive programs
    "pancreatic_enzymes": (r"^(PRSS1|PRSS2|PRSS3|CTRB1|CTRB2|CTRC|CTRL|CELA2A|CELA2B|CELA3A|CELA3B|CPA1|CPA2|CPB1|CLPS|CLPSL1|PNLIP|PNLIPRP1|PNLIPRP2|CEL|SYCN|AMY2A|AMY2B|AMY1A|PLA2G1B)$", None),
    "REG_cluster":     (r"^REG[1-4ABG]", "2"),
    "MUC_11p15":       (r"^MUC(2|5AC|5B|6)$", "11"),
    "DEFA_8p23":       (r"^DEFA", "8"),
    "DEFB_20p13":      (r"^DEFB", "20"),
    "surfactant":      (r"^SFTP[ABCD]", None),
    "APOA_11q23":      (r"^APO(A1|A4|A5|C3)$", "11"),
    "ALB_AFP_4q13":    (r"^(ALB|AFP|AFM|ANG|RASSF6)$", "4"),
    "CYP_19q13":       (r"^CYP2[ABFGS]", "19"),
    "UGT1A_2q37":      (r"^UGT1A", "2"),
    "MT_16q13":        (r"^MT1[A-X]$|^MT2A$", "16"),
    "CEACAM_19q13":    (r"^(CEACAM|PSG)", "19"),
    "chemokine_CCL_17q": (r"^CCL", "17"),
    "chemokine_CXCL_4q": (r"^CXCL", "4"),
    "granzyme":        (r"^GZM", None),
    "complement_C1":   (r"^C1Q[ABC]$", "1"),
    "HOXA_7p15":       (r"^HOXA", "7"),
    "HOXB_17q21":      (r"^HOXB", "17"),
    "HOXC_12q13":      (r"^HOXC", "12"),
    "HOXD_2q31":       (r"^HOXD", "2"),
    "olfactory_recept":(r"^OR[0-9]+[A-Z]", None),
    "ZNF_19q13":       (r"^ZNF", "19"),
    "LILR_KIR_19q13":  (r"^(LILR|KIR[23]D)", "19"),
    "FCGR_1q23":       (r"^FCGR[23]", "1"),
    "SERPINB_18q21":   (r"^SERPINB", "18"),
    "GAGE_XAGE_Xp11":  (r"^(GAGE|XAGE|MAGEA)", "X"),
    "LGALS_lectins":   (r"^LGALS", None),
    # dispersed programs (contrast: co-regulated but NOT clustered)
    "ribosomal_protein": (r"^(RPL|RPS)[0-9]", None),
    "interferon_stim":   (r"^(ISG15|ISG20|IFIT[1-5]|IFI6|IFI27|IFI44|IFI44L|IFIH1|IFITM[1-3]|MX1|MX2|OAS[1-3]|OASL|STAT1|STAT2|IRF7|IRF9|RSAD2|XAF1|EIF2AK2|BST2|CMPK2|HERC5|HERC6|LY6E|SIGLEC1|USP18|SAMD9|SAMD9L|EPSTI1|PARP9|DDX60|SPATS2L)$", None),
    "cell_cycle_S":      (r"^(MCM[2-7]|PCNA|TYMS|RRM1|RRM2|UNG|GINS2|CDC45|CDC6|CLSPN|DTL|E2F8|EXO1|FEN1|GMNN|HELLS|MSH2|NASP|POLA1|POLD3|PRIM1|RFC2|RPA2|SLBP|UBR7|USP1|WDR76|CHAF1B|BRIP1|CASP8AP2|ATAD2|BLM|CDCA7|CENPU|DSCC1|NBN|PCNA|POLD3|RAD51|RAD51AP1|UHRF1)$", None),
    "cell_cycle_G2M":    (r"^(CCNB1|CCNB2|CDK1|TOP2A|MKI67|BUB1|BUB1B|AURKA|AURKB|PLK1|CENPA|CENPE|CENPF|TPX2|UBE2C|BIRC5|KIF11|KIF20B|KIF23|KIF2C|NUSAP1|ANLN|ASPM|CKAP2|CKAP2L|CKAP5|CDC20|CDC25C|CDCA2|CDCA3|CDCA8|DLGAP5|ECT2|G2E3|GAS2L3|HJURP|HMGB2|HMMR|LBR|MZT1|NDC80|NEK2|PSRC1|RANGAP1|SMC4|TACC3|TTK|TUBB4B)$", None),
    "mitochondrial_enc": (r"^MT-", None),
    "HSP_stress":        (r"^(HSPA1A|HSPA1B|HSPA6|HSPB1|HSPH1|DNAJB1|DNAJA1|HSP90AA1|HSPE1|HSPD1|BAG3|AHSA1|CHORDC1)$", None),
}

def _tight_cluster(sub, window=10_000_000, min_n=4):
    """Largest number of family genes falling inside any `window`-wide interval on one chromosome."""
    best = (0, None, None, None)
    for ch, grp in sub.groupby("chrom"):
        pos = np.sort(grp["mid"].values)
        j = 0
        for i in range(len(pos)):
            while pos[i] - pos[j] > window:
                j += 1
            if i - j + 1 > best[0]:
                best = (i - j + 1, ch, int(pos[j]), int(pos[i]))
    return best if best[0] >= min_n else (best[0], best[1], best[2], best[3])

def build(gene_pos="gene_pos_biomart.tsv", out_genes="loci_gene_annotation.csv",
          out_fam="loci_families.csv"):
    g = gene_table(gene_pos)
    fam_rows, gene_fam = [], []
    for name, (rx, ch) in CURATED.items():
        m = g["symbol"].str.match(rx, na=False)
        if ch is not None:
            m &= (g["chrom"] == ch)
        sub = g[m]
        if len(sub) == 0:
            continue
        n, cch, a, b = _tight_cluster(sub)
        frac = n / len(sub)
        fam_rows.append({"family": name, "kind": "curated", "n_genes": len(sub),
                         "cluster_chrom": cch, "cluster_start": a, "cluster_end": b,
                         "n_in_cluster": n, "frac_clustered": round(frac, 3),
                         "span_mb": round(((b - a) / 1e6) if a is not None else np.nan, 2),
                         "clustered": bool(frac >= 0.6 and n >= 4)})
        for _, r in sub.iterrows():
            gene_fam.append({"ensembl": r["ensembl"], "symbol": r["symbol"], "chrom": r["chrom"],
                             "start": r["start"], "end": r["end"], "mid": r["mid"],
                             "family": name, "kind": "curated"})
    # --- data-driven: symbol root families that are physically clustered
    root = g["symbol"].str.replace(r"[0-9\-].*$", "", regex=True)
    g2 = g.assign(root=root)
    seen = {r["family"] for r in fam_rows}
    for rt, grp in g2.groupby("root"):
        if len(rt) < 2 or len(grp) < 5:
            continue
        n, cch, a, b = _tight_cluster(grp, window=5_000_000, min_n=5)
        frac = n / len(grp)
        if n >= 5 and frac >= 0.6:
            name = f"dd_{rt}_{cch}"
            if name in seen:
                continue
            fam_rows.append({"family": name, "kind": "data_driven", "n_genes": len(grp),
                             "cluster_chrom": cch, "cluster_start": a, "cluster_end": b,
                             "n_in_cluster": n, "frac_clustered": round(frac, 3),
                             "span_mb": round((b - a) / 1e6, 2), "clustered": True})
            for _, r in grp.iterrows():
                gene_fam.append({"ensembl": r["ensembl"], "symbol": r["symbol"], "chrom": r["chrom"],
                                 "start": r["start"], "end": r["end"], "mid": r["mid"],
                                 "family": name, "kind": "data_driven"})
    fam = pd.DataFrame(fam_rows).sort_values(["kind", "n_genes"], ascending=[True, False])
    gf = pd.DataFrame(gene_fam)
    g.to_csv("gene_positions.csv", index=False)
    fam.to_csv(out_fam, index=False)
    gf.to_csv(out_genes, index=False)
    return g, fam, gf
