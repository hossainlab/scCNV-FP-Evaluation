"""Derive manuscript/main_pandoc.tex from manuscript/main.tex.

pandoc carries \\cite keys through as bare keys and drops authblk's \\affil, so the
DOCX/HTML conversions lose the citations and the affiliation block. This rewrite
resolves every \\cite to the number of its \\bibitem, turns thebibliography into an
ordered list (so the numbering in the text and in the reference list are the same
sequence), restates the affiliation and correspondence lines as plain paragraphs,
and points \\includegraphics at the raster copies pandoc can embed.
"""
import re

# WeasyPrint has no MathML support, so pandoc's MathML output reaches the PDF as
# the MathML fallback text followed by the raw TeX annotation ("$\rho = 0.54$"
# printed twice over). The manuscript's maths is 27 distinct inline spans of
# Greek letters, subscripts and comparisons, all of which have exact Unicode
# spellings, so each span is rewritten as text before conversion. Coverage is
# asserted: an unlisted span stops the build rather than leaking TeX.
MATH = {
    r"\alpha = 5\%": "\u03b1 = 5\\%", r"\alpha = 1\%": "\u03b1 = 1\\%",
    r"\alpha = 10\%": "\u03b1 = 10\\%", r"\alpha \le 5\%": "\u03b1 \u2264 5\\%",
    r"1-\alpha": "1 \u2212 \u03b1", r"\alpha": "\u03b1",
    r"\rho = -0.53": "\u03c1 = \u22120.53", r"\rho = 0.54": "\u03c1 = 0.54",
    r"\rho": "\u03c1", r"P = 2\times10^{-21}": "\\emph{P} = 2\u00d710\u207b\u00b2\u00b9",
    r"\rightarrow": "\u2192", r"'": "\u2032", r"\pm": "\u00b1",
    r"10\%": "10\\%", r"1.5": "1.5", r"0.5": "0.5",
    r"t": "\\emph{t}", r"b": "\\emph{b}", r"s": "\\emph{s}",
    r"d": "\\emph{d}", r"r": "\\emph{r}",
    r"b \in \{1,2,3,5,8,12,18\}": "\\emph{b} \u2208 {1, 2, 3, 5, 8, 12, 18}",
    r"s_i": "\\emph{s}\\textsubscript{i}", r"d_i": "\\emph{d}\\textsubscript{i}",
    r"F_t": "\\emph{F}\\textsubscript{t}",
    r"r_i = s_i - \hat{s}(d_i)":
        "\\emph{r}\\textsubscript{i} = \\emph{s}\\textsubscript{i} \u2212 "
        "\u015d(\\emph{d}\\textsubscript{i})",
    r"p_i = 1 - F_t(r_i)":
        "\\emph{p}\\textsubscript{i} = 1 \u2212 \\emph{F}\\textsubscript{t}"
        "(\\emph{r}\\textsubscript{i})",
    r"p_i \le \alpha": "\\emph{p}\\textsubscript{i} \u2264 \u03b1",
}

FIGS = ["Fig1_false_positives", "Fig2_spikein", "Fig3_mechanism", "Fig4_calibration"]

t = open("manuscript/main.tex", encoding="utf-8").read()
keys = re.findall(r"\\bibitem\{([^}]+)\}", t)
num = {k: i + 1 for i, k in enumerate(keys)}


def cite_repl(m):
    ks = [k.strip() for k in m.group(1).split(",")]
    missing = [k for k in ks if k not in num]
    assert not missing, f"citation with no bibitem: {missing}"
    return "[" + ",".join(str(num[k]) for k in ks) + "]"


p, n_cite = re.subn(r"\\cite[pt]?\{([^}]+)\}", cite_repl, t)
p = p.replace("\\maketitle\n",
              "\\maketitle\n\n\\noindent Affiliation: [Institution, Department, City, Country]\\\\\n"
              "\\noindent Correspondence: \\texttt{[corresponding author e-mail]}\n")
p = p.replace("\\begin{thebibliography}{99}\n\\small", "\\begin{enumerate}\n\\itemsep2pt")
p = re.sub(r"\\bibitem\{[^}]+\}\s*", "\\\\item ", p)
p = p.replace("\\end{thebibliography}", "\\end{enumerate}")
p = p.replace("\\linenumbers", "")


def math_repl(m):
    body = m.group(1).strip()
    assert body in MATH, f"inline math with no Unicode spelling: {body!r}"
    return MATH[body]


p, n_math = re.subn(r"\$([^$]+)\$", math_repl, p)
assert "$" not in p, "unpaired $ left in conversion source"
for v in FIGS:
    p = p.replace("{%s.pdf}" % v, "{%s.png}" % v)
open("manuscript/main_pandoc.tex", "w", encoding="utf-8").write(p)
print(f"citations resolved: {n_cite} | references: {len(keys)} | "
      f"bibitem remaining: {p.count('bibitem')} | math spans as text: {n_math}")
