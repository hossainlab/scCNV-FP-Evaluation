# Null-control audit of malignant-cell calling — submission bundle

## Manuscript
- `OPERON_null_control_audit_preprint.pdf` — typeset preprint (bioRxiv upload file)
- `OPERON_null_control_audit_preprint.docx` — same content, editable
- `OPERON_null_control_audit_preprint.tex` — LaTeX source; figures are referenced as
  FIG1-FIG4 and resolve to `figures/Fig*.pdf`

Author names and affiliation are placeholders and must be filled in before submission.

## Figures
Fig1 false positives; Fig2 spike-in positive controls; Fig3 mechanism; Fig4 calibration.
PNG (330 dpi) and vector PDF versions of each are in `figures/`.

## Supplement
`OPERON_supplement_code_and_tables.zip` contains every analysis script, the result
tables the manuscript numbers are read from (`results/T*.csv`,
`results/headline_numbers.json`), the figure scripts and the LaTeX source.

## Reproducing
Callers: infercnvpy 0.6.1 (Python), ikarus 0.0.3 (Python), CopyKAT 1.2.5 (R 4.4.3).
Entry points: `run_infercnv.py`, `run_ikarus.py`, `export_for_r.py` + `run_copykat_all.R`.
Calibration layer: `calibrate.py` (type-conditional empirical null with generalized-Pareto
tail). Figures: `fig1.py` ... `fig4.py`.
