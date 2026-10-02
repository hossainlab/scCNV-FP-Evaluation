# Null-control audit of malignant-cell calling

Code, result tables and figures for the preprint
*Malignant-cell callers assign tumour labels to tumour-free tissue: a null-control
audit and a calibrated decision procedure.*

Three malignant-cell callers (inferCNV via infercnvpy, CopyKAT, ikarus) are run on
14 tumour-free scRNA-seq datasets from CELLxGENE (36,528 cells, 76 donors), so every
malignant call is a false positive by construction. In-silico arm-level spike-ins
supply positive controls, and a cell-type-conditional empirical null
(`scripts/calibrate.py`) restores nominal false-positive control.

## Repository layout

```
config/          panel_def.json — the 14 panel datasets (CELLxGENE URLs, assay filter)
metadata/        dataset discovery and citation metadata (CELLxGENE manifests, refs)
scripts/         Python pipeline, analysis modules and figure scripts
R/               CopyKAT runner and R package installers
results/         result tables T0–T10, headline_numbers.json, panel_summary.json
calls/           per-cell caller output (only TS_Lung inferCNV is shipped)
figures/         published Figures 1–4 (PNG 330 dpi + vector PDF)
manuscript/      LaTeX source, DOCX, print CSS
docs/            bioRxiv submission sheet and the original supplement README
```

All scripts use paths relative to the repository root — **run everything from the
root**, e.g. `python scripts/fig1.py`.

## Quick start: regenerate the figures (no raw data needed)

```bash
python -m venv .venv && . .venv/bin/activate    # Windows: .venv\Scripts\activate
pip install -r requirements-figures.txt
make figures       # or: python scripts/fig1.py ... python scripts/fig4.py
```

Figures are rebuilt from `results/*.csv` and `calls/TS_Lung_infercnv.parquet` and
overwrite `figures/`. The plotted data are identical to the published figures;
spacing may differ slightly because the shared style helpers (`scripts/figstyle.py`)
were reconstructed — they were not part of the original supplement.

## Full pipeline

Environment: `conda env create -f environment.yml` (Python 3.11, infercnvpy 0.6.1,
ikarus 0.0.3, R 4.4.3). Install CopyKAT 1.2.5 in R with `Rscript R/install_callers.R`
(or `R/install_callers_binary.R` on Windows without a compiler; set `R_BUILDTOOLS` to
a conda env containing `m2w64-toolchain` if you need to compile).

| Step | Command | Output |
|---|---|---|
| 1. Gene coordinates | `python scripts/gene_table.py` | `gene_positions.csv`, `loci_*.csv` |
| 2. Panel assembly (~15 GB downloaded one file at a time; peak disk ~3.2 GB) | `python scripts/assemble.py` | `panel/*.h5ad`, `results/panel_summary.json` |
| 3. inferCNV scoring + decision rules | `python scripts/run_infercnv.py` | `calls/*_infercnv.{parquet,json}`, `*_profile.npz` |
| 4. ikarus | `python scripts/run_ikarus.py` | `calls/*_ikarus.{csv,json}` |
| 5. CopyKAT | `python scripts/export_for_r.py` then `Rscript R/run_copykat_all.R` | `rmat/`, `calls/*_copykat.{csv,json}` |
| 6. Spike-in positive controls | `python scripts/build_spikein.py`, then steps 3–5 with `spike/*.h5ad` as arguments | `spike/` |
| 7. Aggregate into tables | *not included — see Known gaps* | `results/T*.csv` |
| 8. Figures | `make figures` | `figures/` |
| 9. Manuscript | `python scripts/make_pandoc_tex.py` then pandoc | `manuscript/main_pandoc.tex` |

`scripts/regen_lung.py` is a self-contained check: it rebuilds the TS_Lung panel and
inferCNV scores from scratch and compares the R3_refSD false-positive rate with the
published value in `results/T1_infercnv_fpr.csv`.

Library modules (imported, not run directly): `prep.py`, `dl.py`, `loci.py`,
`spikein.py`, `calibrate.py` (empirical-null calibration with GPD tail),
`attribution.py` (cell-type and genomic-locus attribution), `archetype.py`
(label-free conditioning).

## Known gaps

The original supplement does not contain everything needed to rerun end to end:

- **Table aggregation (step 7).** No script writes `results/T*.csv` or
  `headline_numbers.json` from `calls/`; the modules `calibrate.py`, `attribution.py`
  and `archetype.py` provide the functions but not the driver.
- **`chrom_arms.csv`** (chromosome-arm boundaries, columns
  `chrom,p_start,p_end,q_start,q_end`) required by `spikein.py` is not shipped.
- **ikarus assets.** `run_ikarus.py` expects `pyshim/` and
  `ikarus_assets/{signatures.gmt,core_model.joblib}` in the root.
- **Ensembl release** for the gene coordinates was not recorded.
- **Dependency versions** other than those in `environment.yml` were not recorded.

## Data

All input data are public CELLxGENE Discover datasets; URLs are in
`config/panel_def.json` and provenance in `metadata/`.

## License and citation

To be added.
