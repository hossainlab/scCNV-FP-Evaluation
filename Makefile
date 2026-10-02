# Run every target from the repository root.
PY ?= python
RS ?= Rscript

.PHONY: figures genes panel infercnv ikarus copykat spikein manuscript-tex

figures:            ## Fig 1-4 from results/ (no raw data needed)
	$(PY) scripts/fig1.py
	$(PY) scripts/fig2.py
	$(PY) scripts/fig3.py
	$(PY) scripts/fig4.py

genes:              ## Ensembl coordinates -> gene_positions.csv, loci_*.csv
	$(PY) scripts/gene_table.py

panel:              ## download + subsample CELLxGENE datasets -> panel/*.h5ad
	$(PY) scripts/assemble.py

infercnv: genes panel    ## inferCNV scores + decision rules -> calls/
	$(PY) scripts/run_infercnv.py

ikarus: panel       ## needs pyshim/ and ikarus_assets/ (see README)
	$(PY) scripts/run_ikarus.py

copykat: panel      ## export to Matrix Market, then CopyKAT in R
	$(PY) scripts/export_for_r.py
	$(RS) R/run_copykat_all.R

spikein: genes panel     ## needs chrom_arms.csv (see README)
	$(PY) scripts/build_spikein.py

manuscript-tex:     ## pandoc-friendly copy of manuscript/main.tex
	$(PY) scripts/make_pandoc_tex.py
