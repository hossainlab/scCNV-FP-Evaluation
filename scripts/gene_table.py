"""Fetch Ensembl gene coordinates from BioMart and build gene_positions.csv plus the
clustered gene-family annotation (loci.build). Same query as step 1 of regen_lung.py.

Note: www.ensembl.org serves the *current* Ensembl release, which was not recorded
for the published run; coordinates may shift slightly between releases.
"""
import os, sys, urllib.request, urllib.parse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import loci

if not os.path.exists("gene_pos_biomart.tsv"):
    xml = ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
           '<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1"'
           ' count="" datasetConfigVersion="0.6">'
           '<Dataset name="hsapiens_gene_ensembl" interface="default">'
           '<Attribute name="ensembl_gene_id"/><Attribute name="external_gene_name"/>'
           '<Attribute name="chromosome_name"/><Attribute name="start_position"/>'
           '<Attribute name="end_position"/><Attribute name="gene_biotype"/>'
           '</Dataset></Query>')
    u = "https://www.ensembl.org/biomart/martservice?" + urllib.parse.urlencode({"query": xml})
    with urllib.request.urlopen(u, timeout=900) as r:
        txt = r.read().decode()
    open("gene_pos_biomart.tsv", "w", newline="").write(txt)
    print("biomart rows", txt.count("\n") - 1, flush=True)
g, fam, gf = loci.build()
print("genes", len(g), "families", len(fam), flush=True)
