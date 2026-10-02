# bioRxiv submission sheet — null-control audit of malignant-cell calling

Prepared to accompany `OPERON_null_control_audit_preprint.pdf`. Everything below is
either ready to paste into the submission form or flagged as needing your input.

## Files to upload

| Form field | File | Note |
|---|---|---|
| Manuscript | `OPERON_null_control_audit_preprint.pdf` | Single PDF, text + 4 inline figures + captions, 14 pages, 1.2 MB. Within the 40 MB limit that the preprint servers impose. |
| Supplemental material | `OPERON_supplement_code_and_tables.zip` | Posted as-is; 46 entries (result tables, analysis and figure scripts, regenerated per-cell scores, figure files). |
| Supplemental material (optional) | `OPERON_null_control_audit_preprint.tex` | LaTeX source. bioRxiv requires TeX to be converted to PDF before submission, but the source may be attached as supplemental material. |

The DOCX route is an alternative, not an addition: upload either the PDF, or the
Word file plus separate figure files for server-side conversion. Do not upload both.

## Metadata fields

**Title**
Malignant-cell callers assign tumour labels to tumour-free tissue: a null-control
audit and a calibrated decision procedure

**Abstract for the metadata box** (270 words — the in-PDF abstract is 391 words, which is
fine in the manuscript but long for the form and over the limit of most journals
you would submit to next)

Single-cell cancer studies partition cells into malignant and non-malignant compartments using copy-number inference or supervised classifiers, and every downstream claim inherits that partition's errors. They are rarely evaluated on tissue containing no tumour, so their false-positive rate is unknown. We assembled 14 tumour-free scRNA-seq datasets (36,528 cells, 76 donors, 10 healthy organs and 4 inflammatory-disease tissues) and ran three widely used malignant-cell callers under five published decision rules; every positive call in this panel is false by construction. False malignant calls occurred in every tissue: 0.5-19.3% for the two-criterion rule of Tirosh and colleagues and 4.5-78.0% for a forced two-group split. CopyKAT called 3.0-63.0% of cells aneuploid and ikarus 0-100%; pooled over the panel, each assigned roughly 31% of tumour-free cells to the malignant class. In-silico positive controls imposing arm-level dosage changes show per-cell discrimination near chance for single events (AUC 0.53), reaching 0.85 only at 18 arm events; at matched false-positive rate, ranking cells by the continuous score detects about twice as many aneuploid cells as the published rules, which therefore choose a poor operating point rather than adding information. The false calls share a reproducible pseudo-karyotype, correlated across independent tissues (median Spearman rho = 0.54) and predicted by local gene density (rho = -0.53, P = 2x10-21). A cell-type-conditional empirical null fitted on tumour-free reference cells holds the realized false-positive rate at or below its nominal target in 13 of 14 datasets (median 0.03% at alpha = 5%); the single failure, a different-protocol dataset, is repaired by drawing the null from held-out donors (27.1% to 0.03%). We provide the calibration as a caller-agnostic layer consuming any per-cell score.

**Keywords**
single-cell RNA-seq; copy-number inference; malignant-cell classification; negative controls; false-positive rate; calibration

**Suggested subject category**
Bioinformatics (primary). Cancer Biology or Genomics are defensible alternatives;
the paper is a methods audit, so the bioinformatics category matches the audience
most likely to act on it.

**License**
CC BY 4.0 if you want the calibration layer reused and redistributed freely, which
is the point of releasing it. CC BY-NC-ND is the most common choice on the server
but forbids derivative works, which would sit awkwardly with a paper whose
contribution is a reusable software layer.

## Declarations — drafted, confirm before posting

These are now written into the manuscript. They are statements only you can make,
so read them and correct anything that is not true of your situation:

- **Competing interests**: "The authors declare no competing interests."
- **Funding**: "This work received no specific funding. All input data are public
  and all analyses were run on commodity hardware."
- **Author contributions**: all authors contributed to conception, analysis,
  interpretation, drafting and revision, and approved the submitted version.
- **AI assistance** (in Methods, Software and availability): data processing, figure
  generation and drafting were done with the assistance of an AI coding assistant
  under author direction and review; every reported number is reproducible from the
  accompanying code, and the authors are responsible for the content.

## Still needed from you

1. **Author names** — currently the bracketed placeholder in the title block.
2. **Affiliation** — institution, department, city, country.
3. **Corresponding author e-mail** — the form also requires contact details for the
   corresponding author separately from the manuscript.
4. **ORCID iD per author** — optional but worth linking, since it populates each
   author's ORCID record automatically.

Give me 1-3 and I will rebuild all three formats in one step.

## Before you click submit

- Posting is one-way: a preprint cannot be removed from bioRxiv once posted.
  Revisions are supported and keep the original DOI, but the first version stays
  visible.
- Every submission is screened in-house for completeness, scope and article type,
  and automatically for plagiarism and prior online appearance. Nothing in this
  manuscript has been posted elsewhere.
- The audit's weakest evidence is its positive controls: real-tumour benchmarks were
  not run, so sensitivity rests on in-silico spike-ins. This is stated in
  Limitations. It is the objection to expect, and it is a scientific decision rather
  than a formatting one.
