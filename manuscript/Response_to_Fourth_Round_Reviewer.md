# Response to Fourth-Round Reviewer Report

Manuscript: *NT5DC2 is an independently prognostic host gene of differentially expressed
circRNAs in hepatocellular carcinoma: a cross-cohort circRNA–miRNA–mRNA network analysis*
Journal: Scientific Reports · Recommendation received: Accept after technical check

All four residual items are addressed. RI3 required re-running the univariable analyses and
changed several reported values, so we describe it in full.

---

## RI1 — Abstract one word over the limit

**Corrected: the abstract is now 194 words** (was 201). We trimmed three phrases without
dropping any reported quantity:

- "evaluated the prognostic value of their host-gene mRNAs" → "evaluated their host-gene
  mRNAs";
- "Because neither cohort contains circRNA data …, the markers identified are …" →
  "Neither cohort contains circRNA data …, so the markers identified are …";
- "Comparison with ENCORI/starBase showed minimal overlap … indicating that …" →
  "ENCORI/starBase comparison showed minimal overlap … so …".

The Additional information block now reads **194 words (limit 200)**. We also keep
`count_abstract_words.py` in the archive so that the quoted figure can be regenerated from
the manuscript file rather than maintained by hand — the disagreement between 195, 199, 201
and 216 in earlier rounds came precisely from a hand-maintained number.

## RI2 — Reference 5: wrong issue and article number

**Corrected.** We verified the publisher record independently before changing anything: the
article is **Cell Death Dis 2018;9(2):175** (PMID 29415990, published 7 February 2018), not
9(3):306. The DOI 10.1038/s41419-017-0204-3 is unchanged and resolves correctly. The
reference list now reads "*Cell Death Dis.* 2018;9:175".

## RI3 — Univariable sample size (367 exceeded the funnel value of 365)

**The reviewer was right, and this required re-running the univariable analyses.**

The number 367 came from `step2_tcga_survival.py`, whose sample-type test was
`seg[0] == '0'` on the TCGA barcode. That test treats every `0X` code as tumour, so three
**recurrent** tumour samples (sample type 02) were counted alongside the primary tumours.
The univariable Cox and log-rank analyses had therefore never been re-run on the
tumour-only export we prepared for the multivariable models in the previous round — the
file timestamp on the resulting table predated every correction. This is the same
barcode rule that produced the "374 tumours" figure in an earlier round, and we have
corrected it in both places.

**Correction applied.** The barcode test now requires the `01` sample-type code specifically.
After re-running the univariable analyses on primary tumours only:

| | Previously reported | Corrected (n = 364, 130 events) |
|---|---|---|
| Samples analysed | 367 | **364** |
| NT5DC2 | HR 1.509, P = 1.62 × 10⁻⁵, log-rank 0.013 | **HR 1.536, P = 7.49 × 10⁻⁶, log-rank 0.0076** |
| HSD17B6 | 0.756, 1.23 × 10⁻³, 0.067 | **0.762, 1.58 × 10⁻³, 0.075** |
| BCLAF1 | 1.133, 0.086, 0.248 | **1.141, 0.074, 0.225** |
| SREBF2 | 1.077, 0.201, 0.427 | **1.100, 0.143, 0.289** |
| LRCH1 | 0.893, 0.092, 0.344 | **0.904, 0.120, 0.320** |
| CD5 | 0.877, 0.063, 0.0019 | **0.865, 0.045, 0.0016** |
| Genes with Cox P < 0.05 in TCGA | 37 | **40** |
| Genes significant by both Cox and log-rank | 25 | **26** |

**Why 364 rather than the 365 in the funnel.** The funnel counts primary tumours that have
survival data (365). The univariable models additionally rebuild survival time per sample from
`days_to_death` or `days_to_last_follow_up` and exclude cases whose vital status is not
reported or whose corresponding time field is missing, rather than reading the pre-computed
overall-survival field; one primary tumour is excluded on that basis. We now state this
explicitly in Section 4.1 rather than leaving the one-sample discrepancy unexplained.

**One substantive consequence, which we flag rather than bury.** CD5's univariable Cox P value
in TCGA-LIHC falls from 0.063 to **0.045**, so CD5 now crosses the 0.05 threshold and is
significant by Cox regression in both cohorts before adjustment. We have revised the
consequently affected statements in Section 2.2, Section 2.3, Section 3.1, Table 1 and the
Fig. 5 legend to say so. CD5 remains non-significant after multivariable adjustment
(adjusted HR 0.88, P = 0.228), so the graded conclusion of the paper is unchanged: NT5DC2
is the only gene independently prognostic in both cohorts.

## RI4 — Table order at the end of the manuscript

**Corrected.** The cause was in our own conversion script rather than the manuscript text:
the table-collection step recognised table captions only when they were written as markdown
headings (`## Table 1 | …`), whereas we had reformatted all three captions as bold
paragraphs (`**Table 1 | …**`). The captions were therefore not captured, and the tables were
emitted in order of first appearance in the text — Table 2 (Section 2.4), Table 3 (Section
4.1), then Table 1 (Section 2.2) — giving the observed 2, 3, 1.

We have made two changes:

1. the converter now recognises bold-paragraph captions and captures the table number, and
2. the collected tables are sorted by number before being written, so that main Tables 1–3
   precede Supplementary Tables S1–S6 regardless of where they appear in the text.

We verified the generated file: the tables block now reads Table 1, Table 2, Table 3, S1, S2,
S3, S4, S5, S6. All three main tables, including the Methods data-sources table, are now in
that block.

---

## Archive

The Zenodo archive has been updated to **v1.3.0** with the re-run univariable Cox and
log-rank results (`results/tcga_candidate_gene_stats.tsv`), the corrected
`step2_tcga_survival.py`, the corrected converter, and regenerated manifest checksums. The
manuscript and all supporting files continue to cite the concept DOI
(10.5281/zenodo.23293107), which resolves to the most recent version.

We are grateful for RI3 in particular. It was the third instance of the same underlying
cause — a cohort-selection rule applied inconsistently between analysis steps — and it was
surfaced only because the reviewer compared two numbers arithmetically rather than reading
either in isolation. We have now audited all three analysis entry points (`step2`,
`step3`, `step14`) and they use the same tumour-only definition.
