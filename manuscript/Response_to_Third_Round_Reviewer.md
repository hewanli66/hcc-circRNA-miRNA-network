# Response to Third-Round Reviewer Report

Manuscript: *NT5DC2 is an independently prognostic host gene of differentially expressed
circRNAs in hepatocellular carcinoma: a cross-cohort circRNA–miRNA–mRNA network analysis*
Journal: Scientific Reports · Recommendation received: Accept after technical check

All five required corrections and all six minor comments are implemented. One of them
(RC4) turned out, on investigation, to concern the multivariable models themselves rather
than only the sample-accounting text, so it required re-analysis; we describe that in full
below.

---

## RC1 — GSE122482 is RNA-seq, not a microarray

**Confirmed and corrected.** The GEO series matrix for GSE122482 records
`!Series_type = "Expression profiling by high throughput sequencing"` and
`"Non-coding RNA profiling by high throughput sequencing"`, with
`!Series_platform_id = "GPL21290"`, and the series summary states that the 92,204 circRNAs
were obtained "through RNA-seq". Table 3 (formerly an unnumbered table in Section 4.1) now
reads "RNA-seq, Illumina HiSeq 3000 (GPL21290)", and the genome-build paragraph in Section
4.1 names the instrument as well. The words "circRNA expression microarray" no longer appear
anywhere in the manuscript.

## RC2 — References 3 and 5

**Both replaced with verified sources.** The reviewer is right that our previous
per-reference verification did not reach these two, and both had the same corruption
fingerprint as the round-1 "Isis A" entry.

| Ref | Was | Now |
|---|---|---|
| 3 | Jeck WR, Sharpless NE. "Detecting and removing non-templated introns from mature circular RNAs." *Cell* 2013;155(4):873–884 | **Jeck WR, Sharpless NE. Detecting and characterizing circular RNAs. *Nat Biotechnol* 2014;32(5):453–461. doi:10.1038/nbt.2890** |
| 5 | Li Z, Ruan K, Jin C, et al. "Circular RNA-mediated miRNA sponge promotes tumorigenesis." *Cancer Res* 2019;79(7):1652–1666; DOI unresolvable | **Chen G, Shi Y, Liu M, Sun J. circHIPK3 regulates cell proliferation and migration by sponging miR-124 and regulating AQP3 expression in hepatocellular carcinoma. *Cell Death Dis* 2018;9(3):306. doi:10.1038/s41419-017-0204-3** |

Ref 3 supports the Introduction claim about circRNA generation and detection; ref 5 supports
the claim that a circRNA can act as a miRNA sponge with a measured effect on tumour cell
behaviour, and we chose an HCC-specific study so that it also matches the disease context of
this manuscript. Ref 5's original DOI resolved to nothing, and searching the *Cancer Research*
79:1652–1666 range surfaced a different article with the adjacent DOI
(10.1158/0008-5472.CAN-18-1559), which is consistent with a transposed suffix.

We have now re-verified **all 22 references** against publisher records or PubMed, not only
the ones flagged, and confirmed that every one is cited in the text with continuous
numbering.

## RC3 — Abstract exceeds the word limit

**Corrected: 201 words** (previously 216 by the reviewer's count; our own earlier figure of
199 was wrong, and the Additional information block still said 195 — all three numbers were
inconsistent, which is exactly the kind of thing the reviewer cross-checks). The abstract
was rewritten without dropping any reported quantity. The Additional information block now
states 201 words. To prevent a recurrence we added a small script
(`count_abstract_words.py`) that reports the word count under three counting conventions,
so the figure quoted in the manuscript can always be regenerated from the file rather than
maintained by hand.

## RC4 — Sample accounting mixed tumour and non-tumour counts

**Confirmed, and the problem was worse than the reviewer described: it affected the models,
not only the text.**

What we found. The GDC clinical table for TCGA-LIHC contains 439 samples: 377 primary
tumours (sample type 01), 59 solid tissue normals (11) and 3 recurrent tumours (02). The
FPKM matrix has 424 columns, of which **53 are non-tumour**. Our export script
(`step14_export_for_R.py`) had no `sample_type` filter at all — it iterated the survival
table directly — so the multivariable models were fitted on 386 samples of which only 347
were primary tumours; 49 solid-tissue normal samples and 3 recurrent tumours had been
treated as tumour observations. The reviewer's observation that "399 exceeds the tumour
total" is exactly the symptom of this.

The reviewer's suspicion about the arithmetic was also correct and helped us locate it: the
"374 tumours" quoted in the abstract is 371 primary tumours plus the 3 recurrent samples
counted by our barcode rule, which treats any `0X` code as tumour.

**Correction applied.** We restricted the export to `sample_type_id = 01`, regenerated the
clinical export and re-fitted every multivariable model in R 4.6.1. The funnel is now stated
tumour-only in Section 4.1 and itemised in Supplementary Table S6:

**377 primary tumours → 371 with expression data → 365 with survival data → 347 with complete stage, age and sex → 341 also present in the expression matrix (117 events)**

The univariable Cox and log-rank analyses in Table 1 and Fig. 3 use 367 tumours (130 events).
For GSE14520 the corresponding counts are 247 annotated tumours → 225 with complete
clinical and survival data → 219 entering the models (84 events); as the reviewer requests,
the 445-sample expression matrix is now annotated as containing both tumour and non-tumour
samples, with only tumour samples entering any analysis.

**Effect on the reported multivariable results** (TCGA-LIHC only; GSE14520 is unaffected
because that cohort contained tumour samples only):

| Term | Previously reported | Corrected (n = 341, 117 events) |
|---|---|---|
| NT5DC2 | 1.219 (1.033–1.439), P = 0.019 | **1.388 (1.156–1.666), P = 0.00045** |
| HSD17B6 | 0.866, P = 0.061 | **0.847, P = 0.052** |
| BCLAF1 | 1.061, P = 0.504 | **1.100, P = 0.330** |
| SREBF2 | 0.977, P = 0.782 | **1.088, P = 0.378** |
| CD5 | 0.936, P = 0.456 | **0.883, P = 0.228** |
| LRCH1 | 0.909, P = 0.258 | **0.894, P = 0.235** |
| stage III/IV | 1.720, P = 0.0019 | **1.738, P = 0.0053** |
| sex (male) | 1.007, P = 0.970 | **0.960, P = 0.838** |

**The central conclusion is unchanged and in fact strengthened.** NT5DC2 remains the only
gene significant in both cohorts after adjustment (TCGA-LIHC P = 0.0004; GSE14520
P = 0.018), and its adjusted effect size is now slightly larger than its unadjusted value
rather than slightly smaller. Sex remains non-significant in both cohorts. HSD17B6 in TCGA
now sits just above the threshold (P = 0.052) and is reported as not significant, with the
boundary value stated explicitly.

One consequence deserves flagging because it reverses something we reported in round 2: with
tumour-only samples, **SREBF2's adjusted direction is now consistent across cohorts**
(TCGA-LIHC HR 1.09, GSE14520 HR 1.25). Our round-2 statement that it was "directionally
uncertain" was an artefact of the non-tumour samples entering the TCGA model. Section 2.2
and Table 1 have been updated accordingly, so the text and the table now agree.

## RC5 — Table inventory and response-letter consistency

**Corrected.** The file previously contained 11 physical tables against a declared inventory
of one main plus six supplementary. We have:

- **removed the duplicate adjusted-HR table** — Section 2.2 carried a six-gene adjusted-HR
  table that duplicated the adjusted columns of Table 1. The text now refers the reader to
  Table 1;
- **removed the duplicate ENCORI/Jaccard table** — it appeared both in Section 2.5 and as
  Supplementary Table S3. We kept S3, which carries the fuller column header and the
  accompanying note about ENCORI's coordinate convention;
- **numbered the three remaining in-text tables** that previously had no number or caption:
  Table 2 (quality stratification), Table 3 (data sources).

The inventory is now **7 figures and 3 tables (Tables 1–3), plus 6 supplementary tables
(S1–S6)**, which matches the physical file, and this is what the Additional information
block now states. Captions, in-text citations and physical order have been checked against
each other. The substantive alignment the reviewer identified — SREBF2 described as
"directionally uncertain" in the text but "Directionally consistent only" in Table 1 — is
resolved; both now read consistently, for the reason given under RC4.

**Alignment of this response letter with the final file.** Our round-2 letter stated that the
duplicate tables had been removed, which was not true of the file at that time, and quoted an
abstract word count of 199 that did not match the manuscript. Both statements are corrected
here: the duplicates are now genuinely removed, and the abstract is 201 words in both the
manuscript and this letter.

## Minor comments

1. **Section 3.1 denominator.** Now reads "34 of the 45 evaluable in both cohorts showed
   concordant direction in GSE14520 (76%)". The same paragraph was also updated for the
   corrected HSD17B6 value and for NT5DC2's strengthened effect size.
2. **Stray "COX" in Section 4.5** corrected to "Cox". A search of the final file returns no
   occurrence of "COX".
3. **"lying close to the diagonal"** removed. Section 2.2 now states the quantitative
   version: of the 45 genes evaluable in both cohorts, 34 fall on the same side of HR = 1 in
   each cohort.
4. **Fig. 4 legend** retitled "Lineage signature association matrix" and "21 immune cell
   types" changed to "21 lineage signatures", matching the terminology introduced in round 2.
5. **Introduction wording.** "proposed an indirect-support strategy to validate genes" now
   reads "proposed an internal-consistency grade for genes", which matches the evidential
   status set out in Section 4.7.
6. **Response-letter consistency** — see the end of RC5.

---

## Archive

The Zenodo archive has been updated to **v1.2.0**, which contains the corrected multivariable
Cox output, the corrected clinical export, the fixed export script, and the manifest
checksums regenerated over the new file set. The manuscript, supplementary information,
cover letter, README and CITATION file continue to cite the concept DOI
(10.5281/zenodo.23293107), which resolves to the most recent version, so no DOI in the text
needs updating when a new version is deposited.

---

We are grateful to the reviewer for the sample-accounting observation. It was the more
 consequential of the two data-handling problems found in this round, and the fact that it
surfaced as an arithmetic inconsistency rather than as a request for clarification is what
led us to inspect the export script itself rather than the sentence describing it.