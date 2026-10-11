# Response to Second-Round Reviewer Report

Manuscript: *NT5DC2 is an independently prognostic host gene of differentially expressed
circRNAs in hepatocellular carcinoma: a cross-cohort circRNA-miRNA-mRNA network analysis*
Journal: Scientific Reports · Recommendation received: Minor Revision

We thank the reviewer for the careful second-round assessment. The scientific concerns
raised in round 1 were addressed with new analyses, and we have now also resolved the
editorial and numerical items listed below. One of them — Major comment 4 — turned out to
be a genuine error in our code, not a presentational issue, and we describe it in full
because it changes several reported numbers.

---

## Major comment 4 — GSE14520 clinical variable parsing (the "all-female" cohort)

**The reviewer was right, and this was a coding error on our side.**

Our parsing code read the `Gender` column of `GSE14520_Extra_Supplement.txt` with the test
`value.lower() == 'male'`. That file codes sex as the single letters `M` and `F`, so the
test was never true and all 225 tumour samples were encoded as female. The reviewer correctly
identified that an entirely female cohort is implausible for HBV-related HCC.

**What we found on re-inspection of the source file:**

| | Value |
|---|---|
| Coding in the source file | `M` / `F` (single letter) |
| All 247 tumour samples | 211 M, 31 F, 5 missing |
| The 219 samples entering the multivariable models | 189 M, 30 F |

This is the sex distribution expected for an HBV-related cohort (85% male), which confirms
the diagnosis. The reviewer also asked us to check the stage column, as it comes from the same
file: stage is taken from `TNM staging` (values I, II, IIIA, IIIB, IIIC, IV) and grouped as
I/II versus III/IV. That parser was already guarded against missing values and we found no
defect there, but we have restated the coding explicitly in Section 4.4 and added
Supplementary Table S6 giving the sex and stage distribution of both cohorts.

**Correction applied.** We fixed the parser, regenerated the clinical export, and re-fitted
all multivariable models in R 4.6.1 with `survival::coxph` (Efron ties). Sex is now a
covariate in every model rather than being dropped. The affected numbers:

| Quantity | Previously reported | Corrected |
|---|---|---|
| GSE14520 NT5DC2 adjusted HR (95% CI), P | 1.271 (1.044–1.546), 0.0167 | **1.267 (1.041–1.542), 0.0184** |
| GSE14520 HSD17B6 adjusted HR, P | 0.800, 0.030 | **0.801, 0.033** |
| GSE14520 BCLAF1 adjusted HR, P | 1.186, 0.133 | **1.200, 0.114** |
| GSE14520 SREBF2 adjusted HR, P | 1.237, 0.058 | **1.247, 0.050** |
| GSE14520 CD5 adjusted HR, P | 0.807, 0.052 | **0.815, 0.069** |
| GSE14520 LRCH1 adjusted HR, P | 0.834, 0.134 | **0.846, 0.179** |
| GSE14520 stage III/IV adjusted HR, P | 3.237, 4.85 × 10⁻⁷ | **3.143, 1.20 × 10⁻⁶** |
| GSE14520 sex (male) term | not estimable | **HR 1.347 (0.643–2.822), P = 0.430** |

**The central conclusion is unchanged.** NT5DC2 remains the only gene significant in both
cohorts after adjustment (TCGA-LIHC P = 0.019; GSE14520 P = 0.018), and its effect size is
essentially unaffected. Sex is not a significant covariate in either cohort (TCGA P = 0.970;
GSE14520 P = 0.430), so no prognostic signal in this study is attributable to sex.

One new boundary case deserves our own flagging: **SREBF2 in GSE14520 now sits exactly at the
significance threshold (P = 0.050, 95% CI 1.00–1.55)**. We report it as not significant and
say so explicitly in Section 2.2 and in Table 1, rather than letting a rounded P value imply
otherwise.

---

## Major comment 1 — Residual numerical inconsistencies

| Item | Resolution |
|---|---|
| Section 2.2 "Effect direction" | **Corrected.** The paragraph previously assigned four genes to "poorer outcome", which contradicted both Table 1 and Section 2.3 for CD5. It now states that HSD17B6, CD5 and LRCH1 are protective and NT5DC2 and BCLAF1 are deleterious, and adds that SREBF2 is the one gene whose adjusted direction is *not* consistent between cohorts (TCGA HR 0.98, GSE14520 HR 1.25). This last point was not in the reviewer's list; we found it while checking the direction statements against the corrected model. |
| Fig. 6 legend edge count | **Corrected** from 3,960 to **4,540** (3,126 high + 1,414 medium), recomputed from the edge table. The node counts in the legend (37 circRNAs, 34 host genes) were also brought in line with the values printed inside the figure. |
| 95 matched vs 81 vs 79 circRNAs | **Re-derived from the data, and the reviewer's arithmetic identified a stale number.** Our own counts are 34 exact + 15 loose + 30 gene-fallback = 79, and 103 − 79 = 24 candidates absent from the network. The "95 matched entries" figure that appeared in the previous version was a leftover from an earlier matching pass and could not be reproduced from any current result table; we removed it rather than reconcile it. Section 2.6 now gives the three-tier breakdown with the arithmetic made explicit, and states that the 24 absent candidates are missing because circBank declares either no comparable entry or no miRNA for any entry it declares. |
| Fig. 7c counts 4/1/6 | **Addressed in the legend:** genes are counted once per tier, so a gene with edges in more than one tier appears in each, and the three values need not sum to six. |
| Top eight vs sixteen hub miRNAs | **Reconciled explicitly.** These are two different rankings: Supplementary Table S1 ranks hubs over the complete network (all tiers), whereas Fig. 6 ranks within the high and medium tiers only. Section 2.4 now says so and no longer implies they are the same set. |
| 45 evaluable genes | **Defined.** Of the 87 host genes analysed in TCGA-LIHC, 45 carry a probe on the GPL3921 array and can be evaluated in GSE14520. We verified the platform annotation directly: **none of the 42 remaining genes has any probe on that array**, so their absence is a platform coverage limitation rather than missing data. Section 2.2 and the Fig. 3 legend now state this and give the denominator (34/45). |
| Duplicate tables and counts | **Resolved.** Only one adjusted-HR table (Table 1) and one ENCORI/Jaccard table now remain; the duplicates have been removed. "Additional information" now reads 7 figures, 1 main table and 6 supplementary tables, and 22 references. |
| Section 4.4 verification sentence | **Rewritten.** The comparison is now stated as like-for-like on identical sample sets, and the three disagreements are attributed to P values between 0.03 and 0.09 changing category when the analysis set narrowed from 225 to 219 samples. |

## Major comment 2 — Leftover first-draft material and references

The "must-fix before submission" block has been deleted. On its content: the note asked for
reference 24 to be completed, citing GSE122482 with PMID 30537115. **There is no reference
24 to add** — PMID 30537115 is Qiu and colleagues 2019, which is already reference 19 in
this manuscript and is the study that generated GSE122482. Section 4.1 cites it as such. The
block appears to have been an earlier plan item that outlived its purpose.

Reference corrections, each verified against the publisher record:

| Ref | Correction |
|---|---|
| 4 | Was listed as Enuka TF, Beukman CJ, Aoki M, *Genes Dev* 2016, "Circular RNAs are a large class of animal RNAs…". The DOI resolved to an unrelated article. Now: Enuka Y, Lauriola M, Feldman ME, Sas-Chen A, Ulitsky I, Yarden Y. "Circular RNAs are long-lived and display only minimal early alterations in response to a growth factor." *Nucleic Acids Res* 2016;44(3):1370–1383. doi:10.1093/nar/gkv1367 |
| 8 | Author "Isis A" was a mangled entry and the journal/year were wrong. Now: Salmena L, Poliseno L, Tay Y, Kats L, Pandolfi PP. "A ceRNA hypothesis: the Rosetta Stone of a hidden RNA language." *Cell* 2011;146(3):353–358. doi:10.1016/j.cell.2011.07.013 |
| 11 | Second author corrected from "Mohn F" to Mlecnik B, per *Genome Biol* 2013;14(7):R115. |
| 21 | Duplicated bibliographic entry and the Chinese-language annotation removed; the DOI 10.7150/jca.45890 does not resolve, so the entry now cites the verifiable PubMed identifier 34094907. |
| 22 | Repeated "Liu S" and the "et al." after nine authors removed; now Zhou K, Huang J, Liu S, et al. |
| Section 4.6 | circBank citation corrected from ref 17 to **ref 15**. |

All 22 references are cited in the text and the numbering is continuous; we verified this
programmatically rather than by eye.

## Major comment 3 — Methods details

| Item | Resolution |
|---|---|
| Genome build | **Added, with an honest caveat.** Neither GSE122482 nor circBank states its build in the files we retrieved, so we could not verify that they agree. We now say so explicitly, state that **no liftOver was performed** and that matching relied on chromosome-prefix normalisation alone, and identify an unverified build difference as one possible explanation both for the residual unmatched candidates and for the MGC50722/ADAMTS13 annotation conflict. |
| Expression preprocessing | **Added.** GSE122482: per-sample circRNA RPM tables from the GEO supplementary files, used as delivered, no further normalisation. GSE14520: probe-level values from the GEO series matrix as delivered (log-transformed by the submitters; observed value range 3.1–12.2 confirms the log scale), collapsed to gene level by averaging probes per symbol, no batch correction as the cohort is a single series. TCGA-LIHC: FPKM from the GDC harmonised STAR-Counts data set, mapped from Ensembl gene IDs, log2(FPKM + 1). All models use z-scores. |
| GSE122482 platform | **Stated:** GPL21290, with the series title and the 20 per-sample RPM tables described. |
| Section 4.3 | **Rewritten.** The previous sentence conflated two things. It now reports that 48,109 of 92,204 records (52.2%) carry a gene assignment spanning 11,048 distinct genes, and separately that the 103 candidates map to 91 host genes, of which nine host more than one candidate circRNA. |
| Section 4.6 "chr12" | Corrected to `chrN`. |

## Major comment 5 — Availability statements and archive versioning

- **Code availability** rewritten. It previously claimed the Python standard library was used
  throughout, which contradicted Section 4.9. It now describes the two-stage pipeline: Python
  3.13.12 standard library for data handling, R 4.6.1 with `stats` and `survival` for all
  statistical inference, with `Pillow` and a headless browser needed for figures only.
- **Zenodo archive updated** to v1.1.1, containing the corrected multivariable Cox output,
  the full 92,204-feature FDR table, the ENCORI comparison and the R scripts. The manifest
  checksums were regenerated. The manuscript, supplementary information, cover letter,
  README and CITATION file now cite the concept DOI (10.5281/zenodo.23293107), which
  resolves to the most recent version and therefore needs no further updating.
- **TCGA-LIHC 374 vs 386** reconciled in Section 4.1 and itemised in Supplementary Table S6:
  424 samples with expression data → 399 with survival and covariate data → 386 entering
  multivariable models (145 events). The attrition reflects missing clinical fields, not
  missing expression.

## Minor comments

1. **ENCORI access date** changed from the future date 11 November 2026 to **11 October 2026**
   in both places where it appeared.
2. **Ref 22** verified and corrected as described above. We did not add a separate starBase
   citation: the ENCORI resource is described by the same paper, and citing it twice would
   duplicate one source.
3. **Abstract** now reads "34 of the 45 evaluable in both cohorts showed concordant direction
   (76%); six of these met the combined validation criterion", and the FDR sentence now gives
   both denominators: 3 features at FDR < 0.05 across all 92,204 tested features, of which 2
   belonged to the 103 candidates, against 5,649 at uncorrected P < 0.05. The abstract is 199
   words.
4. **Title** changed. The plural form overstated the evidence once only NT5DC2 survived
   adjustment. The new title foregrounds NT5DC2 while retaining the network-analysis
   descriptor.
5. **"chromatic harmonisation"** → **"chromosomal harmonisation"**.
6. **r-value notation** made consistent: two-cohort values are given as a pair throughout
   (e.g. CD5 with Treg cells r = 0.812 in TCGA-LIHC and 0.463 in GSE14520; CD5 with CD4+ T
   cells r = 0.877 and 0.378).
7. **"COX"** → **"Cox"** throughout (13 occurrences); **"Charoentaneath"** → **"Charoentong"**
   (2 occurrences). The "close to the diagonal" phrasing in Section 2.2 is replaced by the
   explicit count of 34 genes with concordant direction out of 45 evaluable. The Fig. 5 legend
   now reports the non-significant TCGA log-rank for HSD17B6 (P = 0.067) and adds a sentence
   reconciling CD5's log-rank (0.0019) with its Cox P (0.063) results in terms of what the two
   tests measure. On the MGC50722/ADAMTS13 annotation conflict: Section 2.1 already explains it
   as a transcript-level annotation difference between resources, and we now cross-reference
   the genome-build caveat added to Section 4.1 so the two explanations are visible together.
8. **"ssGSEA-style"** removed. The score is now described accurately as a rank-average
   signature score, with the statement that no ECDF-based enrichment weighting is applied and
   that we therefore do not describe it as ssGSEA; the citation now points to Bindea and
   colleagues rather than GSVA. We also state that the 27 sets are **lineage** signatures and
   not exclusively immune: they include endothelial cells, cancer-associated fibroblasts and
   fibroblasts, and we use "lineage signatures" throughout.

---

## Summary of changes to reported numbers

Only one set of numbers changed as a result of re-analysis, all traceable to the corrected
sex parsing: the six GSE14520 adjusted hazard ratios, the GSE14520 stage hazard ratio, and
the addition of the GSE14520 sex term. TCGA-LIHC results are unchanged. All derived counts
quoted in the text (4,540 edges in Fig. 6, the 79 circRNAs and their 34/15/30 tier
breakdown, the 45 evaluable genes) were recomputed from the archived result tables rather
than carried over.

We believe the manuscript now addresses every point raised. Thank you again — the
all-female observation was decisive, and we would not have found that parsing error without
it.