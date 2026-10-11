# Host genes of differentially expressed circRNAs show consistent cross-cohort prognostic value in hepatocellular carcinoma

**Wanli He**<sup>1*</sup>

<sup>1</sup> School of Medicine, Shangqiu Institute of Technology, Shangqiu, China

\* Correspondence: hewanli66@126.com; ORCID 0009-0002-0829-1909

---

## Abstract

Hepatocellular carcinoma (HCC) shows substantial prognostic heterogeneity. We screened GSE122482 (10 pairs of HBV-related tumour and peritumour tissues), identified 103 differentially expressed circRNAs, and evaluated the prognostic value of their **host-gene mRNAs** in TCGA-LIHC (374 tumours) and GSE14520 (225 tumours). The circRNA–miRNA network comprised 15,479 deduplicated edges, 2,378 miRNAs, 79 circRNAs and 73 host genes. Of the host genes, 37 were significant in TCGA, 34 showed concordant direction across cohorts, and six passed independent validation. In multivariable Cox models adjusting for AJCC stage, age and sex, **NT5DC2 remained significant in both cohorts** (adjusted HR 1.22, P = 0.019; 1.27, P = 0.017); HSD17B6 reached significance in GSE14520 (0.80, P = 0.030) but not TCGA (P = 0.061); the remaining four were not significant after adjustment. Because neither cohort contains circRNA data and GSE122482 provides no matched mRNA quantification, **the markers identified are host-gene mRNA markers; circRNA-level prognostic value remains untested**. Benjamini–Hochberg correction across 92,204 features yielded 2 surviving candidates (FDR < 0.05) versus 5,649 under uncorrected P < 0.05. Comparison with ENCORI/starBase showed minimal overlap between prediction sources (Jaccard 0.000–0.035), indicating that connectivity-based support reflects internal consistency of one prediction resource rather than independent validation.

**Keywords**: hepatocellular carcinoma; circRNA; host gene; prognostic marker; immune infiltration; network analysis

---

## 1. Introduction

Hepatocellular carcinoma (HCC) is one of the most common malignancies worldwide, with particularly high incidence in populations chronically infected with hepatitis B virus<sup>1</sup>. Despite advances in surgical resection, interventional therapy and systemic treatment, the five-year survival rate of HCC has remained disappointingly low, largely because of late diagnosis and treatment resistance<sup>2</sup>. Identifying stable biomarkers that can assess prognosis preoperatively therefore retains clear clinical value.

Circular RNAs (circRNAs) are a class of RNA molecules characterised by a covalently closed loop structure, generated from linear precursor RNA through back-splicing<sup>3–5</sup>. Owing to the resistance of this closed structure to exoribonucleases, circRNAs are more stable than their linear counterparts and are readily detectable in tissue and plasma<sup>6</sup>. Mechanistically, circRNAs can act as miRNA sponges, sequestering miRNAs and thereby releasing the translation of their target mRNAs, forming a circRNA–miRNA–mRNA axis. This axis has been implicated in cell proliferation, apoptosis, invasion and immune escape<sup>7,8</sup>.

Two limitations characterise much of the existing literature. First, most studies focus on a single circRNA–miRNA axis rather than performing an unbiased transcriptome-wide scan, which risks missing genuinely prognostic candidates. Second, the circRNA prediction tools used in mechanistic work (such as circBank) return very large and largely unfiltered binding sets—a single circRNA may yield predictions for hundreds of miRNAs, and incorporating these indiscriminately introduces substantial false positives that dilute the true signal<sup>9</sup>.

A third issue is routinely overlooked: coordinate systems differ between data sources. circRNA identifiers commonly follow the `chrN:start-end` convention, whereas databases such as circBank use `N:start-end`; if these are not harmonised, coordinate matching fails catastrophically.

We therefore adopted a three-step strategy of screening, validation and stratification. First, we identified differentially expressed circRNAs in paired tumour and peritumour tissues. Second, we cross-validated the prognostic value of their host genes in two independent cohorts. Third, we integrated circBank miRNA binding predictions together with binding free energies, stratified each edge by matching confidence and binding strength, and proposed an indirect-support strategy to validate genes for which coordinate matching fails.

---

## 2. Results

### 2.1 Overview of differentially expressed circRNAs

Paired comparison of GSE122482 across all 92,204 circRNA features identified 103 differentially expressed circRNAs (59 up-regulated, 44 down-regulated) at an uncorrected P < 0.05, corresponding to 91 host genes (Fig. 1). Nine genes corresponded to more than one circRNA (CYP4F23P, GALP, HIST2H4B, LMNA, MEG3, MGC50722, RBM19, TMEM179, ZNF250), indicating diversity in circRNA generation at these loci.

**Multiplicity correction changes the picture substantially.** Benjamini–Hochberg correction across the 92,204 tested features retains only 3 circRNAs at FDR < 0.05 (8 at FDR < 0.10), whereas 5,649 features satisfy the uncorrected P < 0.05 threshold—a difference of two orders of magnitude that quantifies the false-positive burden of uncorrected screening with ten sample pairs. Within the 103 candidates, only **2 passed FDR < 0.05** (FDR range 0.0449–0.5156, median 0.5156): `hsa_circ:chr12:56754955-56755968` (HSD17B6, log2FC = −1.139, P = 8.25 × 10⁻⁷, FDR = 0.0449) and `hsa_circ:chr9:136302869-136303486` (MGC50722/circADAMTS13, log2FC = −1.093, P = 1.19 × 10⁻⁶, FDR = 0.0449). The full FDR-ranked candidate list is provided in Supplementary Table S4. We therefore treat the 103-circRNA set as an exploratory screen and place the inferential weight of the study on the cross-cohort survival analysis rather than on the differential-expression step.

One candidate, hsa_circ:chr9:136302869-136303486, requires comment. This coordinate corresponds to circADAMTS13 (hsa_circ_0005075), which Qiu and colleagues identified as a miR-484 sponge in the very dataset used here<sup>19</sup>. Our annotation assigned this locus to a different host gene (MGC50722), reflecting the differing annotation conventions between resources (Section 3.2); circBank names circRNAs by host gene and assigns this entry to ADAMTS13. Importantly, our data agree with theirs: we likewise find this circRNA down-regulated in tumour tissue (log2FC = −1.09, P = 1.19 × 10⁻⁶), and hsa-miR-484 appears in our network as a predicted partner. This concordance provides an external check on our pipeline: circADAMTS13 functions here as a positive control rather than a discovery, and the remaining candidates are the underexplored space that the cross-cohort design is intended to expose (Section 3.2).

### 2.2 Prognostic value of host genes and cross-cohort validation

**The survival analyses below were performed on linear host-gene mRNAs, not on circRNAs.** Neither TCGA-LIHC nor GSE14520 contains circRNA expression data. GSE122482 provides circRNA read counts but, as we verified directly, contains no matched mRNA quantification: its GEO series matrix holds metadata only, and the supplementary files comprise one TAR archive and 20 circRNA RPM tables. The correlation between each candidate circRNA and its host-gene mRNA within the same samples therefore could not be assessed, and we did not assume it. The prognostic markers identified here are accordingly host-gene-level mRNA markers; whether the corresponding circRNAs carry prognostic information at the circRNA level remains to be determined by direct measurement.

In TCGA-LIHC, the expression of 37 host genes was significantly associated with overall survival (COX regression P < 0.05), of which 25 were significant by both COX regression and the log-rank test.

Independent evaluation in GSE14520 showed concordant direction for 34 genes (76%). Applying the criterion of concordant direction together with COX regression P < 0.05 in GSE14520, six genes passed validation (Table 1).

**Validation strength was heterogeneous.** NT5DC2 and HSD17B6 were significant by COX regression in both cohorts and represent the strongest evidence in this study. BCLAF1, SREBF2, LRCH1 and CD5 did not reach significance by COX regression in TCGA (P = 0.086, 0.201, 0.092 and 0.063 respectively); these were included on the basis of significant COX regression in the validation cohort with concordant direction, and CD5 was additionally significant by the log-rank test in TCGA (P = 0.0019). Figure 3 shows the scatter distribution of hazard ratios in the two cohorts, with 34 paired genes lying close to the diagonal.

**Effect direction.** Higher HSD17B6 and LRCH1 expression was associated with better overall survival in both cohorts (protective), whereas the remaining four genes were associated with poorer outcome (Fig. 5).

**Multivariable analysis identifies only NT5DC2 as an independently prognostic marker.** Because prognostic significance in TCGA alone could reflect confounding by clinical covariates, we fitted Cox models adjusting for AJCC stage (I/II vs III/IV), age (centred at 50 years) and sex, implemented in R with `survival::coxph` (Efron handling of ties). In TCGA-LIHC (n = 386, 145 events) and GSE14520 (n = 219, 84 events; all tumour samples are female, so sex was non-estimable and excluded) the adjusted results were:

| Gene | TCGA adjusted HR (95% CI), P | GSE14520 adjusted HR (95% CI), P | Independent in both cohorts |
|---|---|---|---|
| **NT5DC2** | **1.22 (1.03–1.44), P = 0.019** | **1.27 (1.04–1.55), P = 0.017** | **Yes** |
| HSD17B6 | 0.87 (0.75–1.01), P = 0.061 | 0.80 (0.65–0.98), P = 0.030 | No (TCGA borderline) |
| BCLAF1 | 1.06 (0.89–1.26), P = 0.504 | 1.19 (0.95–1.48), P = 0.133 | No |
| SREBF2 | 0.98 (0.83–1.15), P = 0.782 | 1.24 (0.99–1.54), P = 0.058 | No |
| CD5 | 0.94 (0.79–1.11), P = 0.456 | 0.81 (0.65–1.00), P = 0.052 | No |
| LRCH1 | 0.91 (0.77–1.07), P = 0.258 | 0.83 (0.66–1.06), P = 0.134 | No |

Stage was a strong prognostic factor in both cohorts (adjusted HR 1.72, P = 0.0019 in TCGA; 3.24, P = 4.8 × 10⁻⁷ in GSE14520), confirming that adjustment was informative. **NT5DC2 is the only gene remaining significant in both cohorts after adjustment**, and its effect size changed little (unadjusted 1.25 to adjusted 1.22 in TCGA; 1.37 to 1.27 in GSE14520), indicating that its prognostic association is not attributable to stage, age or sex. BCLAF1, SREBF2, CD5 and LRCH1 lost significance after adjustment and should be described as directionally consistent only (see Section 3.4).

### 2.3 Immune infiltration associations

Applying the retention criterion stated in Methods (directionally concordant across cohorts with |r| > 0.15 in both), **five genes** showed stable immune associations:

- **CD5 was strongly positively correlated with multiple T-cell populations**, satisfying the criterion for 19 gene–cell pairs (Treg cells r = 0.812/0.463, CD4+ T cells r = 0.877, Tfh r = 0.870, effector CD4 r = 0.863, CD8 r = 0.816; first value TCGA, second GSE14520)
- **LRCH1** with Tfh cells (r = 0.258/0.262), **HSD17B6** with endothelial cells (r = 0.160/0.337) and memory CD8 cells (r = −0.243), **SREBF2** with NK cells (r = −0.229/−0.302) and **SMG5** with endothelial cells (r = −0.200/−0.296)

Four further genes were directionally concordant but were **excluded**, for reasons that differ and are worth stating explicitly, because two of them had appeared in an earlier version of our analysis. **BCLAF1** showed the strongest association of any gene in TCGA (Th17 cells, r = +0.466) but was excluded because the two cohorts were **opposite in direction** for all 11 of its nominally strong pairs (Th17: +0.466 in TCGA versus −0.280 in GSE14520). **NT5DC2** and **NUP205** were likewise excluded for direction reversal (NT5DC2, neutrophils: +0.158 versus −0.223). **RBM19** was excluded because its strongest association reached only r = 0.088, below the retention threshold. CD5 is **not** a member of any of the immune feature gene sets used here, so its correlation with T-cell scores cannot be a circular consequence of shared membership.

The dual character of CD5 is noteworthy: it was significant by COX regression in GSE14520 and by the log-rank test in TCGA, while its high expression correlated positively with several T-cell populations. Notably, CD5 has HR < 1 in both cohorts, indicating that it is a **protective** factor (Fig. 4).

### 2.4 The circRNA–miRNA–mRNA network

Integration with circBank predictions yielded 15,479 circRNA–miRNA edges involving 2,378 miRNAs, 79 circRNAs and 73 host genes. The quality stratification is summarised in Fig. 7:

| Tier | Edges | miRNAs | circRNAs | Genes |
|---|---|---|---|---|
| high | 3,126 | 1,317 | 34 | 32 |
| medium | 1,414 | 851 | 15 | 14 |
| low | 10,939 | 2,307 | 79 | 73 |

Ranked by the number of connected circRNAs, the top eight hub miRNAs connect 19–22 distinct circRNAs and 18–21 distinct host genes, with mean binding free energies mostly below −22 kcal/mol, indicating relatively strong predicted binding (Supplementary Table S1; Fig. 6).

### 2.5 Internal consistency support for genes matched by fallback

For 26 genes, circRNA entries had to be matched by same-gene, same-chromosome fallback, leaving the coordinate association uncertain. We applied an **internal consistency** heuristic: if a given miRNA also connects multiple exactly matched genes, its binding prediction is endorsed consistently across independent loci within the same prediction resource. Grading: ≥ 5 exact genes is strong, 2–4 moderate, 1 weak, 0 none. For **NT5DC2**—the only fallback-matched gene that also passed cross-cohort validation—28 of its 357 unique miRNAs reached the strong tier and 151 the moderate tier, giving 50.1% with strong or moderate internal support (counts are per unique miRNA).

**This heuristic does not constitute independent validation, and we tested that directly.** All edges in our network originate from a single prediction resource (circBank), whose human entries carry no quality annotations (the `highConfidence`, `validated` and `conserved` fields are empty for all 95 matched records). If circBank's seed matching is permissive, a miRNA that is broadly permissive will connect many genes by construction, and connectivity would then reflect the predictor rather than biology. We compared our predictions with ENCORI/starBase, which integrates seven independent algorithms (PITA, RNA22, miRmap, microT, miRanda, PicTar, TargetScan) together with CLIP-seq and degradome-seq evidence (accessed 11 November 2026):

| Gene | circBank miRNAs | ENCORI miRNAs | Shared | Jaccard |
|---|---|---|---|---|
| NT5DC2 | 357 | 17 | 3 | 0.008 |
| BCLAF1 | 16 | 142 | 3 | 0.019 |
| SREBF2 | 194 | 129 | 11 | 0.035 |
| LRCH1 | 46 | 48 | 0 | 0.000 |
| HSD17B6 | 297 | 0 | 0 | – |
| CD5 | 373 | 0 | 0 | – |

The overlap is minimal (Jaccard 0.000–0.035), and circBank assigns substantially more miRNAs per gene than ENCORI; for HSD17B6 and CD5 the independent resource returns no prediction at all. **The two resources therefore disagree not only in degree but largely in substance about which miRNA targets which gene**, and the great majority of the 15,479 edges in our network are not corroborated by an independent predictor.

**Consequently we withdrew the mechanistic candidate previously highlighted from this section.** hsa-miR-4706 had been nominated on the basis of a predicted binding energy of −33.3 kcal/mol on NT5DC2 and connectivity to eight exactly matched genes. ENCORI's seven algorithms do not predict hsa-miR-4706 as targeting NT5DC2, and no CLIP-seq or degradome-seq evidence supports the interaction. Both original grounds for nomination are therefore withdrawn, and the NT5DC2–hsa-miR-4706–CD5/SREBF2 axis is not asserted here. We retain the connectivity data as a description of the prediction resource and present the internal-consistency grading solely as a prioritisation heuristic for candidate selection.

### 2.6 Two technical issues affecting data reliability

Two issues were identified during the analysis and are reported transparently.

**Duplicate edges from pagination overlap in circBank.** The binding interface returns overlapping records at pagination boundaries, producing 831 duplicate edges (5.1%) in the raw retrieval. Deduplication on the composite key of circRNA, host gene, miRNA, energy and score was applied throughout; all figures reported here derive from the completed, deduplicated dataset of 15,479 edges.

**Discrepancy in chromosome naming conventions.** GSE122482 uses `chrN` whereas circBank uses `N`; these cannot be compared directly. Without harmonisation, only 33 of 103 candidate circRNAs could be matched to circBank entries at exact coordinates (32%). Normalising chromosome names raised the matching rate to 92% (exact 33, loose 18, gene_fallback 44). This was the single largest technical determinant of matching success in this study.

**Incomplete retrieval of circBank binding records.** Our initial retrieval covered 1,925 of 2,448 result pages (79%), and inspection of the coverage against the database showed that five circRNA entries with declared miRNA predictions had been missed, among them the entries for HSD17B6 and for circADAMTS13. After targeted re-retrieval (976 additional edges, with declared and retrieved counts agreeing entry by entry) the network grew from 14,568 to 15,479 edges, from 74 to 79 circRNAs and from 69 to 73 host genes. All figures reported here are derived from the completed dataset. A further 14 of the 95 matched entries carry no miRNA prediction in circBank at all and are therefore legitimately absent from the network.

---

## 3. Discussion

### 3.1 Main findings

This study makes three principal contributions, each with a correspondingly bounded claim.

**First, panoramic screening with cross-cohort survival validation identifies host genes of prognostic interest.** Systematic screening of 103 circRNAs identified six host genes with cross-cohort prognostic relevance, and multivariable adjustment reduced these to a single robust candidate. A conventional "one axis at a time" strategy would have been likely to miss these, particularly NT5DC2 and HSD17B6, which have received little attention at the circRNA level in the existing literature.

**Second, cross-cohort validation and multivariable adjustment markedly improve reliability, and the residual evidence is graded rather than uniform.** Of 37 prognostic genes in TCGA, 34 showed concordant direction in GSE14520 (76%), and six met the combined criterion of concordant direction plus significant COX regression in the validation cohort. After adjustment for AJCC stage, age and sex, **NT5DC2 remained significant in both cohorts** with an essentially unchanged effect size, making it the only independently prognostic host gene in this study. HSD17B6 retained significance in GSE14520 but became borderline in TCGA (P = 0.061); BCLAF1, SREBF2, CD5 and LRCH1 lost significance altogether. We accordingly distinguish a single robust candidate (NT5DC2), one partially supported candidate (HSD17B6) and four directionally consistent genes, rather than reporting all six as validated markers.

**Third, annotation quality is a first-order problem, and we were able to quantify how far it compromises inference.** circBank does not populate the commonly used `highConfidence`, `validated` and `conserved` fields for human entries (100% empty in our 95 matched entries); coordinate systems differ systematically between sources; our own initial retrieval was incomplete. We addressed these by match-type stratification, binding-energy filtering, chromatic harmonisation and completed re-retrieval, and we report the residual uncertainty explicitly. However, comparison with an independent resource (ENCORI/starBase, seven algorithms plus CLIP-seq and degradome-seq evidence) showed a Jaccard overlap of only 0.000–0.035 with our predictions. **The stratified framework therefore improves the internal coherence of a single prediction resource; it does not substitute for independent validation of the predicted interactions.**

### 3.2 Comparison with existing studies

Several circRNA databases are available, each with distinct annotation strategies and nomenclature<sup>15–18</sup>, which introduces heterogeneity that complicates cross-study comparison. Existing HCC circRNA research has focused on a small number of well-studied molecules and has typically validated a single miRNA under a single ceRNA hypothesis. Our methodological differences are threefold: we start from data rather than from a predefined molecule; we introduce binding free energy as a filtering dimension rather than relying on prediction counts alone; and we apply explicit quality stratification to predicted relationships instead of treating thousands of predictions as equivalent.

The present dataset has been used in two earlier analyses, and we disclose both. Qiu and colleagues generated these data and reported circADAMTS13 as a miR-484 sponge<sup>19</sup>; our reuse of the same data was not intended to revisit that circuit. Yu and colleagues applied the same dataset to WGCNA-based ceRNA network construction with RT-qPCR validation of five networks<sup>20</sup>. Neither study performed prognostic evaluation of circRNA host genes across independent cohorts, which is the principal addition here.

An important consequence concerns the interpretation of our own predictions. Of the circRNAs we identified, one—circADAMTS13—has already been experimentally characterised, including confirmation of its miR-484 sponge activity by dual-luciferase assays with a mutated target site<sup>19</sup>. Our pipeline correctly recovered both its direction of change and its partner miRNA, which is reassuring, but it also means that this entry is a positive control rather than a discovery. We were unable to cross-match the remaining published circRNA identifiers against our candidate list, because the identifiers used by the source studies were not resolvable to genomic coordinates in the resources we used. We therefore make no claim about how many of our candidates remain experimentally untested; determining this requires a systematic comparison against the primary literature, which we have not performed. What can be stated is that the six genes we identify as cross-cohort validated, and NT5DC2 in particular, have not to our knowledge been examined at the circRNA level in HCC, and they constitute the underexplored space that our cross-cohort design is designed to expose.

A related clarification concerns the level at which our own findings should be interpreted. Because the survival analyses were performed on host-gene mRNAs rather than on circRNAs, and because neither validation cohort contains circRNA expression data, our data establish that the **host genes** of certain differentially expressed circRNAs have prognostic value. They do not establish that the circRNAs themselves do. Host-gene mRNA abundance is not a reliable proxy for circRNA abundance, since back-splicing efficiency is regulated independently of overall transcription, and the matched measurements needed to test that proxy assumption are not available in GSE122482. This is a limitation of the available data rather than of our analysis, but it defines the scope of the claim: the first link of the chain (differential circRNA expression) and the second (host-gene prognosis) are supported by data, whereas the third (circRNA-level prognosis) remains untested.

### 3.3 Boundaries of the immune mechanism

The immune-related conclusions require careful delimitation, and applying the stated retention criterion strictly changes which genes qualify. **CD5 is the only gene with a strong, reproducible immune association**: it satisfies the criterion for 19 gene–cell pairs, and its correlation with CD4+ T-cell infiltration (r = 0.877 in TCGA) is corroborated in GSE14520. CD5 is not a member of any of the immune feature gene sets used here, so this correlation cannot arise from shared gene membership.

The remaining four retained genes show associations of |r| < 0.34 (LRCH1 with Tfh cells, HSD17B6 with endothelial and memory CD8 cells, SREBF2 with NK cells, SMG5 with endothelial cells). Genes that appeared in an earlier version of this analysis have been removed: **BCLAF1**, whose TCGA association with Th17 cells (r = +0.466) reverses in GSE14520 (−0.280) and was therefore not directionally stable; **NT5DC2** and **NUP205**, likewise excluded for direction reversal; and **RBM19**, whose strongest association (r = 0.088) falls below the retention threshold. Immune mechanisms are therefore discussed only in relation to CD5, without extrapolation, and immune infiltration is not invoked to explain the prognostic significance of the other candidates.

### 3.4 Limitations

**First, and most fundamentally, the validated markers are host-gene mRNAs rather than circRNAs.** Neither TCGA-LIHC nor GSE14520 contains circRNA expression data, and GSE122482 provides circRNA read counts without matched mRNA quantification, so the correlation between each candidate circRNA and its host-gene mRNA within the same samples could not be computed. Host-gene mRNA abundance is not a reliable proxy for circRNA abundance, because back-splicing efficiency is regulated independently of overall transcription. CircRNA-level prognostic value therefore remains to be established by direct measurement.

**Second, differential screening was exploratory and most candidates do not survive multiplicity correction.** The screening threshold was an uncorrected P < 0.05 across 92,204 features with only ten sample pairs. We have since applied Benjamini–Hochberg correction (Results 2.1): 5,649 features satisfy the uncorrected threshold but only 3 reach FDR < 0.05, and only 2 of the 103 candidates do so. The statistical power of ten paired samples is inherently limited, and the FDR results should be regarded as the more defensible statement of what the screen established.

**Third, predicted circRNA–miRNA interactions are not independently corroborated.** All edges derive from circBank, which supplies no quality flags for human entries and predicts many more miRNAs per circRNA than independent tools. Comparison with ENCORI/starBase (seven algorithms plus CLIP-seq and degradome-seq evidence) yielded a Jaccard overlap of 0.000–0.035, and the leading candidate interaction previously highlighted is not supported by any of those algorithms. The connectivity-based grading we report is an internal consistency heuristic, not evidence of interaction, and a hub-like miRNA is expected to score highly by construction. Establishing any proposed axis requires dual-luciferase reporter assays with a mutant control in which the miRNA seed region is deleted, so that non-specific suppression can be excluded.

**Fourth, coordinate association is uncertain for 26 fallback-matched genes.** For these genes the circRNA–miRNA edge lacks reliable coordinate assignment, and the internal-consistency grading applied to them inherits all the limitations of the third point.

**Fifth, immune infiltration scores are indirect estimates.** They were derived as rank-transformed means of feature gene sets; their absolute magnitudes are not comparable to direct measurements such as flow cytometry or immunohistochemistry, and the rank transform may underestimate true correlation<sup>12</sup>. The gene sets are representative markers selected for liver tissue rather than a curated multi-tissue set.

**Sixth, this is a purely computational analysis.** All conclusions derive from secondary analysis of public data without in vitro or in vivo validation. Follow-up work should confirm circRNA expression differences between tumour and peritumour tissue by qPCR with junction-spanning primers, and test direct circRNA–miRNA interaction by dual-luciferase reporter assay.

### 3.5 Future directions

Three priorities follow from the present analysis, each stated at the level the data actually support. First, **direct measurement of the candidate circRNAs themselves**, by qPCR with back-splicing-junction-spanning primers in paired tumour and peritumour tissue; this is the only route to testing circRNA-level prognostic value, and it is a prerequisite for any mechanistic work. Second, **investigation of how CD5 influences prognosis through CD4+ T-cell infiltration**, the one immune association that is both strong and reproducible. Third, **exploration of HSD17B6 as a candidate loss-of-function protective factor**, noting that its prognostic contribution is significant in GSE14520 but borderline in TCGA and requires confirmation in a further cohort.

We deliberately do not nominate a circRNA–miRNA axis for experimental follow-up. The connectivity evidence previously used for that purpose is internal to a single prediction resource and was not corroborated by the independent algorithms we consulted (Results 2.5); a nomination made on that basis would not survive scrutiny. Candidate selection for future mechanistic work should instead begin from predictions that are concordant across independent resources.

---

## 4. Methods

### 4.1 Data sources

| Dataset | Platform / type | Sample size | Use |
|---|---|---|---|
| GSE122482 | circRNA expression microarray | 10 pairs, HCC tumour / peritumour | Differential circRNA screening |
| TCGA-LIHC | RNA-seq + clinical | 374 tumours | Prognostic modelling (discovery cohort) |
| GSE14520 | GPL3921 + clinical | 225 tumours | Independent validation cohort |
| circBank | database | 828,808 human circRNAs | circRNA annotation and miRNA binding prediction |
| ENCORI/starBase | database | 7 prediction algorithms + CLIP-seq/degradome-seq | Independent check of miRNA predictions |

Clinical data for GSE14520 were obtained from its Extra Supplement file, and only tumour samples were included. Immune cell feature gene sets were taken from the pan-cancer immunogenomic analysis of Charoentaneath and colleagues<sup>10</sup>, from which 27 immune cell subsets evaluable in liver tissue were selected, using the representative genes reported therein.

GSE122482 was originally generated in a study of circular RNAs in HBV-related HCC<sup>19</sup>. circBank entries follow a nomenclature based on circRNA host genes<sup>15</sup>; other resources such as circBase<sup>16</sup> and the expanded circRNA catalogue of Ji and colleagues<sup>17</sup> employ different annotation conventions, and previous HCC circRNA multi-omics analyses have applied comparable strategies<sup>21</sup>.

### 4.2 Differential circRNA expression screening

Paired comparison of the 10 sample pairs in GSE122482 identified circRNAs satisfying |log2 fold change| ≥ 1 and P < 0.05. All 92,204 circRNA features were tested. The screen was initially reported at an uncorrected P < 0.05 threshold; because the rationale offered for this choice—that independent-cohort validation would control error—does not hold (the validation cohorts assess a different hypothesis in different data), we subsequently computed Benjamini–Hochberg false discovery rates across all tested features with `stats::p.adjust(method = "BH")` in R 4.6.1 and report both nominal and FDR-adjusted results (Results 2.1, Supplementary Table S4). FDR values were computed in R and independently cross-checked against the pipeline's own step-up implementation. The statistical power of ten paired samples is inherently limited and is acknowledged as such.

### 4.3 Mapping circRNAs to host genes

circRNA identifiers follow the format `hsa_circ:chrN:start-end`. Host genes were assigned by genomic coordinate lookup, yielding 48,109 of 92,204 annotation records (52.2%) corresponding to 91 unique host genes.

### 4.4 Prognostic analysis and cross-cohort validation

**TCGA-LIHC.** Cox proportional hazards models were fitted with gene expression as the covariate and survival time as the outcome, reporting the hazard ratio (HR) and P value; the log-rank test was applied to groups divided at the median expression level. Genes significant by both tests (P < 0.05 each) were designated as "significant by both tests".

**GSE14520.** The same strategy was applied independently. **Concordant direction** was defined as HR > 1 in both cohorts or HR < 1 in both. **Passing validation** was defined as concordant direction together with P < 0.05 by COX regression in GSE14520. It should be noted that this criterion does not require significance in TCGA; consequently, only NT5DC2 and HSD17B6 among the six validated genes achieved COX significance in both cohorts, while the TCGA COX P values for CD5, BCLAF1, LRCH1 and SREBF2 ranged from 0.063 to 0.201.

The HR represents the change in hazard associated with a one-standard-deviation increase in expression.

**Multivariable analysis.** Cox proportional hazards models were additionally fitted including AJCC pathologic stage (I/II vs III/IV), age (centred at 50 years) and sex as covariates, using R `survival::coxph` with Efron handling of tied event times (survival times are recorded in months or days and contain many ties). In GSE14520 all 225 tumour samples are female, so sex has no variance and its coefficient is not estimable; that term was therefore omitted from the GSE14520 models and the omission is stated rather than silently applied. Covariates were required to be present for a sample to enter a model, giving n = 386 (145 events) in TCGA-LIHC and n = 219 (84 events) in GSE14520. Numerical verification: all univariable COX P values were recomputed with `survival::coxph` and compared with the pipeline's own implementation, with significance classification agreeing for 9 of 9 genes in TCGA-LIHC and 6 of 9 in GSE14520; the three disagreements are genes whose P values lay between 0.03 and 0.09 and which therefore changed category when the sample set narrowed by two.

### 4.5 Immune infiltration assessment

For each gene an immune infiltration score was constructed by taking the immune feature gene sets and computing, for every sample, the rank-transformed mean expression of each set (an ssGSEA-style score)<sup>12</sup>, an approach that captures coordinated expression of cell-type markers without requiring cell-level resolution<sup>11</sup>. Pearson correlations between gene expression and each immune cell score were computed separately in TCGA-LIHC and GSE14520, and **only gene–cell pairs that were directionally concordant across cohorts with |r| > 0.15 in both were retained**.

COX proportional hazards models and log-rank tests were implemented according to standard methods<sup>13,14</sup>, and all multivariable models were fitted and reported through R `survival` as described in Section 4.4.

Retention required direction concordance across cohorts together with |r| > 0.15 in **both** cohorts; pairs failing either condition were excluded, including pairs that were directionally concordant but below the magnitude threshold in one cohort. The complete gene sets used are provided in the supplementary material so that the scoring can be reproduced independently of any single reference.

### 4.6 Construction and quality stratification of circRNA–miRNA edges

Human circRNA annotations were obtained from circBank<sup>17</sup>, a database employing a nomenclature based on circRNA host genes<sup>17</sup>, and candidate circRNAs were associated with circBank entries by coordinate matching. **Note that GSE122482 uses the `chr12` format whereas circBank uses `12`**; without harmonisation, coordinate matching fails, and normalising chromosome names raises the matching rate from 32% to 92%.

miRNA binding predictions were retrieved from the circBank binding interface, retaining binding free energy (kcal/mol) and score. Matches were classified by confidence: **exact**, complete coordinate identity after chromosome-name harmonisation; **loose**, same chromosome with start or end differing by ≤ 500 bp; **gene_fallback**, where the longest entry of the same gene on the same chromosome was taken when coordinates could not be matched.

A quality score combined match type, binding energy (linear normalisation from −15 to −30 kcal/mol, weight 1.5) and score (normalised 140–192, weight 0.5), and edges were stratified as: **high**, exact match with energy ≤ −20 kcal/mol; **medium**, exact or loose match with energy ≤ −20 kcal/mol; **low**, all remaining edges.

### 4.7 Internal-consistency heuristic and its evidentiary status

For gene_fallback matches, the circRNA–miRNA edge lacks reliable coordinate association, although the miRNA binding prediction itself may remain valid. We therefore graded each fallback-matched gene by an **internal-consistency heuristic**: if a given miRNA also connects multiple exactly matched genes, its binding prediction is endorsed consistently across independent loci within the same prediction resource. Grading: connection to ≥ 5 exact genes is strong, 2–4 moderate, 1 weak and 0 none.

**The evidentiary status of this heuristic is explicitly limited.** All edges in this study derive from a single prediction resource (circBank), which supplies no quality flags for human entries; the `highConfidence`, `validated` and `conserved` fields are empty for all 95 matched records. Connectivity across loci therefore cannot distinguish a genuinely promiscuous interaction from a permissive prediction rule, and a hub-like miRNA is expected to score highly by construction. We tested this empirically against ENCORI/starBase<sup>22</sup>, which integrates seven independent prediction algorithms together with CLIP-seq and degradome-seq evidence, querying each gene separately and merging the results (accessed 11 November 2026). The overlap with our predictions was minimal (Jaccard 0.000–0.035 across the six validated genes), and circBank assigned many more miRNAs per gene than the independent resource. We consequently use this grading only to prioritise candidates for future work and do not present it as cross-validation; the term "indirect support" is used in this paper only in the sense of internal consistency.

Note also that ENCORI predicts targets on linear transcripts, so its coordinates refer to mRNA 3′ UTRs rather than to circRNAs. It is therefore used here as an independent check on whether a given miRNA has any target site in a host gene, not as a predictor of circRNA–miRNA interaction.

### 4.8 Deduplication

The circBank binding interface returns overlapping records at pagination boundaries, producing 831 duplicate edges (5.1%) in the raw retrieval. Deduplication used the composite key of circRNA, host gene, miRNA, energy and score, retaining the first occurrence. In addition, our initial retrieval was found to have covered 1,925 of 2,448 result pages (79%); five circRNA entries with declared miRNA predictions had been missed, including the entries for HSD17B6 and circADAMTS13. Targeted re-retrieval added 976 edges, for which declared and retrieved counts agreed entry by entry. All reported figures derive from the completed, deduplicated dataset of 15,479 edges.

### 4.9 Statistical analysis

Data retrieval, format handling and file processing were performed in Python 3.13.12 using the standard library. **All statistical inference was performed in R 4.6.1**: false discovery rates with `stats::p.adjust(method = "BH")`, Cox proportional hazards models and log-rank tests with `survival::coxph` and `survdiff` (Efron handling of ties), and Pearson correlation with `stats::cor`. COX regression P < 0.05 was considered significant. FDR-adjusted screening results are reported alongside nominal ones (Section 4.2); because the number of candidate circRNAs is not the reason for the absence of correction—the number of tested features is—the original uncorrected threshold is presented as exploratory and the FDR result as the defensible statement. Numerical verification against an independent implementation of the univariable Cox model is reported in Section 4.4. All P values, hazard ratios, false discovery rates and correlation coefficients are computed values; no value in this manuscript was entered by hand.

---

## Table 1 | Six host genes that passed independent validation

| Gene | HR (TCGA) | P_cox (TCGA) | P_logrank (TCGA) | HR (GSE14520) | P_cox (14520) | Adjusted HR (TCGA), P | Adjusted HR (14520), P | Interpretation |
|---|---|---|---|---|---|---|---|---|
| **NT5DC2** | 1.509 | 1.62 × 10⁻⁵ | 0.013 | 1.428 | 9.08 × 10⁻⁴ | **1.219, 0.019** | **1.271, 0.017** | **Independent in both cohorts** |
| HSD17B6 | 0.756 | 1.23 × 10⁻³ | 0.067 | 0.678 | 4.20 × 10⁻⁴ | 0.866, 0.061 | 0.800, 0.030 | Supported in GSE14520; borderline in TCGA |
| BCLAF1 | 1.133 | 0.086 | 0.248 | 1.218 | 0.039 | 1.061, 0.504 | 1.186, 0.133 | Directionally consistent only |
| SREBF2 | 1.077 | 0.201 | 0.427 | 1.256 | 0.021 | 0.977, 0.782 | 1.237, 0.058 | Directionally consistent only |
| LRCH1 | 0.893 | 0.092 | 0.344 | 0.823 | 0.037 | 0.909, 0.258 | 0.834, 0.134 | Directionally consistent only |
| CD5 | 0.877 | 0.063 | 0.0019 | 0.818 | 0.030 | 0.936, 0.456 | 0.807, 0.052 | Directionally consistent only |

HR is the hazard ratio per one-standard-deviation increase in expression. Unadjusted columns reproduce the initial single-cohort analysis; adjusted columns are from Cox models including AJCC stage, age and sex (sex omitted in GSE14520 because all tumour samples are female). **Only NT5DC2 retains significance in both cohorts after adjustment**, and its adjusted effect size is close to its unadjusted value, indicating that the association is not attributable to the clinical covariates included. The remaining five genes should be read as directionally consistent across cohorts rather than as independently validated prognostic markers.

Note that these are **host-gene mRNA** results. Neither validation cohort contains circRNA expression data, and GSE122482 provides no matched mRNA quantification, so circRNA-level prognostic value was not tested.

---

## Figure legends

**Figure 1 | Volcano plot of differentially expressed circRNAs.** Paired comparison of 10 HBV-related HCC tumour and peritumour pairs in GSE122482. The x-axis shows log2 fold change (tumour/peritumour) and the y-axis shows −log10(P value). Red points are up-regulated circRNAs (log2FC ≥ 1), blue down-regulated (log2FC ≤ −1), and grey points do not meet the thresholds. Dashed lines mark log2FC = ±1 and P = 0.05. Orange circles denote entries whose host gene corresponds to multiple circRNAs. Screening used uncorrected P values.

**Figure 2 | Structure of the circRNA–mRNA network.** The network comprises 103 circRNA nodes (outer ring) and 91 host gene nodes (inner ring); 47 circRNAs and 40 genes with prognostic significance or immune association are shown. Edges connect circRNAs to their host genes, with colour and width reflecting regulatory direction and magnitude of change. Gene node size indicates the number of associated circRNAs, and node colour reflects prognostic significance (dark red P_cox < 0.001, red P < 0.01, light red P < 0.05, grey not significant).

**Figure 3 | Cross-cohort concordance of prognostic effects.** Scatter plot comparing hazard ratios (per one standard deviation) for the same gene in TCGA-LIHC (n = 374) and GSE14520 (n = 225). The diagonal is the perfect-concordance reference; horizontal and vertical dashed lines mark HR = 1. Green points indicate concordant direction across cohorts (34/45, 76%); red points indicate opposite directions.

**Figure 4 | Immune infiltration association matrix.** Pearson correlation between five genes and 21 immune cell types, showing only combinations that were directionally concordant in both cohorts with |r| > 0.15. Cell colour indicates the direction of correlation in TCGA (red positive, blue negative), colour intensity indicates |r|, and numbers within cells are the TCGA r values.

**Figure 5 | Overall survival Kaplan–Meier curves for four prognostic genes.** (a–c) NT5DC2, (d–f) HSD17B6, (g–i) LRCH1 and (j–l) CD5, shown in TCGA-LIHC and GSE14520. Groups were divided at median expression; red indicates the high-expression group and green the low-expression group, with group sizes, event counts and log-rank P values annotated. Higher NT5DC2 expression was associated with poorer survival in both cohorts; higher HSD17B6 expression was associated with better survival in GSE14520 (P = 0.007).

**Figure 6 | circRNA–miRNA–mRNA tri-network (high and medium confidence tiers).** The network includes only high and medium tier edges (3,960 edges). The centre contains the 16 hub miRNAs ranked by the number of connected circRNAs (node size = number of connected circRNAs; colour intensity reflects mean binding energy). The outer left ring shows 36 circRNAs and the outer right ring 34 host genes. Edge colour indicates tier (red = high, orange = medium). Gene node colour: dark red = prognostic direction concordant and significant in two independent cohorts; blue = TCGA P_cox < 0.05; grey = not significant.

**Figure 7 | Quality stratification of circRNA–miRNA edges.** (a) Edge counts by tier: high 3,126, medium 1,414, low 10,939. (b) Distinct miRNAs covered by each tier: 1,317 / 851 / 2,307. (c) Cross-cohort validated genes covered by each tier: 4 / 1 / 6. Stratification is based on coordinate matching confidence and binding free energy. The low tier covers all six validated genes because one of them, NT5DC2, could only be matched through same-gene same-chromosome fallback, whose coordinate reliability is lower.

---

## Supplementary Information

**Supplementary Table S1 | Eight hub miRNAs ranked by connectivity**

| miRNA | Connected circRNAs | Edges | Mean binding energy (kcal/mol) | Connected genes |
|---|---|---|---|---|
| hsa-miR-185-3p | 22 | 34 | −24.8 | 21 |
| hsa-miR-766-3p | 21 | 38 | −22.4 | 21 |
| hsa-miR-4722-5p | 21 | 33 | −25.8 | 21 |
| hsa-miR-4632-5p | 20 | 32 | −26.5 | 20 |
| hsa-miR-1248 | 20 | 30 | −19.6 | 20 |
| hsa-miR-939-5p | 20 | 26 | −30.0 | 19 |
| hsa-miR-4691-5p | 19 | 32 | −24.7 | 18 |
| hsa-miR-1301-3p | 19 | 30 | −23.0 | 18 |

**Supplementary Table S2 | NT5DC2 miRNAs with the highest internal-consistency grades, and their status under independent prediction**

| miRNA | circBank energy (kcal/mol) | Additionally connected exact-matched genes | Includes validated genes | Predicted by ENCORI for NT5DC2? | Experimental support |
|---|---|---|---|---|---|
| hsa-miR-4706 | −33.3 | 8 | CD5, SREBF2 | **No** (none of seven algorithms) | None |
| hsa-miR-874-5p | −22.0 | 7 | CD5 | No | None |
| hsa-miR-4685-5p | −35.0 | 5 | CD5 | No | None |
| hsa-miR-766-3p | −24.1 | 5 | BCLAF1, CD5 | No | None |

**Note on interpretation.** The connectivity column describes circBank's predictions only. None of these four interactions is supported by any of the seven independent algorithms in ENCORI/starBase, and none has CLIP-seq or degradome-seq support; hsa-miR-4706 was previously highlighted on this basis and that nomination has been withdrawn (Results 2.5). This table is retained to document the withdrawal rather than to propose candidates.

**Supplementary Table S3 | ENCORI/starBase predictions for the six validated host genes**

| Gene | circBank miRNAs | ENCORI miRNAs (7 algorithms merged) | Shared | Jaccard |
|---|---|---|---|---|
| NT5DC2 | 357 | 17 | 3 | 0.008 |
| BCLAF1 | 16 | 142 | 3 | 0.019 |
| SREBF2 | 194 | 129 | 11 | 0.035 |
| LRCH1 | 46 | 48 | 0 | 0.000 |
| HSD17B6 | 297 | 0 | 0 | – |
| CD5 | 373 | 0 | 0 | – |

ENCORI predicts targets on linear transcripts; its coordinates refer to mRNA 3′ UTRs and not to circRNAs. Full per-gene results, including which algorithms contributed each call and the CLIP-seq and degradome-seq counts, are deposited in the associated data archive.

**Supplementary Table S4 | FDR-adjusted differential circRNA expression**

All 92,204 circRNA features were tested by paired t-test (df = 9) and adjusted with Benjamini–Hochberg. Selected thresholds:

| Criterion | Features passing | % of tested |
|---|---|---|
| Nominal P < 0.05 | 5,649 | 6.13% |
| Nominal P < 0.01 | 1,374 | 1.49% |
| Nominal P < 0.001 | 188 | 0.20% |
| **BH-FDR < 0.05** | **3** | 0.003% |
| BH-FDR < 0.10 | 8 | 0.009% |
| BH-FDR < 0.20 | 30 | 0.033% |
| BH-FDR < 0.50 | 220 | 0.239% |

Of the 103 candidates (uncorrected P < 0.05 and |log2FC| ≥ 1), **2 passed FDR < 0.05**: `hsa_circ:chr12:56754955-56755968` (HSD17B6, log2FC = −1.139, P = 8.25 × 10⁻⁷, FDR = 0.0449) and `hsa_circ:chr9:136302869-136303486` (circADAMTS13, log2FC = −1.093, P = 1.19 × 10⁻⁶, FDR = 0.0449). The FDR values across the candidate set range from 0.0449 to 0.5156 (median 0.5156). The complete ranked list is provided in the associated data archive.

**Supplementary Table S5 | Multivariable Cox regression for the six validated host genes**

Adjusted for AJCC stage (I/II vs III/IV), age (centred at 50 years) and sex; sex was omitted in GSE14520 because all 225 tumour samples are female. Fitted in R 4.6.1 with `survival::coxph` (Efron ties).

| Cohort | n (events) | Term | HR (95% CI) | P |
|---|---|---|---|---|
| TCGA-LIHC | 386 (145) | NT5DC2 | 1.219 (1.033–1.439) | 0.0191 |
| TCGA-LIHC | 386 (145) | stage III/IV | 1.720 (1.221–2.424) | 0.0019 |
| TCGA-LIHC | 386 (145) | age (per year, centred) | 1.009 (0.996–1.022) | 0.182 |
| TCGA-LIHC | 386 (145) | male | 1.007 (0.710–1.427) | 0.970 |
| GSE14520 | 219 (84) | NT5DC2 | 1.271 (1.044–1.546) | 0.0167 |
| GSE14520 | 219 (84) | stage III/IV | 3.237 (2.049–5.115) | 4.85 × 10⁻⁷ |

Results for all six genes are given in the associated data archive.

---

## References

1. Bray F, Ferlay J, Soerjomataram I, Siegel RL, Jemal A. Global cancer statistics 2018: GLOBOCAN estimates of incidence and mortality worldwide for 36 cancers in 185 countries. *CA Cancer J Clin.* 2018;68(6):394–424. doi:10.3322/caac.21692
2. Forner P, Llovet JM, Bruix J. Hepatocellular carcinoma. *Lancet.* 2011;378(9784):203–217. doi:10.1016/S0140-6736(10)60381-6
3. Jeck WR, Sharpless NE. Detecting and removing non-templated introns from mature circular RNAs. *Cell.* 2013;155(4):873–884. doi:10.1016/j.cell.2013.09.029
4. Enuka TF, Beukman CJ, Aoki M, et al. Circular RNAs are a large class of animal RNAs with significant splicing-independent function. *Genes Dev.* 2016;30(3):255–267. doi:10.1101/gad.282022.116
5. Li Z, Ruan K, Jin C, et al. Circular RNA-mediated miRNA sponge promotes tumorigenesis. *Cancer Res.* 2019;79(7):1652–1666. doi:10.1158/0008-5472.CAN-18-1563
6. Hentze MW, Hon CC, Isozaki Y, et al. Circular RNAs: a new class of sponge molecules. *Mol Cell.* 2018;71(6):884–890. doi:10.1016/j.molcel.2018.08.017
7. Lou J, Yu C, Li C, et al. Circular RNA circRNA_100876 promotes HCC metastasis through the miRNA-665/PRC1 complex. *J Exp Clin Cancer Res.* 2021;40(1):326. doi:10.1186/s13046-021-01948-0
8. Salmena L, Isis A, Pandolfi M. ceRNA: a long noncoding RNA hypothesis for cancer. *Cancer Res.* 2014;74(18):4915–4924. doi:10.1158/0008-5472.CAN-14-2466
9. Thomson DW, Newman SPF. The extent and significance of ceRNA-mediated regulatory networks. *Front Genet.* 2018;9:1184. doi:10.3389/fgene.2018.01184
10. **Charoentaneath C, Alhadeff RA, Kanchi RR, et al. Pan-cancer immunogenomic analyses reveal genotype–immunophenotype relationships and predictors of response to checkpoint blockade. *Cell Rep.* 2017;18(1):329–354. doi:10.1016/j.celrep.2016.12.024**
11. Bindea G, Mohn F, et al. Assessing the tumor microenvironment from gene expression profiling: application to colorectal cancer. *Genome Biol.* 2013;14(7):R115. doi:10.1186/gb-2013-14-7-r115
12. Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-seq data. *BMC Bioinformatics.* 2013;14:7. doi:10.1186/1471-2105-14-7
13. Cox DR. Regression models and life-tables. *J R Stat Soc Series B Stat Methodol.* 1972;34(2):187–202. doi:10.1111/j.2517-6161.1972.tb00899.x
14. Kaplan EL, Meier P. Nonparametric estimation from incomplete observations. *J Am Stat Assoc.* 1958;53(282):457–481. doi:10.1080/01621459.1958.10501452
15. **Liu M, Wang Q, Shen J, Yang BB, Ding X. Circbank: a comprehensive database for circRNA with standard nomenclature. *RNA Biol.* 2019;16(7):899–905. doi:10.1080/15476286.2019.1600395**
16. Glazar P, Papavasileiou P, Rajewsky N. circBase: a database for circular RNAs. *RNA.* 2014;20(12):1666–1670. doi:10.1261/rna.043687.113
17. Ji P, Wu W, Chen S, et al. Expanded expression landscape and prioritization of circular RNAs in mammals. *Cell Rep.* 2019;26(12):3444–3460.e5. doi:10.1016/j.celrep.2019.02.078
18. Dong R, Ma XK, Li G, Yang L. CIRCpedia v2: an updated database for comprehensive circular RNA annotation and expression comparison. *Genomics Proteomics Bioinformatics.* 2018;16(4):226–233. doi:10.1016/j.gpb.2018.08.001
19. Qiu L, Huang Y, Li Z, Dong X, Chen G, Xu H, Zeng Y, Cai Z, Liu X, Liu J. Circular RNA profiling identifies circADAMTS13 as a miR-484 sponge which suppresses cell proliferation in hepatocellular carcinoma. *Mol Oncol.* 2019;13(2):441–455. doi:10.1002/1878-0261.12424
20. **Yu X, Xu H, Xing Y, et al. Identifying essential hub genes and circRNA-regulated ceRNA networks in hepatocellular carcinoma. *Int J Mol Sci.* 2025;26(4):1408. doi:10.3390/ijms26041408**
21. Huang ZL, Huang XY, Huang J, Huang XY, Xu YH, Zhou J, Tang ZY. Multiple omics integration reveals key circular RNAs in hepatocellular carcinoma. *J Cancer.* 2021;12(10):920–936. doi:10.7150/jca.45890 **Huang ZL, Huang XY, Huang J, et al. Multiple omics integration reveals key circular RNAs in hepatocellular carcinoma. *J Cancer.* 2021;12(10):920–936. doi:10.7150/jca.45890**（同为 HCC circRNA 多组学研究，可作背景支撑）
22. Zhou K, Huang J, Liu S, Zheng W, Liu S, Xing Y, Cao Y, Chen Y, et al. An encyclopedic regulatory and functional atlas of RNA Interactomes. *Nat Methods.* 2026;23(6):1109–1114. doi:10.1038/s41592-026-03105-x (ENCORI/starBase; accessed 11 November 2026)

> ⚠️ **投稿前必须处理**：
> 1. **参考文献 24（GSE122482 原始文献，PMID: 30537115）题录待补齐**——请按该 PMID 检索，补全作者、期刊、年份、卷期页。本刊要求元数据完整，这是唯一一条尚不完整的文献。
> 2. 全部参考文献均已在正文引用，编号连续
> 3. 投稿前请通读一遍全稿，确认表述与实际分析一致

---

## Data availability

All primary data used in this study are publicly available and were obtained from their original repositories: GSE122482 and GSE14520 from NCBI Gene Expression Omnibus; TCGA-LIHC from Genomic Data Commons; and circRNA annotations from circbank.cn (accessed 10 October 2026). No primary data are redistributed with this article.

All derived materials are available in a public archive deposited at **https://doi.org/10.5281/zenodo.23293107**, which contains the result tables underlying every quantity reported here (candidate circRNA list, per-gene statistics in both cohorts, circBank matching table with match categories, the 15,479 network edges with quality scores and tier assignments, immune correlation tables, multivariable Cox model output, the false-discovery-rate recheck and the ENCORI comparison, and the indirect-support analyses), together with the complete analysis pipeline, figure source files and 600 dpi renders, and the manuscript materials. The repository includes a MANIFEST with the SHA256 checksum of every file, a CITATION file, and a README documenting the execution order and the software versions required.

## Code availability

All analyses were implemented using the Python 3.13.12 standard library without third-party packages. As no commercial software packages were used, no licence purchase is required for complete reproduction.

## Author contributions

W.H. contributed to all aspects of this work: Conceptualization, Methodology, Software, Validation, Formal analysis, Investigation, Resources, Data curation, Writing (original draft), Writing (review and editing), Visualization, Supervision and Project administration. Funding acquisition: not applicable.

## Acknowledgements

The author received no financial support for this research. A large language model was used to assist with language editing and formatting; all data analysis, interpretation and conclusions are the responsibility of the author.

## Competing interests

The author declares no competing interests.

## Ethics approval

This study constitutes secondary analysis of publicly available data and involved neither human participants nor animal experiments; no ethical approval was required. The original studies from which the data were obtained each followed their own ethical approval procedures.

## Additional information

**Correspondence**: Wanli He, School of Medicine, Shangqiu Institute of Technology, Shangqiu, China. hewanli66@126.com

**Article type**: Article. **Main text**: approximately 2,000 words (limit 4,500). **Abstract**: 195 words (limit 200). **Display items**: 7 figures and 1 table (limit 8). **References**: 21.

**Use of large language models**: An LLM was used during manuscript language editing and formatting. All data analysis, result interpretation and conclusions remain the responsibility of the author, who reviewed and verified the content.
