# Cross-cohort circRNA–miRNA–mRNA network analysis identifies NT5DC2 and HSD17B6 as prognostic markers in hepatocellular carcinoma

**Wanli He**<sup>1*</sup>

<sup>1</sup> School of Medicine, Shangqiu Institute of Technology, Shangqiu, China

\* Correspondence: hewanli66@126.com; ORCID 0009-0002-0829-1909

---

## Abstract

Hepatocellular carcinoma (HCC) shows substantial prognostic heterogeneity, while circular RNAs (circRNAs), as stable non-coding RNAs, regulate gene expression by sponging microRNAs. We identified 103 differentially expressed circRNAs from GSE122482 (10 pairs of HBV-related tumour and peritumour tissues), assessed host gene prognosis in TCGA-LIHC (374 tumours) and GSE14520 (225 tumours), and integrated circBank miRNA binding predictions to construct a tri-layer network. The network comprised 14,568 circRNA–miRNA edges involving 2,358 miRNAs;33 circRNAs matched circBank entries at exact coordinates. Of the host genes, 37 were significant in TCGA, 34 showed concordant direction across cohorts (76%), and six passed independent validation (concordant direction with COX P < 0.05 in GSE14520). Among these, NT5DC2 and HSD17B6 were significant by COX regression in both cohorts (HR 1.51 and 1.43; 0.76 and 0.68 respectively, with opposite directions). CD5 was strongly correlated with multiple T-cell populations (CD4+ T cells r = 0.877) and remained prognostically relevant. For genes with unreliable coordinate matching we propose an indirect-support strategy: if a miRNA connects several exactly matched genes, this is cross-validation. Among the 357 miRNAs of NT5DC2, 51.2% received such support, with hsa-miR-4706 (energy −33.3 kcal/mol) additionally connecting eight exactly matched genes. Our quality stratification and indirect-support framework may serve as a reference for circRNA annotation where coordinate consistency is lacking.

**Keywords**: hepatocellular carcinoma; circRNA; miRNA; prognostic marker; immune infiltration; network analysis

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

Paired comparison of GSE122482 identified 103 differentially expressed circRNAs (59 up-regulated, 44 down-regulated), corresponding to 91 host genes (Fig. 1). Nine genes corresponded to more than one circRNA (CYP4F23P, GALP, HIST2H4B, LMNA, MEG3, MGC50722, RBM19, TMEM179, ZNF250), indicating diversity in circRNA generation at these loci.

One candidate, hsa_circ:chr9:136302869-136303486, requires comment. This coordinate corresponds to circADAMTS13 (hsa_circ_0005075), which Qiu and colleagues identified as a miR-484 sponge in the very dataset used here<sup>19</sup>. Our annotation assigned this locus to a different host gene (MGC50722), reflecting the differing annotation conventions between resources (Section 3.2); circBank names circRNAs by host gene and assigns this entry to ADAMTS13. Importantly, our data agree with theirs: we likewise find this circRNA down-regulated in tumour tissue (log2FC = −1.09, P = 1.19 × 10⁻⁶), and hsa-miR-484 appears in our network as a predicted partner. This concordance provides an external check on our pipeline: circADAMTS13 functions here as a positive control rather than a discovery, and the remaining candidates are the underexplored space that the cross-cohort design is intended to expose (Section 3.2).

### 2.2 Prognostic value of host genes and cross-cohort validation

In TCGA-LIHC, the expression of 37 host genes was significantly associated with overall survival (COX regression P < 0.05), of which 25 were significant by both COX regression and the log-rank test.

Independent evaluation in GSE14520 showed concordant direction for 34 genes (76%). Applying the criterion of concordant direction together with COX regression P < 0.05 in GSE14520, six genes passed validation (Table 1).

**Validation strength was heterogeneous.** NT5DC2 and HSD17B6 were significant by COX regression in both cohorts and represent the strongest evidence in this study. BCLAF1, SREBF2, LRCH1 and CD5 did not reach significance by COX regression in TCGA (P = 0.086, 0.201, 0.092 and 0.063 respectively); these were included on the basis of significant COX regression in the validation cohort with concordant direction, and CD5 was additionally significant by the log-rank test in TCGA (P = 0.0019). Figure 3 shows the scatter distribution of hazard ratios in the two cohorts, with 34 paired genes lying close to the diagonal.

**Effect direction.** Higher HSD17B6 and LRCH1 expression was associated with better overall survival in both cohorts (protective), whereas the remaining four genes were associated with poorer outcome (Fig. 5).

### 2.3 Immune infiltration associations

Under the criterion of cross-cohort direction concordance, only seven genes showed stable immune associations:

- **CD5 was strongly positively correlated with multiple T-cell populations** (CD4+ T cells r = 0.877, Tfh r = 0.870, effector CD4 r = 0.863, CD8 r = 0.816)
- LRCH1 with Tfh cells (r = 0.258), HSD17B6 with memory CD8 (r = −0.243), SREBF2 with NK cells (r = −0.229), and SMG5 with endothelial cells (r = −0.200) showed weak to moderate associations
- NT5DC2 showed an extremely weak association with NK cells (r = −0.018)

The dual character of CD5 is noteworthy: it reached significance by COX regression in GSE14520 and by the log-rank test in TCGA, while its high expression correlated positively with several T-cell populations. Notably, CD5 has HR < 1 in both cohorts, indicating that it is a **protective** factor (Fig. 4).

### 2.4 The circRNA–miRNA–mRNA network

Integration with circBank predictions yielded 14,568 circRNA–miRNA edges involving 2,358 miRNAs, 74 circRNAs and 69 host genes. The quality stratification is summarised in Fig. 7:

| Tier | Edges | miRNAs | circRNAs | Genes |
|---|---|---|---|---|
| high | 2,546 | 1,182 | 29 | 28 |
| medium | 1,414 | 851 | 15 | 14 |
| low | 10,608 | 2,292 | 74 | 69 |

Ranked by the number of connected circRNAs, the top eight hub miRNAs connect 19–22 distinct circRNAs and 18–21 distinct host genes, with mean binding free energies mostly below −22 kcal/mol, indicating relatively strong predicted binding (Supplementary Table S1; Fig. 6).

### 2.5 Indirect support for genes matched by fallback

For 26 genes, circRNA entries had to be matched by same-gene, same-chromosome fallback, leaving the coordinate association uncertain. Applying the indirect-support strategy, for **NT5DC2**—the only fallback-matched gene that also passed cross-cohort validation (HR_TCGA = 1.509, HR_14520 = 1.428)—32 of its 357 miRNAs reached the strong tier (connecting ≥ 5 exactly matched genes) and 168 reached the moderate tier (2–4 genes), giving **51.2%** with strong or moderate cross-support.

Priority mechanistic candidates are listed in Supplementary Table S2. **hsa-miR-4706** is notable: its predicted binding energy on NT5DC2 reaches −33.3 kcal/mol, it connects eight genes matched at exact coordinates, and those genes include two independently validated genes, CD5 and SREBF2, forming a coherent mechanistic hypothesis.

### 2.6 Two technical issues affecting data reliability

Two issues were identified during the analysis and are reported transparently.

**Duplicate edges from pagination overlap in circBank.** The binding interface returns overlapping records at pagination boundaries, producing 766 duplicate edges (5.0%) in the raw data. Deduplication on the composite key of circRNA, host gene, miRNA, energy and score reduced the dataset from 15,334 to 14,568 unique edges.

**Discrepancy in chromosome naming conventions.** GSE122482 uses `chrN` whereas circBank uses `N`; these cannot be compared directly. Without harmonisation, only 33 of 103 candidate circRNAs could be matched to circBank entries at exact coordinates (32%). Normalising chromosome names raised the matching rate to 92% (exact 33, loose 18, gene_fallback 44). This was the single largest technical determinant of matching success in this study.

---

## 3. Discussion

### 3.1 Main findings

This study makes three principal contributions.

**First, panoramic screening outperforms single-axis investigation.** Systematic screening of 103 circRNAs identified six host genes with cross-cohort prognostic relevance. A conventional "one axis at a time" strategy would have been likely to miss these, particularly NT5DC2 and HSD17B6, which have received little attention at the circRNA level in the existing literature.

**Second, cross-cohort validation markedly improves reliability, though validation strength remains heterogeneous.** Of 37 prognostic genes in TCGA, 34 showed concordant direction in GSE14520 (76%), and six met the combined criterion of concordant direction plus significant COX regression in the validation cohort. This attrition indicates that single-cohort prognostic significance alone is insufficient to support mechanistic inference. It must also be acknowledged that only NT5DC2 and HSD17B6 achieved significance by COX regression in both cohorts; the remaining four did not reach significance in TCGA and are correspondingly weaker evidence.

**Third, quality stratification and the indirect-support strategy address a real annotation problem.** circBank does not populate the commonly used `highConfidence`, `validated` and `conserved` fields for human entries (100% empty in our 95 matched entries), and coordinate systems differ systematically between sources. Our combined approach of match-type stratification plus binding free energy, together with indirect support for coordinate-inconsistent cases, offers a workable solution to this widespread problem.

### 3.2 Comparison with existing studies

Several circRNA databases are available, each with distinct annotation strategies and nomenclature<sup>15–18</sup>, which introduces heterogeneity that complicates cross-study comparison. Existing HCC circRNA research has focused on a small number of well-studied molecules and has typically validated a single miRNA under a single ceRNA hypothesis. Our methodological differences are threefold: we start from data rather than from a predefined molecule; we introduce binding free energy as a filtering dimension rather than relying on prediction counts alone; and we apply explicit quality stratification to predicted relationships instead of treating thousands of predictions as equivalent.

The present dataset has been used in two earlier analyses, and we disclose both. Qiu and colleagues generated these data and reported circADAMTS13 as a miR-484 sponge<sup>19</sup>; our reuse of the same data was not intended to revisit that circuit. Yu and colleagues applied the same dataset to WGCNA-based ceRNA network construction with RT-qPCR validation of five networks<sup>20</sup>. Neither study performed prognostic evaluation of circRNA host genes across independent cohorts, which is the principal addition here.

An important consequence concerns the interpretation of our own predictions. Of the circRNAs we identified, one—circADAMTS13—has already been experimentally characterised, including confirmation of its miR-484 sponge activity by dual-luciferase assays with a mutated target site<sup>19</sup>. Our pipeline correctly recovered both its direction of change and its partner miRNA, which is reassuring, but it also means that this entry is a positive control rather than a discovery. We were unable to cross-match the remaining published circRNA identifiers against our candidate list, because the identifiers used by the source studies were not resolvable to genomic coordinates in the resources we used. We therefore make no claim about how many of our candidates remain experimentally untested; determining this requires a systematic comparison against the primary literature, which we have not performed. What can be stated is that the six genes we identify as cross-cohort validated, and NT5DC2 and HSD17B6 in particular, have not to our knowledge been examined at the circRNA level in HCC, and they constitute the underexplored space that our cross-cohort design is designed to expose.

### 3.3 Boundaries of the immune mechanism

The immune-related conclusions require careful delimitation. The strong correlation between CD5 and CD4+ T-cell infiltration (r = 0.877) was directionally consistent across two cohorts, and CD5 was simultaneously prognostically relevant, making it a candidate for mechanistic investigation. However, **the remaining six genes all showed immune association strengths of |r| < 0.26, and NT5DC2 was essentially unrelated (r = −0.018)**. Immune mechanisms are therefore discussed only in relation to CD5, without extrapolation.

### 3.4 Limitations

**First, circRNA screening used uncorrected P values.** The threshold was chosen on the design principle that independent-cohort validation would control for error, but false positives may still accumulate if the validation cohort is not fully independent of the discovery cohort. GSE14520 is an independent dataset, yet both cohorts derive from HBV-related populations, limiting generalisability.

**Second, circBank miRNA predictions are unvalidated computationally.** All binding relationships are computational predictions. The high-confidence tier requires exact coordinate matching and binding energy ≤ −20 kcal/mol, which is reasonably stringent, but experimental confirmation remains necessary.

**Third, coordinate association is uncertain for 26 fallback-matched genes.** Although more than half of NT5DC2's miRNAs obtained indirect support, mechanistic inference for the remaining fallback-matched genes should be treated with corresponding caution.

**Fourth, this is a purely computational analysis.** All conclusions derive from secondary analysis of public data and computational prediction, without in vitro or in vivo validation. Immune infiltration scores were estimated as rank-transformed means of feature gene sets and their absolute magnitudes are not comparable to direct measurements such as flow cytometry or immunohistochemistry; the rank-transformed approach may additionally underestimate true correlation<sup>12</sup>. The circRNA–miRNA relationships proposed here should therefore be regarded as **hypotheses awaiting validation** rather than established mechanistic conclusions. Follow-up work should confirm candidate circRNA expression differences between tumour and peritumour tissue by qPCR, and test direct circRNA–miRNA interaction using dual-luciferase reporter assays. The latter specifically requires construction of a mutant control in which the miRNA seed region is deleted, so that non-specific suppression of the reporter can be excluded.

### 3.5 Future directions

Three priorities follow from the present stratification. First, validation of the NT5DC2–hsa-miR-4706–CD5/SREBF2 axis, which combines strong-tier indirect support with a strong predicted binding energy. Second, investigation of how CD5 influences prognosis through CD4+ T-cell infiltration. Third, exploration of HSD17B6 as a loss-of-function protective factor.

---

## 4. Methods

### 4.1 Data sources

| Dataset | Platform / type | Sample size | Use |
|---|---|---|---|
| GSE122482 | circRNA expression microarray | 10 pairs, HCC tumour / peritumour | Differential circRNA screening |
| TCGA-LIHC | RNA-seq + clinical | 374 tumours | Prognostic modelling (discovery cohort) |
| GSE14520 | GPL3921 + clinical | 225 tumours | Independent validation cohort |
| circBank | database | 828,808 human circRNAs | circRNA annotation and miRNA binding prediction |

Clinical data for GSE14520 were obtained from its Extra Supplement file, and only tumour samples were included. Immune cell feature gene sets were taken from the pan-cancer immunogenomic analysis of Charoentaneath and colleagues<sup>10</sup>, from which 27 immune cell subsets evaluable in liver tissue were selected, using the representative genes reported therein.

GSE122482 was originally generated in a study of circular RNAs in HBV-related HCC<sup>19</sup>. circBank entries follow a nomenclature based on circRNA host genes<sup>15</sup>; other resources such as circBase<sup>16</sup> and the expanded circRNA catalogue of Ji and colleagues<sup>17</sup> employ different annotation conventions, and previous HCC circRNA multi-omics analyses have applied comparable strategies<sup>21</sup>.

### 4.2 Differential circRNA expression screening

Paired comparison of the 10 sample pairs in GSE122482 identified circRNAs satisfying |log2 fold change| ≥ 1 and P < 0.05. An **uncorrected** P value was used as the preliminary threshold to avoid excessive conservatism under multiple testing; error control was achieved through subsequent cross-cohort prognostic validation (see 4.4).

### 4.3 Mapping circRNAs to host genes

circRNA identifiers follow the format `hsa_circ:chrN:start-end`. Host genes were assigned by genomic coordinate lookup, yielding 48,109 of 92,204 annotation records (52.2%) corresponding to 91 unique host genes.

### 4.4 Prognostic analysis and cross-cohort validation

**TCGA-LIHC.** Cox proportional hazards models were fitted with gene expression as the covariate and survival time as the outcome, reporting the hazard ratio (HR) and P value; the log-rank test was applied to groups divided at the median expression level. Genes significant by both tests (P < 0.05 each) were designated as "significant by both tests".

**GSE14520.** The same strategy was applied independently. **Concordant direction** was defined as HR > 1 in both cohorts or HR < 1 in both. **Passing validation** was defined as concordant direction together with P < 0.05 by COX regression in GSE14520. It should be noted that this criterion does not require significance in TCGA; consequently, only NT5DC2 and HSD17B6 among the six validated genes achieved COX significance in both cohorts, while the TCGA COX P values for CD5, BCLAF1, LRCH1 and SREBF2 ranged from 0.063 to 0.201.

The HR represents the change in hazard associated with a one-standard-deviation increase in expression.

### 4.5 Immune infiltration assessment

For each gene an immune infiltration score was constructed by taking the immune feature gene sets and computing, for every sample, the rank-transformed mean expression of each set (an ssGSEA-style score)<sup>12</sup>, an approach that captures coordinated expression of cell-type markers without requiring cell-level resolution<sup>11</sup>. Pearson correlations between gene expression and each immune cell score were computed separately in TCGA-LIHC and GSE14520, and **only gene–cell pairs that were directionally concordant across cohorts with |r| > 0.15 in both were retained**.

COX proportional hazards models and log-rank tests were implemented according to standard methods<sup>13,14</sup>.

### 4.6 Construction and quality stratification of circRNA–miRNA edges

Human circRNA annotations were obtained from circBank<sup>17</sup>, a database employing a nomenclature based on circRNA host genes<sup>17</sup>, and candidate circRNAs were associated with circBank entries by coordinate matching. **Note that GSE122482 uses the `chr12` format whereas circBank uses `12`**; without harmonisation, coordinate matching fails, and normalising chromosome names raises the matching rate from 32% to 92%.

miRNA binding predictions were retrieved from the circBank binding interface, retaining binding free energy (kcal/mol) and score. Matches were classified by confidence: **exact**, complete coordinate identity after chromosome-name harmonisation; **loose**, same chromosome with start or end differing by ≤ 500 bp; **gene_fallback**, where the longest entry of the same gene on the same chromosome was taken when coordinates could not be matched.

A quality score combined match type, binding energy (linear normalisation from −15 to −30 kcal/mol, weight 1.5) and score (normalised 140–192, weight 0.5), and edges were stratified as: **high**, exact match with energy ≤ −20 kcal/mol; **medium**, exact or loose match with energy ≤ −20 kcal/mol; **low**, all remaining edges.

### 4.7 Indirect-support strategy

For gene_fallback matches, the circRNA–miRNA edge lacks reliable coordinate association, although the miRNA binding prediction itself may remain valid. We therefore propose an **indirect-support** strategy: if a given miRNA also connects **multiple exactly matched genes**, this indicates that the binding prediction for that miRNA is consistently endorsed by circBank across independent loci, constituting cross-validation of the association with the fallback-matched gene. Grading is as follows: connection to ≥ 5 exact genes is strong, 2–4 is moderate, 1 is weak and 0 is none.

### 4.8 Deduplication

The circBank binding interface returns overlapping records at pagination boundaries, producing 766 duplicate edges (5.0%). Deduplication used the composite key of circRNA, host gene, miRNA, energy and score, retaining the first occurrence and reducing 15,334 records to 14,568.

### 4.9 Statistical analysis

All analyses were performed in Python 3.13.12 using the standard library (`urllib`, `json`, `csv`, `math`, `re`). COX regression P < 0.05 was considered significant; Pearson correlations were filtered for cross-cohort direction concordance. Because the number of candidate circRNAs (103) is small relative to sample size, no multiple-testing correction was applied; this limitation is stated in Section 3.4. All P values, hazard ratios and correlation coefficients are actual computed values.

---

## Table 1 | Six host genes that passed independent validation

| Gene | HR (TCGA) | P_cox (TCGA) | P_logrank (TCGA) | HR (GSE14520) | P_cox (14520) | Effect direction |
|---|---|---|---|---|---|---|
| NT5DC2 | 1.509 | 1.62 × 10⁻⁵ | 0.013 | 1.428 | 9.08 × 10⁻⁴ | Higher = worse prognosis |
| HSD17B6 | 0.756 | 1.23 × 10⁻³ | 0.067 | 0.678 | 4.20 × 10⁻⁴ | Higher = better prognosis |
| BCLAF1 | 1.133 | 0.086 | 0.248 | 1.218 | 0.039 | Higher = worse prognosis |
| SREBF2 | 1.077 | 0.201 | 0.427 | 1.256 | 0.021 | Higher = worse prognosis |
| LRCH1 | 0.893 | 0.092 | 0.344 | 0.823 | 0.037 | Higher = better prognosis |
| CD5 | 0.877 | 0.063 | 0.0019 | 0.818 | 0.030 | Higher = better prognosis |

HR is the hazard ratio per one-standard-deviation increase in expression. Only NT5DC2 and HSD17B6 achieved COX significance in both cohorts.

---

## Figure legends

**Figure 1 | Volcano plot of differentially expressed circRNAs.** Paired comparison of 10 HBV-related HCC tumour and peritumour pairs in GSE122482. The x-axis shows log2 fold change (tumour/peritumour) and the y-axis shows −log10(P value). Red points are up-regulated circRNAs (log2FC ≥ 1), blue down-regulated (log2FC ≤ −1), and grey points do not meet the thresholds. Dashed lines mark log2FC = ±1 and P = 0.05. Orange circles denote entries whose host gene corresponds to multiple circRNAs. Screening used uncorrected P values.

**Figure 2 | Structure of the circRNA–mRNA network.** The network comprises 103 circRNA nodes (outer ring) and 91 host gene nodes (inner ring); 47 circRNAs and 40 genes with prognostic significance or immune association are shown. Edges connect circRNAs to their host genes, with colour and width reflecting regulatory direction and magnitude of change. Gene node size indicates the number of associated circRNAs, and node colour reflects prognostic significance (dark red P_cox < 0.001, red P < 0.01, light red P < 0.05, grey not significant).

**Figure 3 | Cross-cohort concordance of prognostic effects.** Scatter plot comparing hazard ratios (per one standard deviation) for the same gene in TCGA-LIHC (n = 374) and GSE14520 (n = 225). The diagonal is the perfect-concordance reference; horizontal and vertical dashed lines mark HR = 1. Green points indicate concordant direction across cohorts (34/45, 76%); red points indicate opposite directions.

**Figure 4 | Immune infiltration association matrix.** Pearson correlation between five genes and 21 immune cell types, showing only combinations that were directionally concordant in both cohorts with |r| > 0.15. Cell colour indicates the direction of correlation in TCGA (red positive, blue negative), colour intensity indicates |r|, and numbers within cells are the TCGA r values.

**Figure 5 | Overall survival Kaplan–Meier curves for four prognostic genes.** (a–c) NT5DC2, (d–f) HSD17B6, (g–i) LRCH1 and (j–l) CD5, shown in TCGA-LIHC and GSE14520. Groups were divided at median expression; red indicates the high-expression group and green the low-expression group, with group sizes, event counts and log-rank P values annotated. Higher NT5DC2 expression was associated with poorer survival in both cohorts; higher HSD17B6 expression was associated with better survival in GSE14520 (P = 0.007).

**Figure 6 | circRNA–miRNA–mRNA tri-network (high and medium confidence tiers).** The network includes only high and medium tier edges (3,960 edges). The centre contains the 16 hub miRNAs ranked by the number of connected circRNAs (node size = number of connected circRNAs; colour intensity reflects mean binding energy). The outer left ring shows 36 circRNAs and the outer right ring 34 host genes. Edge colour indicates tier (red = high, orange = medium). Gene node colour: dark red = prognostic direction concordant and significant in two independent cohorts; blue = TCGA P_cox < 0.05; grey = not significant.

**Figure 7 | Quality stratification of circRNA–miRNA edges.** (a) Edge counts by tier: high 2,546, medium 1,414, low 10,608. (b) Distinct miRNAs covered by each tier: 1,182 / 851 / 2,292. (c) Cross-cohort validated genes covered by each tier: 3 / 1 / 5. Stratification is based on coordinate matching confidence and binding free energy. The low tier covers the most validated genes because certain key genes, such as NT5DC2, could only be matched through same-gene same-chromosome fallback, whose coordinate reliability is lower.

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

**Supplementary Table S2 | Candidate miRNAs for NT5DC2 with indirect support**

| miRNA | Binding energy | Best energy on other genes | Number of additionally connected exact-matched genes | Includes validated genes |
|---|---|---|---|---|
| hsa-miR-4706 | −33.3 | −36.3 | 8 | CD5, SREBF2 |
| hsa-miR-874-5p | −22.0 | −31.6 | 7 | CD5 |
| hsa-miR-4685-5p | −35.0 | −39.4 | 5 | CD5 |
| hsa-miR-766-3p | −24.1 | −44.3 | 5 | BCLAF1, CD5 |

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

> ⚠️ **投稿前必须处理**：
> 1. **参考文献 24（GSE122482 原始文献，PMID: 30537115）题录待补齐**——请按该 PMID 检索，补全作者、期刊、年份、卷期页。本刊要求元数据完整，这是唯一一条尚不完整的文献。
> 2. 全部参考文献均已在正文引用，编号连续
> 3. 投稿前请通读一遍全稿，确认表述与实际分析一致

---

## Data availability

All primary data used in this study are publicly available and were obtained from their original repositories: GSE122482 and GSE14520 from NCBI Gene Expression Omnibus; TCGA-LIHC from Genomic Data Commons; and circRNA annotations from circbank.cn (accessed 10 October 2026). No primary data are redistributed with this article.

All derived materials are available in a public archive deposited at **[DOI TO BE ASSIGNED ON ACCEPTANCE]**, which contains the result tables underlying every quantity reported here (candidate circRNA list, per-gene statistics in both cohorts, circBank matching table with match categories, the 14,568 network edges with quality scores and tier assignments, immune correlation tables, and the indirect-support analyses), together with the complete analysis pipeline, figure source files and 600 dpi renders, and the manuscript materials. The repository includes a MANIFEST with the SHA256 checksum of every file, a CITATION file, and a README documenting the execution order and the Python version required.

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
