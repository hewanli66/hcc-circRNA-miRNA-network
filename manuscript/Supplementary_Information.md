# Supplementary Information

**Cross-cohort circRNA–miRNA–mRNA network analysis identifies NT5DC2 and HSD17B6 as prognostic markers in hepatocellular carcinoma**

Wanli He
School of Medicine, Shangqiu Institute of Technology, Shangqiu, China
hewanli66@126.com｜ORCID 0009-0002-0829-1909

---

## Contents

| Item | Description |
|---|---|
| Supplementary Table S1 | Eight hub miRNAs ranked by connectivity |
| Supplementary Table S2 | Candidate miRNAs for NT5DC2 with indirect support |
| Supplementary Note S1 | Matching confidence categories and their assignment |
| Supplementary Note S2 | Quality scoring formula |
| Supplementary Note S3 | Duplicate edge handling |
| Supplementary Note S4 | Complete definition of indirect support tiers |
| Supplementary Note S5 | Composition of the immune cell feature gene sets |
| Supplementary Note S6 | Full numbers for the six validated genes |
| Supplementary Note S7 | Known limitation of the fallback matching strategy |
| Supplementary Table S3 | Files deposited in the associated repository |

---

## Supplementary Table S1 | Eight hub miRNAs ranked by connectivity

Hub miRNAs were ranked by the number of distinct circRNAs they connect (ranked first), with edges and mean binding energy reported for reference.

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

**Note.** All eight hub miRNAs connect 19–22 distinct circRNAs and 18–21 distinct host genes, with mean binding energies mostly below −22 kcal/mol. A previous version of this table ranked hub miRNAs by edge count rather than by the number of connected circRNAs; the two criteria give different orderings, and the present table uses the latter consistently with Figure 6.

---

## Supplementary Table S2 | Candidate miRNAs for NT5DC2 with indirect support

NT5DC2 was matched to circBank entries only by the gene_fallback route, so its circRNA–miRNA edges were not independently anchored by coordinate identity. The table lists miRNAs that received strong or moderate indirect support (Supplementary Note S4).

| miRNA | Binding energy | Best energy on other genes | Number of additionally connected exact-matched genes | Includes validated genes |
|---|---|---|---|---|
| hsa-miR-4706 | −33.3 | −36.3 | 8 | CD5, SREBF2 |
| hsa-miR-874-5p | −22.0 | −31.6 | 7 | CD5 |
| hsa-miR-4685-5p | −35.0 | −39.4 | 5 | CD5 |
| hsa-miR-766-3p | −24.1 | −44.3 | 5 | BCLAF1, CD5 |

**Note.** hsa-miR-4706 is the strongest candidate: it combines a strong predicted binding energy on NT5DC2, eight additionally connected exactly matched genes, and inclusion of two independently validated genes (CD5, SREBF2) among those supporting genes. The complete list of 357 NT5DC2 miRNAs with support tier assignments is deposited in the associated repository.

---

## Supplementary Note S1 | Matching confidence categories

Candidate circRNAs from GSE122482 were associated with circBank entries by coordinate comparison after harmonising chromosome naming. Three categories were used.

| Category | Definition | Count | Share |
|---|---|---|---|
| exact | Complete coordinate identity after harmonisation | 33 | 32.0% |
| loose | Same chromosome, start or end differing by ≤ 500 bp | 18 | 17.5% |
| gene_fallback | Coordinates could not be matched; the longest entry of the same gene on the same chromosome was taken | 44 | 42.7% |
| **Total matched** | | **95** | **92.2%** |

Eight of the 103 candidates could not be matched to any circBank entry and were excluded from network construction.

**Critical note on chromosome nomenclature.** GSE122482 uses the `chrN:start-end` convention whereas circBank uses `N:start-end`. Without harmonisation, only the 33 exact matches would have been found (32% rather than 92%). This was the single largest technical determinant of matching success in this study.

---

## Supplementary Note S2 | Quality scoring formula

Each edge received a quality score combining three normalised components:

```
quality = w_match + 1.5 × w_energy + 0.5 × w_score
```

where

- **w_match**: exact = 3.0, loose = 2.0, gene_fallback = 1.0
- **w_energy**: linear normalisation of binding energy from −15 kcal/mol (0) to −30 kcal/mol (1), clipped to [0, 1]
- **w_score**: linear normalisation of the circBank score from 140 (0) to 192 (1), clipped to [0, 1]

The score ranges from 1.0 to 5.0. Tier assignment used the match category and binding energy rather than the composite score:

| Tier | Rule |
|---|---|
| high | exact match and energy ≤ −20 kcal/mol |
| medium | exact or loose match and energy ≤ −20 kcal/mol |
| low | all remaining edges (including all gene_fallback edges) |

---

## Supplementary Note S3 | Duplicate edge handling

The circBank binding interface returns overlapping records at pagination boundaries. Across the 15,334 raw records this produced 766 duplicate edges, representing 5.0% of the dataset. Duplicates were removed using the composite key of circRNA, host gene, miRNA, energy and score, retaining the first occurrence and reducing the dataset to 14,568 unique edges.

A single primary key was insufficient, because the same miRNA can legitimately bind the same circRNA at different energies, and those records are distinct rather than duplicated.

---

## Supplementary Note S4 | Definition of indirect support tiers

For edges assigned by gene_fallback, the coordinate association is uncertain, but the underlying miRNA binding prediction may still be valid. Indirect support was defined as follows.

For a given miRNA m and a fallback-matched gene g:

```
support(m, g) = |{ g' : m connects to g', match_type(g') = exact, g' ≠ g }|
```

| Tier | Rule | Count (NT5DC2) |
|---|---|---|
| strong | support ≥ 5 | 28 |
| moderate | 2 ≤ support ≤ 4 | 151 |
| weak | support = 1 | 91 |
| none | support = 0 | 87 |

For NT5DC2, strong plus moderate accounts for 50.1% of its 357 miRNAs. Counts here are per unique miRNA; a small number of miRNAs appear more than once in the edge-level file because they bind more than one circRNA of the same gene. The logic is that a miRNA independently connected to several exactly matched genes has received consistent endorsement from circBank across separate loci, which corroborates the prediction for the fallback-matched gene.

---

## Supplementary Note S5 | Immune cell feature gene sets

Immune cell feature gene sets were taken from the pan-cancer immunogenomic analysis of Charoentaneath and colleagues (Cell Rep 2017;18:329–354). Twenty-seven immune cell subsets evaluable in liver tissue were selected, using the representative genes reported therein. For each gene, an immune infiltration score was computed as the rank-transformed mean expression of each gene set (an ssGSEA-style score), separately in TCGA-LIHC and GSE14520.

Only gene–cell pairs that were directionally concordant across cohorts with |r| > 0.15 in both cohorts were retained. Applying this filter left seven genes with stable immune associations:

| Gene | Cell type | r (TCGA) | Interpretation |
|---|---|---|---|
| CD5 | CD4+ T cell | +0.877 | strong |
| CD5 | Effector CD4 | +0.863 | strong |
| CD5 | CD8+ T cell | +0.816 | strong |
| LRCH1 | Tfh cell | +0.258 | weak–moderate |
| HSD17B6 | Memory CD8 | −0.243 | weak–moderate |
| SREBF2 | NK cell | −0.229 | weak |
| SMG5 | Endothelial | −0.200 | weak |
| NT5DC2 | NK cell | −0.018 | none |
| NUP205 | NK cell | −0.006 | none |

**Interpretation limit.** Except for CD5, all associations are weak (|r| < 0.26), and for NT5DC2 the association is effectively absent. Immune mechanisms are therefore discussed only in relation to CD5. Because the scores are rank-transformed gene-set means, their absolute magnitudes are not comparable to direct measurements such as flow cytometry or immunohistochemistry, and the rank transformation may underestimate true correlation.

---

## Supplementary Note S6 | Full numbers for the six validated genes

| Gene | HR (TCGA) | P_cox (TCGA) | P_logrank (TCGA) | HR (GSE14520) | P_logrank (14520) | Effect direction |
|---|---|---|---|---|---|---|
| NT5DC2 | 1.509 | 1.62 × 10⁻⁵ | 0.013 | 1.428 | 0.0053 | Higher = worse |
| HSD17B6 | 0.756 | 1.23 × 10⁻³ | 0.067 | 0.678 | 0.0066 | Higher = better |
| BCLAF1 | 1.133 | 0.086 | 0.248 | 1.218 | 0.175 | Higher = worse |
| SREBF2 | 1.077 | 0.201 | 0.427 | 1.256 | 0.189 | Higher = worse |
| LRCH1 | 0.893 | 0.092 | 0.344 | 0.823 | 0.187 | Higher = better |
| CD5 | 0.877 | 0.063 | 0.0019 | 0.818 | 0.066 | Higher = better |

Hazard ratios are per one-standard-deviation increase in gene expression. **Only NT5DC2 and HSD17B6 reached COX significance in both cohorts**; the remaining four were significant by COX regression only in the validation cohort, with TCGA COX P values ranging from 0.063 to 0.201, and their evidence is correspondingly weaker. CD5 is the only one of the six significant by the log-rank test in TCGA.

---

## Supplementary Note S7 | Known limitation of the fallback matching strategy

Twenty-six genes were matched by the gene_fallback route, which assigns a circRNA to the longest entry of the same gene on the same chromosome when exact coordinates fail. Two limitations follow.

First, the assigned entry may be a different circRNA from the one originally interrogated. The convergent case is instructive: the candidate hsa_circ:chr9:136302869-136303486 was annotated to host gene MGC50722 in the source data but named hsa_ADAMTS13_0000800 by circBank, which names circRNAs by host gene. Qiu and colleagues identified this locus as circADAMTS13 and experimentally confirmed its miR-484 sponge activity (Mol Oncol 2019;13:441–455). Our pipeline recovered both its direction of change and its miR-484 partner, which we regard as a positive control; but the entry cannot be presented as a discovery.

Second, for the remaining 25 fallback-matched genes the coordinate association remains unverified. The indirect-support strategy of Supplementary Note S4 mitigates but does not eliminate this problem, and mechanistic inference for these genes should be treated with corresponding caution.

We were unable to cross-match the circRNA identifiers used in the two prior studies of this dataset against our candidate list, because those identifiers could not be resolved to genomic coordinates in the resources used here. We therefore make no claim about how many of our candidates have previously been experimentally tested.

---

## Supplementary Table S3 | Files deposited in the associated repository

All files are deposited at **https://doi.org/10.5281/zenodo.23293107** (concept DOI, resolving to the most recent version; this package corresponds to v1.1.1), as cited in the Data availability statement of the main manuscript. No public data are redistributed; those are obtained from their original repositories using the accessions given in Methods 4.1.

**Results**

| File | Description |
|---|---|
| `candidate_circRNA.tsv` | 103 differentially expressed circRNAs with host gene, log2FC, P value and direction |
| `circbank_id_map.tsv` | Mapping of each candidate to its circBank entry, with match category |
| `network_edges_qc.tsv` | 14,568 circRNA–miRNA edges with quality score and tier assignment |
| `network_genes.tsv` | Host gene statistics in both cohorts |
| `trinet_nodes.tsv` | Tri-network node table (miRNA and circRNA layers) |
| `trinet_summary.tsv` | Tier summary counts |
| `immune_correlation_2queues.tsv` | Gene–immune cell correlations retaining cross-cohort concordant pairs |
| `gse14520_validation.tsv` | Per-gene validation statistics in GSE14520 |
| `NT5DC2_mirna_indirect_support.tsv` | Indirect support tier for each of the 357 NT5DC2 miRNAs |
| `fallback_gene_mirna_support.tsv` | Edge-level indirect support for all fallback-matched genes |
| `fallback_genes_indirect_support.tsv` | Per-gene summary of indirect support |

**Scripts**

| Script | Purpose |
|---|---|
| `step1_map_filter.py` | Differential circRNA screening and host gene mapping |
| `step2_tcga_survival.py` | COX and log-rank analysis in TCGA-LIHC |
| `step3_validate14520.py` | Independent validation in GSE14520 |
| `step4_immune.py`, `step4b_immune_tcga.py` | Immune infiltration scoring and correlation |
| `step6_network.py` | circRNA–miRNA–mRNA network assembly |
| `step7_miRNA_binding_v2.py` | circBank retrieval, coordinate harmonisation and edge construction |
| `step8_network_qc.py` | Deduplication, quality scoring and tier assignment |
| `step9_nt5dc2_support.py`, `step9b_fallback_support.py` | Indirect support analysis |
| `make_figs3.py` | Figure generation |

**Environment.** All analyses were performed in Python 3.13.12 using the standard library only. Scripts therefore require no third-party packages, except `make_figs3.py` for figure rendering.

---

*All quantities reported in this document and in the main manuscript are actual computed values derived from the deposited result files.*