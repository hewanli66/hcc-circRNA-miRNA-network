# Cover Letter — Scientific Reports

**To**: Editorial Board, *Scientific Reports*

**From**: Wanli He
School of Medicine, Shangqiu Institute of Technology, Shangqiu, China
hewanli66@126.com｜ORCID 0009-0002-0829-1909

**Date**: 11 October 2026

---

Dear Editors,

We submit for your consideration the manuscript **"Host genes of differentially expressed circRNAs show consistent cross-cohort prognostic value in hepatocellular carcinoma"** for publication as an Article in *Scientific Reports*.

## Why this work suits Scientific Reports

Your journal assesses manuscripts on **scientific validity and technical soundness rather than perceived importance or impact**, and explicitly welcomes computational and bioinformatics analyses. This study is a computational analysis of public data, reported with that degree of methodological transparency: data retrieval and file handling were performed in Python, and **all statistical inference was performed in R 4.6.1** (`survival::coxph` with Efron ties for the survival models, `stats::p.adjust(method = "BH")` for false discovery rates). Every quantity reported in the manuscript is a computed value, and the survival statistics were independently recomputed with a second implementation as a numerical check. We have taken care to state our limitations plainly rather than to frame computational predictions as established mechanism (Section 3.4).

## What the study shows

We screened 103 differentially expressed circRNAs in paired HBV-related HCC tissue (GSE122482), evaluated the prognostic value of their host genes in TCGA-LIHC (374 tumours) and in the independent GSE14520 cohort (225 tumours), and integrated circBank miRNA binding predictions to build a circRNA–miRNA–mRNA network of 15,479 edges.

Four points may interest your readership:

1. **The claims are graded by what the data actually support.** The survival analyses were performed on linear host-gene mRNAs, not on circRNAs: neither validation cohort contains circRNA expression data, and GSE122482 provides circRNA read counts without matched mRNA quantification. We state this explicitly in the title, abstract and Results, and describe the markers accordingly as host-gene-level mRNA markers, with circRNA-level prognostic value identified as untested.

2. **Multivariable adjustment leaves a single robust candidate.** Of 37 host genes significant in TCGA, 34 showed concordant direction in GSE14520 (76%) and six passed the combined validation criterion. After adjustment for AJCC stage, age and sex, **only NT5DC2 remained significant in both cohorts** (adjusted HR 1.22, P = 0.019; 1.27, P = 0.017) with an essentially unchanged effect size, indicating that its association is not attributable to the clinical covariates. HSD17B6 reached significance in GSE14520 but was borderline in TCGA (P = 0.061); the remaining four lost significance. We report this graded outcome rather than presenting all six as validated markers.

3. **We report what multiplicity correction does to the screen.** Across 92,204 circRNA features, 5,649 pass an uncorrected P < 0.05 but only 3 reach BH-FDR < 0.05, and only 2 of our 103 candidates do so. The original rationale—that cross-cohort validation would control the false-positive rate—does not hold, because the validation cohorts test a different hypothesis in different data; we say so and treat the screen as exploratory.

4. **We tested our own annotation heuristic and report the negative result.** All circRNA–miRNA edges derive from a single resource (circBank), which provides no quality flags for human entries. Comparing with ENCORI/starBase (seven independent algorithms plus CLIP-seq and degradome-seq evidence) gave a Jaccard overlap of only 0.000–0.035, and none of those algorithms supports the candidate interaction we had highlighted in an earlier version. We accordingly withdraw that nomination, reframe connectivity-based grading as internal consistency rather than cross-validation, and state the evidentiary limits in Methods and Discussion.

## Disclosure of prior use of the same dataset

We wish to disclose this proactively. The circRNA expression dataset (GSE122482) has been used in two earlier studies: the study that generated it, which characterised circADAMTS13 as a miR-484 sponge, and a 2025 analysis that applied WGCNA-based network construction with RT-qPCR validation. Neither performed prognostic evaluation of circRNA host genes across independent cohorts, which is the principal addition here. We also note that circADAMTS13 appears among our candidates: our pipeline recovered both its direction of change and its miR-484 partner, which we regard as a positive control on the pipeline, and we have not presented it as a discovery. Neither study's circRNA identifiers could be resolved to genomic coordinates in the resources we used, so we make no claim about how many of our remaining candidates have been experimentally tested; that would require a systematic literature comparison we have not performed.

## Statements

The manuscript contains no competing interests. The author received no financial support for this work. All primary data are publicly available from their original repositories and are not redistributed. All derived materials — result tables, the complete analysis pipeline, figure files and manuscript materials — are archived at Zenodo (https://doi.org/10.5281/zenodo.23293108); the archive includes a SHA256 manifest and a README documenting the execution order and the software versions used.

We have not previously submitted this manuscript to any journal, and it is not under consideration elsewhere. All authors approve the submission.

We hope you find the manuscript suitable for *Scientific Reports* and look forward to your response.

Yours sincerely,

**Wanli He**
School of Medicine, Shangqiu Institute of Technology, Shangqiu, China
hewanli66@126.com｜ORCID 0009-0002-0829-1909

---

## Suggested reviewers

*（以下为建议名单，投稿时请在系统内填写机构邮箱；此处仅列研究方向供作者参考）*

| 领域 | 方向 |
|---|---|
| circRNA / ceRNA 网络分析 | 熟悉 circRNA 数据库注释与网络构建 |
| HCC 预后建模 | 生存分析、多队列验证 |
| 免疫浸润计算 | ssGSEA / 肿瘤免疫微环境 |
| 公共数据二次分析 | 统计方法与可重复性 |

## 备选送审人排除名单

（投稿时如需排除特定审稿人，此处填写；本稿无需排除）

---

## 投稿前自查

- [ ] 全文核对所有 P 值、HR、r 值与源文件一致
- [ ] 补充参考文献 19 之外的必读文献（若审稿意见要求）
- [ ] 图件转 TIFF/EPS，线图 ≥ 300 dpi（当前为 SVG+PNG，需按期刊要求转换）
- [ ] 补充材料表格 S1、S2 单独成文件
- [ ] 在投稿系统填写推荐审稿人及其机构邮箱
- [ ] 确认数据可用性声明中的 Zenodo/figshare DOI（若已上传）
