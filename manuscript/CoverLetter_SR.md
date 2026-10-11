# Cover Letter — Scientific Reports

**To**: Editorial Board, *Scientific Reports*

**From**: Wanli He
School of Medicine, Shangqiu Institute of Technology, Shangqiu, China
hewanli66@126.com｜ORCID 0009-0002-0829-1909

**Date**: 11 October 2026

---

Dear Editors,

We submit for your consideration the manuscript **"Cross-cohort circRNA–miRNA–mRNA network analysis identifies NT5DC2 and HSD17B6 as prognostic markers in hepatocellular carcinoma"** for publication as an Article in *Scientific Reports*.

## Why this work suits Scientific Reports

Your journal assesses manuscripts on **scientific validity and technical soundness rather than perceived importance or impact**, and explicitly welcomes computational and bioinformatics analyses. This study is a computational analysis of public data, reported with that degree of methodological transparency: all code was written using the Python standard library with no third-party packages, and every quantity reported in the manuscript is an actual computed value. We have taken care to state our limitations plainly rather than to frame computational predictions as established mechanism (Section 3.4).

## What the study shows

We screened 103 differentially expressed circRNAs in paired HBV-related HCC tissue (GSE122482), evaluated the prognostic value of their host genes in TCGA-LIHC (374 tumours) and in the independent GSE14520 cohort (225 tumours), and integrated circBank miRNA binding predictions to build a circRNA–miRNA–mRNA network of 14,568 edges.

Three points may interest your readership:

1. **Cross-cohort validation materially changes the picture.** Of 37 host genes significant in TCGA, 34 showed concordant direction in GSE14520 (76%), but only six passed the combined criterion of concordant direction plus significant COX regression in the validation cohort. We report this explicitly rather than presenting the discovery-cohort result alone.

2. **Two genes are robust in both cohorts.** NT5DC2 and HSD17B6 were significant by COX regression in both TCGA and GSE14520 (HR 1.51 and 1.43; and 0.76 and 0.68, with opposite directions). We deliberately distinguish these two from the remaining four validated genes, whose TCGA COX P values range from 0.063 to 0.201 and whose evidence is correspondingly weaker.

3. **A practical annotation issue is addressed.** circBank leaves `highConfidence`, `validated` and `conserved` empty for human entries (100% empty in our matched set), and coordinate conventions differ systematically between data sources — GSE122482 uses `chrN` whereas circBank uses `N`. Harmonising chromosome names raised our matching rate from 32% to 92%. For genes matched only by fallback we propose an indirect-support criterion, which certifies a predicted miRNA when that miRNA is independently connected to several exactly matched genes. This certifies 51.2% of the miRNAs attached to NT5DC2.

## Disclosure of prior use of the same dataset

We wish to disclose this proactively. The circRNA expression dataset (GSE122482) has been used in two earlier studies: the study that generated it, which characterised circADAMTS13 as a miR-484 sponge, and a 2025 analysis that applied WGCNA-based network construction with RT-qPCR validation. Neither performed prognostic evaluation of circRNA host genes across independent cohorts, which is the principal addition here. We also note that circADAMTS13 appears among our candidates: our pipeline recovered both its direction of change and its miR-484 partner, which we regard as a positive control on the pipeline, and we have not presented it as a discovery. Neither study's circRNA identifiers could be resolved to genomic coordinates in the resources we used, so we make no claim about how many of our remaining candidates have been experimentally tested; that would require a systematic literature comparison we have not performed.

## Statements

The manuscript contains no competing interests. The author received no financial support for this work. All primary data are publicly available from their original repositories and are not redistributed. All derived materials — result tables, the complete analysis pipeline, figure files and manuscript materials — are archived at Zenodo (DOI assigned on acceptance); the archive includes a SHA256 manifest and a README documenting the execution order and Python version.

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
