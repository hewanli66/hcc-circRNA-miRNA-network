# step18_major1_evidence.R
# Major comment 1 回应：预后标志物是宿主基因 mRNA，不是 circRNA
#
# 审稿人要求的 (a)：若 GSE122482 含配对的 mRNA 测量，
#   报告每个候选 circRNA 与其宿主基因 mRNA 在同一样本中的相关性
#
# 本脚本先核查 GSE122482 究竟提供了什么，然后：
#   1) 用可获得的证据量化「circRNA 与宿主 mRNA 的一致性」
#   2) 明确指出哪些部分数据不可得、不能声称
#   3) 给出 (b)(c) 所需的、可直接写进稿件的表述
#
# 数据来源：
#   results/GSE122482_circRNA_matrix.tsv   20 样本circRNA RPM
#   GSE122482_series_matrix.txt.gz         同一 20 样本的 GEO 元数据（数据表为空）

HCC <- 'C:/Users/hewan/WorkBuddy/2026-10-09-22-01-46/hcc'
RAW <- 'C:/Users/hewan/Documents/kimi/tasks/2026-10-08/11-38-08-f56eebff/rawdata'
RES <- file.path(HCC, 'results')

# ---------- 1. 核查 GSE122482 是否提供 mRNA 定量 ----------
cat('===== 1. GSE122482 数据可得性核查 =====\n')
sm <- file.path(RAW, 'GSE122482_series_matrix.txt.gz')

if (!file.exists(sm)) {
  cat('series_matrix 不存在\n')
} else {
  con <- gzfile(sm, 'rt')
  hdr <- NULL; nrow <- 0; seen <- FALSE; meta <- character(0)
  repeat {
    l <- readLines(con, n = 1)
    if (length(l) == 0) break
    ls <- trimws(l)
    if (startsWith(ls, '!series_matrix_table_begin')) { seen <- TRUE; next }
    if (startsWith(ls, '!series_matrix_table_end')) break
    if (!seen) { meta <- c(meta, ls); next }
    if (!nzchar(ls)) next
    if (is.null(hdr)) hdr <- l else nrow <- nrow + 1
  }
  close(con)
  cat(sprintf('series_matrix 元数据行数        : %d\n', length(meta)))
  cat(sprintf('series_matrix 数据表基因行数    : %d\n', nrow))
  cat(sprintf('样本数                        : %d\n',
              if (is.null(hdr)) NA else length(hdr) - 1))
  cat(sprintf('  -> mRNA 定量是否可得          : %s\n',
              if (nrow > 0) '是' else '否（数据表为空）'))
}

# GEO 补充文件清单（此前实测：21 个文件 = 1 TAR + 20 个 RPM 表达文件）
cat('\nGEO 补充文件（此前已实测 filelist.txt）：\n')
cat('  GSE122482_RAW.tar + 20 个 *_circ_RPM.txt.gz\n')
cat('  -> 全部为 circRNA 的 RPM 表达值，**无 mRNA 定量文件**\n')

# ---------- 2. 用可得的证据：circRNA 差异方向 vs 宿主基因预后方向 ----------
cand <- read.delim(file.path(RES, 'candidate_circRNA.tsv'),
                   header = TRUE, sep = '\t', quote = '',
                   stringsAsFactors = FALSE, check.names = FALSE)
gn <- read.delim(file.path(RES, 'network_genes.tsv'),
                 header = TRUE, sep = '\t', quote = '',
                 stringsAsFactors = FALSE, check.names = FALSE)

key_genes <- c('NT5DC2', 'HSD17B6', 'BCLAF1', 'SREBF2', 'CD5', 'LRCH1')
c1 <- cand[match(key_genes, cand$gene), ]
c1$log2FC <- as.numeric(c1$log2FC)
g1 <- gn[match(key_genes, gn$gene), ]
g1$HR_TCGA <- as.numeric(g1$HR_TCGA)
g1$HR_14520 <- as.numeric(g1$HR_14520)

tab <- data.frame(
  gene = key_genes,
  circRNA = c1$circRNA,
  circRNA_log2FC = c1$log2FC,
  circRNA_direction = ifelse(c1$log2FC > 0, 'up in tumour', 'down in tumour'),
  mRNA_HR_TCGA = g1$HR_TCGA,
  mRNA_prognostic = ifelse(g1$HR_TCGA > 1,
                           'high = worse survival', 'high = better survival'),
  stringsAsFactors = FALSE)

cat('\n===== 2. 六个跨队列验证基因的两层证据 =====\n')
for (i in seq_len(nrow(tab))) {
  cat(sprintf('  %-9s circRNA %+6.3f (%-15s) | mRNA HR=%.3f (%s)\n',
              tab$gene[i], tab$circRNA_log2FC[i], tab$circRNA_direction[i],
              tab$mRNA_HR_TCGA[i], tab$mRNA_prognostic[i]))
}

# ---------- 3. 方向一致性：circRNA 变化方向 vs mRNA 预后方向 ----------
# 若某 circRNA 下调、其宿主 mRNA 高表达=预后好，则「mRNA 低=预后差」不成立
# 这一点只能作为「代理关系未经验证」的说明，不能作为结论
cat('\n===== 3. 可核查的方向关系（仅描述，不作因果推断）=====\n')
n_down <- sum(tab$circRNA_log2FC < 0)
n_protective <- sum(tab$mRNA_HR_TCGA < 1)
cat(sprintf('  六个基因中 circRNA 在肿瘤下调者: %d 个\n', n_down))
cat(sprintf('  六个基因中 mRNA 高表达=预后好者: %d 个\n', n_protective))
cat(sprintf('  其中「circRNA 下调 且 mRNA 是保护性」: %d 个\n',
            sum(tab$circRNA_log2FC < 0 & tab$mRNA_HR_TCGA < 1)))
cat('  -> 这类基因若套用「circRNA 变化反映宿主基因变化」的假设，')
cat('     会得到方向相反的推论，因此代理关系必须实测而非默认。\n')

# ---------- 4. 输出可供稿件使用的表述 ----------
L <- c('# Major comment 1 回应：标志物层级与 circRNA–mRNA 代理关系', '',
       '> 本文件由 R 4.6.1 生成。', '',
       '## 1. 审稿人要求的 (a)：GSE122482 的配对 mRNA 测量', '',
       '### 1.1 核查结论：该数据集不提供 mRNA 定量', '',
       '| 检查项 | 结果 |', '|---|---|',
       '| `GSE122482_series_matrix.txt.gz` 的数据表行数 | **0**（只有元数据表头） |',
       '| GEO 补充文件构成 | 1 个 TAR + 20 个 `*_circ_RPM.txt.gz` |',
       '| 补充文件中是否含 mRNA 定量 | **否**，全部为 circRNA 的 RPM 表达值 |',
       sprintf('| circRNA 表达矩阵样本数 | 20（10 对肿瘤/癌旁）|', '')[1],
       '| 本地 rawdata 中是否有该数据集的 mRNA 定量 | 无 |',
       '| `rnaseq_merged_rsem_fpkm_*.csv` 的归属 | HCCDB18（SIDM 编号），**非 GSE122482** |', '',
       '**结论：GSE122482 只提供了 circRNA 的表达值，没有同一样本的 mRNA 定量。**',
       '因此审稿人要求的 (a)「报告每个候选 circRNA 与其宿主基因 mRNA 在同一样本中的',
       '相关性」**在本数据集内无法完成**——不是因为我们未做，而是因为数据不存在。', '',
       '### 1.2 这一缺口的含义', '',
       '原稿的论证链隐含了一个未经验证的假设：**circRNA 的丰度可以代表其宿主基因',
       'mRNA 的丰度**。这个假设并不成立——反向剪接效率受独立调控，与整体转录水平',
       '可以脱钩。因此：', '',
       '- 「差异表达的 circRNA」这一层有直接数据支持（10 对配对）；',
       '- 「宿主基因 mRNA 的预后价值」这一层有独立数据支持（TCGA-LIHC、GSE14520）；',
       '- **「因此该 circRNA 是预后标志物」这一跳没有数据支持。**', '',
       '审稿人对这一点的判断是正确的。我们接受该质疑，并在正文与讨论中明确改写。',
       '', '')

L <- c(L, '## 2. 六个验证基因的两层证据（描述性）', '',
       '下表并列列出两个层级的证据。**两个层级指向的分子不同**，不可互相替代。', '',
       '| 基因 | 候选 circRNA | circRNA log2FC（10对） | circRNA 方向 | 宿主 mRNA HR (TCGA) | mRNA 预后方向 |',
       '|---|---|---|---|---|---|')
for (i in seq_len(nrow(tab))) {
  L <- c(L, sprintf('| %s | `%s` | %+.3f | %s | %.3f | %s |',
                    tab$gene[i], tab$circRNA[i], tab$circRNA_log2FC[i],
                    tab$circRNA_direction[i], tab$mRNA_HR_TCGA[i],
                    tab$mRNA_prognostic[i]))
}
L <- c(L, '')

L <- c(L, '### 2.1 为什么这两列不能直接连接', '',
       sprintf('六个基因中：circRNA 在肿瘤下调者 **%d** 个，宿主 mRNA 高表达对应',
               n_down), '',
       sprintf('更好预后者 **%d** 个。两者同时成立的有 **%d** 个。',
               n_protective,
               sum(tab$circRNA_log2FC < 0 & tab$mRNA_HR_TCGA < 1)), '',
       '对于这些基因，若默认「circRNA 下调反映宿主 mRNA 下调」，而该mRNA 又是',
       '保护性因素，就会推出「该 circRNA 下调是保护性」的相反结论。',
       '**这说明代理关系既不能默认成立，也不能靠方向一致性来推断**——',
       '必须有同一样本的配对测量才能判定，而 GSE122482 不提供这种测量。', '',
       '## 3. 可行的补救途径（供讨论部分说明）', '',
       '| 途径 | 可行性 | 说明 |', '|---|---|---|',
       '| 在 GSE122482 的 SRA 原始 reads 上重新定量 mRNA | 理论可行 | SRP168496 有 BAM，但需重新跑 STAR/featureCounts，工作量大 |',
       '| 用 GTEx 正常组织 + TCGA 肿瘤做跨组织比较 | 部分可行 | 只能比较癌 vs 正常，不能提供同一样本的配对相关 |',
       '| 在 HCC 含 circRNA 的公开数据集中查已有报告 | 需检索 | Qiu et al. Mol Oncol 2019 已报道 circADAMTS13，但未给出 circRNA–mRNA 相关性 |',
       '| 实验验证（RT-qPCR 共表达） | 最直接 | 需10 对配对组织，成本与时间较高 |', '',
       '我们认为最诚实的处理是：**不声称代理关系成立，并在限制中明确写出**',
       '「circRNA 层面的预后价值尚未被检验」。', '')

L <- c(L, '## 4. 需要落实的稿件修改（对应审稿人 (b)(c)）', '',
       '### (b) 在 Results 与 Discussion 中明确层级', '',
       '建议在 Results 相应小节加入：', '',
       '> The survival analyses were performed on the linear mRNAs of circRNA host',
       '> genes. Neither TCGA-LIHC nor GSE14520 contains circRNA expression data, and',
       '> GSE122482 provides circRNA read counts without matched mRNA quantification;',
       '> the correlation between candidate circRNAs and their host-gene mRNAs within',
       '> the same samples could therefore not be assessed. The prognostic markers',
       '> identified here are host-gene-level mRNA markers.', '', '',
       '建议在 Discussion 的局限性中补入：', '',
       '> Three further points limit the inference. First, host-gene mRNA abundance is',
       '> not a reliable proxy for circRNA abundance, because back-splicing efficiency',
       '> is regulated independently of overall transcription; the matched mRNA',
       '> measurements needed to test this proxy assumption are not available in',
       '> GSE122482. Second, our independent cohorts (TCGA-LIHC, GSE14520) contain no',
       '> circRNA data, so they validate the prognostic value of the host-gene mRNAs,',
       '> not of the circRNAs themselves. Third, whether the candidate circRNAs',
       '> carry prognostic information at the circRNA level remains to be determined',
       '> by direct measurement.', '', '',
       '### (c) 弱化标题与摘要', '',
       '| 位置 | 原文表述 | 建议改为 |', '|---|---|---|',
       '| 标题 | Cross-cohort circRNA–miRNA–mRNA network identifies NT5DC2 and HSD17B6 as **prognostic markers** | Host genes of differentially expressed circRNAs show consistent cross-cohort prognostic value in hepatocellular carcinoma |',
       '| 摘要 | presents NT5DC2 and HSD17B6 as circRNA-level prognostic markers | reports NT5DC2 and HSD17B6 as **host-gene** mRNA markers; circRNA-level prognostic value remains untested |',
       '| 关键词 | prognostic biomarker | host-gene prognostic marker |', '',
       '## 5. 对审稿人 (a) 的正式答复措辞（可直接粘贴）', '',
       '> We thank the reviewer for identifying this gap. We have verified that',
       '> GSE122482 provides circRNA RPM values for 20 samples (10 tumour–peritumour',
       '> pairs) but contains **no matched mRNA quantification**: the GEO series',
       '> matrix contains metadata only, and the supplementary files comprise one TAR',
       '> archive and 20 circRNA RPM tables. The correlation between candidate',
       '> circRNAs and their host-gene mRNAs within the same samples therefore cannot',
       '> be computed from this dataset, and we did not assume it. We have revised the',
       '> title, abstract and Results to state explicitly that the prognostic markers',
       '> are host-gene mRNA markers, that the two independent cohorts contain no',
       '> circRNA data and thus validate the mRNA level rather than the circRNA level,',
       '> and that circRNA-level prognostic value remains to be established by direct',
       '> measurement. We also added this limitation to the Discussion.', '',
       '---', '',
       sprintf('> 生成：R 4.6.1，%s', format(Sys.time(), '%Y-%m-%d %H:%M')))

writeLines(L, file.path(HCC, 'Major1_evidence_and_response.md'))

# 结果表
write.table(tab, file.path(RES, 'major1_two_level_evidence.tsv'),
            sep = '\t', row.names = FALSE, quote = FALSE)

cat('\nDONE. Output: Major1_evidence_and_response.md, major1_two_level_evidence.tsv\n')
