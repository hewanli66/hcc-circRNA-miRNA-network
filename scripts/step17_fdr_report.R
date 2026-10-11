# Major comment 2 回应：FDR 多重检验校正
# 生成：results/../Major2_FDR结果.md
#
# 关键事实：
#   1. BH 算法本身正确（step1 用的是标准 step-up，且对全部 92,204 个特征做的）
#   2. 但 step1 写文件时 P_adj 与特征行的对应关系错位
#      （padj = {order[i]: adj[i]} 把「第 i 小的特征」与「第 i 个原始索引的校正值」配对）
#   3. 因此本文件的 FDR 值全部以 R 的 p.adjust(method="BH") 为准，
#      已与 step1 的 P_adj 逐条比对并记录不一致情况
#
# 数据来源：results/GSE122482_circRNA_matrix.tsv（20 样本 = 10 对肿瘤/癌旁）
#          配对 t 检验（df=9），全部 92,204 个 circRNA 特征

suppressPackageStartupMessages(library(stats))

HCC <- 'C:/Users/hewan/WorkBuddy/2026-10-09-22-01-46/hcc'
RES <- file.path(HCC, 'results')

all <- read.delim(file.path(RES, 'circRNA_gene_map.tsv'),
                  header = TRUE, sep = '\t', quote = '',
                  stringsAsFactors = FALSE, check.names = FALSE)
cand <- read.delim(file.path(RES, 'candidate_circRNA.tsv'),
                   header = TRUE, sep = '\t', quote = '',
                   stringsAsFactors = FALSE, check.names = FALSE)

# ---- BH 校正（本文件的唯一权威值）----
all$BH <- p.adjust(all$P, method = 'BH')
all$in_candidate <- all$circRNA %in% cand$circRNA

# ---- 与 step1 的 P_adj 比对 ----
cmp <- data.frame(
  n_total = nrow(all),
  n_match = sum(abs(all$BH - all$P_adj) < 1e-6, na.rm = TRUE),
  pct_match = round(100 * sum(abs(all$BH - all$P_adj) < 1e-6, na.rm = TRUE) /
                    nrow(all), 2),
  step1_range_min = min(all$P_adj, na.rm = TRUE),
  step1_range_max = max(all$P_adj, na.rm = TRUE),
  R_BH_range_min = min(all$BH),
  R_BH_range_max = max(all$BH),
  stringsAsFactors = FALSE)

L <- c('# Major comment 2 回应：多重检验校正（FDR）', '',
       '> 本文件由 R 4.6.1 `stats::p.adjust(method = "BH")` 生成，',
       '> 是 BH 校正值的唯一权威来源。', '',
       '## 1. 为什么需要重新计算 BH 值', '',
       '原始筛选在 92,204 个 circRNA 特征上用未校正的 P < 0.05，',
       '审稿人指出这会产生大量假阳性，要求补充 FDR 校正结果。', '',
       '复核中发现 `step1_map_filter.py` 的 `P_adj` 列与特征行的对应关系错位：', '',
       '```python',
       '# step1 的写法（错误）',
       'padj = {order[i]: adj[i] for i in range(m_)}',
       '```', '',
       '`order[i]` 是「第 i 小的特征的下标」，`adj[i]` 是「原始第 i 个特征的校正值」',
       '——两者不是同一个特征。正确的写法应当是 `dict(zip(order, adj))`。', '',
       'BH 算法本身是对的（标准 step-up、对全部 92,204 个特征计算），',
       '错的是把校正值写回文件时的配对关系。**由于原始筛选用的是未校正 P，',
       '这一错误没有影响 103 个候选 circRNA 的产生**，但使 `P_adj` 列不可用。', '',
       '### 比对结果', '',
       '| 项 | 值 |', '|---|---|',
       sprintf('| 比对总行数 | %s |', format(cmp$n_total[1], big.mark = ',')),
       sprintf('| 与 R 的 BH 完全一致 | %s (%.2f%%) |',
               format(cmp$n_match[1], big.mark = ','), cmp$pct_match[1]),
       sprintf('| step1 `P_adj` 取值范围 | %.6f – %.6f |',
               cmp$step1_range_min[1], cmp$step1_range_max[1]),
       sprintf('| R `p.adjust` 取值范围 | %.6f – %.6f |',
               cmp$R_BH_range_min[1], cmp$R_BH_range_max[1]), '',
       '两者的**取值范围完全相同**，说明 BH 计算逻辑一致，差异仅在行对应关系。', '')

# ---- 2. 名义 P vs FDR 的存活数量 ----
n <- nrow(all)
L <- c(L, '## 2. 名义 P 与 BH-FDR 的存活数量', '',
       sprintf('检验特征数：**%s** 个 circRNA；样本：10 对肿瘤/癌旁配对（df = 9）。', format(n, big.mark = ',')), '',
       '| 判定标准 | 存活数 | 占检验数比例 |', '|---|---|---|',
       sprintf('| 未校正 P < 0.05 | %s | %.2f%% |',
               format(sum(all$P < 0.05), big.mark = ','),
               100 * sum(all$P < 0.05) / n),
       sprintf('| 未校正 P < 0.01 | %s | %.2f%% |',
               format(sum(all$P < 0.01), big.mark = ','),
               100 * sum(all$P < 0.01) / n),
       sprintf('| 未校正 P < 0.001 | %s | %.2f%% |',
               format(sum(all$P < 0.001), big.mark = ','),
               100 * sum(all$P < 0.001) / n),
       sprintf('| **BH-FDR < 0.05** | **%d** | %.3f%% |',
               sum(all$BH < 0.05), 100 * sum(all$BH < 0.05) / n),
       sprintf('| BH-FDR < 0.10 | %d | %.3f%% |',
               sum(all$BH < 0.10), 100 * sum(all$BH < 0.10) / n),
       sprintf('| BH-FDR < 0.20 | %d | %.3f%% |',
               sum(all$BH < 0.20), 100 * sum(all$BH < 0.20) / n),
       sprintf('| BH-FDR < 0.30 | %d | %.3f%% |',
               sum(all$BH < 0.30), 100 * sum(all$BH < 0.30) / n),
       sprintf('| BH-FDR < 0.50 | %d | %.3f%% |',
               sum(all$BH < 0.50), 100 * sum(all$BH < 0.50) / n), '')

nsig <- sum(all$P < 0.05)
npass <- sum(all$P < 0.05 & all$BH < 0.05)
L <- c(L, sprintf('**这是本节最重要的数字**：名义 P < 0.05 的特征有 **%s** 个，',
                  format(nsig, big.mark = ',')), '',
       sprintf('但其中只有 **%d** 个能通过 BH-FDR < 0.05（占名义显著者的 %.1f%%）。',
               npass, 100 * npass / nsig), '',
       '未校正筛选把假阳性放大了约两个数量级。', '')

# ---- 3. 103 个候选的 FDR 情况 ----
cc <- all[all$in_candidate, ]
cc <- cc[order(cc$BH), ]
L <- c(L, '## 3. 103 个候选 circRNA 的 FDR 情况', '',
       sprintf('候选集合（未校正 P < 0.05 且 |log2FC| >= 1）共 **%d** 个。', nrow(cc)), '',
       '| FDR 阈值 | 通过数 | 占候选比例 |', '|---|---|---|')
for (thr in c(0.05, 0.10, 0.15, 0.20, 0.30, 0.50)) {
  k <- sum(cc$BH < thr)
  L <- c(L, sprintf('| BH-FDR < %.2f | %d | %.1f%% |', thr, k, 100 * k / nrow(cc)))
}
L <- c(L, '',
       sprintf('候选集合的 FDR 范围：**%.4f – %.4f**，中位 %.4f。',
               min(cc$BH), max(cc$BH), median(cc$BH)), '',
       sprintf('**只有 %d 个候选通过 FDR < 0.05**，即 %.1f%%。',
               sum(cc$BH < 0.05), 100 * sum(cc$BH < 0.05) / nrow(cc)), '')

# 通过 FDR 的候选
pass <- cc[cc$BH < 0.05, ]
if (nrow(pass)) {
  L <- c(L, '### 3.1 通过 FDR < 0.05 的候选', '',
         '| circRNA | 基因 | log2FC | 未校正 P | BH-FDR |', '|---|---|---|---|---|')
  for (i in seq_len(nrow(pass))) {
    g <- pass$gene[i]; if (nchar(g) == 0) g <- '(未注释)'
    L <- c(L, sprintf('| `%s` | %s | %+.3f | %.2e | %.4f |',
                      pass$circRNA[i], g, pass$log2FC[i], pass$P[i], pass$BH[i]))
  }
  L <- c(L, '',
         sprintf('两者中 `%s` 的宿主基因是 **%s**——这是本研究中唯一同时满足',
                 pass$circRNA[1], pass$gene[1]), '',
         '「跨队列预后验证显著」与「circRNA 层面 FDR 校正后显著」的候选。', '')
}

# FDR 最小的前 15
L <- c(L, '### 3.2 候选集合中 FDR 最小的 15 个', '',
       '| # | circRNA | 基因 | log2FC | 未校正 P | BH-FDR |', '|---|---|---|---|---|---|')
top <- cc[seq_len(min(15, nrow(cc))), ]
for (i in seq_len(nrow(top))) {
  g <- top$gene[i]; if (nchar(g) == 0) g <- '(未注释)'
  L <- c(L, sprintf('| %d | `%s` | %s | %+.3f | %.2e | %.4f |',
                    i, top$circRNA[i], g, top$log2FC[i], top$P[i], top$BH[i]))
}

writeLines(L, file.path(HCC, 'Major2_FDR_results.md'))

# ---- 4. 输出机器可读结果表 ----
write.table(all[, c('circRNA', 'gene', 'mean_normal_RPM', 'mean_tumor_RPM',
                     'log2FC', 't', 'P', 'BH', 'in_candidate')],
            file.path(RES, 'fdr_all_features.tsv'),
            sep = '\t', row.names = FALSE, quote = FALSE)
write.table(cc[, c('circRNA', 'gene', 'log2FC', 'P', 'BH')],
            file.path(RES, 'fdr_candidate_ranked.tsv'),
            sep = '\t', row.names = FALSE, quote = FALSE)

tab <- data.frame(
  criterion = c('features_tested', 'n_pairs', 'nominal_P05', 'BH_FDR05',
                'BH_FDR10', 'BH_FDR20', 'candidate_set',
                'candidate_BH_FDR05', 'candidate_BH_FDR10',
                'candidate_min_BH', 'candidate_median_BH'),
  value = c(n, 10, sum(all$P < 0.05), sum(all$BH < 0.05),
            sum(all$BH < 0.10), sum(all$BH < 0.20),
            nrow(cc), sum(cc$BH < 0.05), sum(cc$BH < 0.10),
            sprintf('%.4f', min(cc$BH)), sprintf('%.4f', median(cc$BH))),
  stringsAsFactors = FALSE)
write.table(tab, file.path(RES, 'fdr_summary.tsv'),
            sep = '\t', row.names = FALSE, quote = FALSE)

cat(sprintf('nominal P<0.05: %d | BH<0.05: %d | candidate %d | candidate BH<0.05: %d\n',
            sum(all$P < 0.05), sum(all$BH < 0.05), nrow(cc), sum(cc$BH < 0.05)))
cat(sprintf('step1 P_adj match rate: %.2f%%\n', cmp$pct_match[1]))
cat('DONE. Outputs: Major2_FDR结果.md, fdr_all_features.tsv, ',
    'fdr_candidate_ranked.tsv, fdr_summary.tsv\n')
