# step16_cox_benchmark.R
# Major 7(e) 数值基准验证
#
# 目的：用 R survival::coxph 重算稿件现有的单因素 Cox P 值，
#      与 step2_tcga_survival.py（Python 自实现）的结果逐项比对，
#      证明 P 值可复现、不依赖特定软件实现。
#
# 关键：必须用与稿件完全相同的口径 —— TCGA 用原始 FPKM（未取对数、
#      未标准化），GSE14520 用原始探针强度。这样 HR 与 P 才可直接比较。
#
# 输出：results/cox_benchmark.tsv

suppressPackageStartupMessages(library(survival))

IN <- 'C:/Users/hewan/WorkBuddy/2026-10-09-22-01-46/hcc/results/r_input'
RES <- 'C:/Users/hewan/WorkBuddy/2026-10-09-22-01-46/hcc/results'

read_tsv <- function(p) read.delim(p, header = TRUE, sep = '\t',
                                   quote = '', stringsAsFactors = FALSE,
                                   check.names = FALSE)

tcga_clin <- read_tsv(file.path(IN, 'tcga_clin.tsv'))
tcga_expr <- read_tsv(file.path(IN, 'tcga_expr_raw.tsv'))   # 原始 FPKM
g_clin <- read_tsv(file.path(IN, 'gse14520_clin.tsv'))
g_expr <- read_tsv(file.path(IN, 'gse14520_expr_raw.tsv'))
genes <- read_tsv(file.path(IN, 'gene_list.tsv'))$gene

get_expr <- function(expr_df, clin_df, gene) {
  i <- which(expr_df$gene_symbol == gene)
  if (length(i) == 0) return(NULL)
  v <- as.numeric(expr_df[i, -1])
  names(v) <- colnames(expr_df)[-1]
  v[clin_df$sample]
}

# 稿件口径的单因素 Cox：原始表达值（每单位变化）
fit_raw <- function(clin, gene_val) {
  d <- clin
  v <- as.numeric(gene_val)
  names(v) <- names(gene_val)
  d$x <- unname(v[d$sample])
  d <- d[complete.cases(d[, c('x', 'time', 'event')]), , drop = FALSE]
  f <- coxph(Surv(time, event) ~ x, data = d, ties = 'efron')
  s <- summary(f)
  list(n = nrow(d), ev = sum(d$event),
       HR = s$coefficients[1, 'exp(coef)'],
       lo = s$conf.int[1, 'lower .95'],
       hi = s$conf.int[1, 'upper .95'],
       p = s$coefficients[1, 5])
}

run <- function(tag, clin, expr) {
  rows <- list()
  for (g in genes) {
    gv <- get_expr(expr, clin, g)
    if (is.null(gv)) next
    r <- tryCatch(fit_raw(clin, gv), error = function(e) {
      cat(sprintf('  [%s] %s ERROR: %s\n', tag, g, conditionMessage(e)))
      NULL
    })
    if (!is.null(r)) {
      rows[[length(rows) + 1]] <- data.frame(
        cohort = tag, gene = g, n = r$n, events = r$ev,
        HR_R = r$HR, lo_R = r$lo, hi_R = r$hi, p_R = r$p,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

cat('=== BENCHMARK: R coxph on raw expression (same scale as manuscript) ===\n')
b1 <- run('TCGA-LIHC', tcga_clin, tcga_expr)
b2 <- run('GSE14520', g_clin, g_expr)
bench <- rbind(b1, b2)
write.table(bench, file.path(RES, 'cox_benchmark.tsv'),
            sep = '\t', row.names = FALSE, quote = FALSE)

for (tag in c('TCGA-LIHC', 'GSE14520')) {
  sub <- bench[bench$cohort == tag, , drop = FALSE]
  if (!nrow(sub)) next
  cat(sprintf('\n--- %s (n=%d, events=%d) ---\n', tag, sub$n[1], sub$events[1]))
  for (i in seq_len(nrow(sub))) {
    cat(sprintf('  %-9s HR=%.4f (%.3f-%.3f)  P=%.4e\n',
                sub$gene[i], sub$HR_R[i], sub$lo_R[i], sub$hi_R[i],
                sub$p_R[i]))
  }
}
cat('\nDONE. Output: cox_benchmark.tsv\n')
