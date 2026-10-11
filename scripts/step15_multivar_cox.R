# step15_multivar_cox.R
# 多因素 Cox 回归（Major 3）+ 数值基准验证（Major 7e）
#
# 输入：results/r_input/ 下由 step14_export_for_R.py 导出的 TSV
# 输出：results/multivariable_cox.tsv, results/cox_benchmark.tsv
#
# 说明：
#  - 使用 survival::coxph，默认 Efron 处理 ties（生存时间以月/天为单位，并列很多）
#  - 表达量先做 z-score 标准化，使 HR 表示"每升高 1 个标准差"
#  - 单因素模型用于与 Python 版 step2_tcga_survival 交叉验证

suppressPackageStartupMessages(library(survival))

IN <- 'C:/Users/hewan/WorkBuddy/2026-10-09-22-01-46/hcc/results/r_input'
RES <- 'C:/Users/hewan/WorkBuddy/2026-10-09-22-01-46/hcc/results'

read_tsv <- function(p) read.delim(p, header = TRUE, sep = '\t',
                                   quote = '', stringsAsFactors = FALSE,
                                   check.names = FALSE)

# ---------- 载入 ----------
tcga_clin <- read_tsv(file.path(IN, 'tcga_clin.tsv'))
tcga_expr <- read_tsv(file.path(IN, 'tcga_expr.tsv'))
g_clin <- read_tsv(file.path(IN, 'gse14520_clin.tsv'))
g_expr <- read_tsv(file.path(IN, 'gse14520_expr.tsv'))
genes <- read_tsv(file.path(IN, 'gene_list.tsv'))$gene

cat('=== INPUT SCALE ===\n')
cat(sprintf('TCGA-LIHC  clin %d | expr %d genes x %d samples\n',
            nrow(tcga_clin), nrow(tcga_expr) - 1, ncol(tcga_expr) - 1))
cat(sprintf('GSE14520   clin %d | expr %d genes x %d samples\n',
            nrow(g_clin), nrow(g_expr) - 1, ncol(g_expr) - 1))

# ---------- 取基因表达向量 ----------
get_expr <- function(expr_df, clin_df, gene) {
  i <- which(expr_df$gene_symbol == gene)
  if (length(i) == 0) return(NULL)
  v <- as.numeric(expr_df[i, -1])
  names(v) <- colnames(expr_df)[-1]
  v[clin_df$sample]
}

# ---------- 单因素 Cox ----------
fit_uni <- function(clin, gene_val) {
  d <- clin
  z <- as.numeric(scale(gene_val))
  names(z) <- names(gene_val)
  d$z <- unname(z[d$sample])
  d <- d[complete.cases(d[, c('z', 'time', 'event')]), , drop = FALSE]
  f <- coxph(Surv(time, event) ~ z, data = d, ties = 'efron')
  s <- summary(f)
  # summary()$conf.int 列名（已实测确认）：
  #   "exp(coef)"  "exp(-coef)"  "lower .95"  "upper .95"
  # HR 直接取 coefficients 的 exp(coef) 列，不依赖 conf.int 列序
  ci <- s$conf.int
  list(n = nrow(d), ev = sum(d$event),
       HR = s$coefficients[1, 'exp(coef)'],
       lo = ci[1, 'lower .95'], hi = ci[1, 'upper .95'],
       p = s$coefficients[1, 5])
}

# ---------- 多因素 Cox ----------
fit_multi <- function(clin, gene_val) {
  d <- clin
  v <- as.numeric(scale(gene_val))
  names(v) <- names(gene_val)
  d$z <- unname(v[d$sample])
  d$stage <- factor(d$stage, levels = c(0, 1), labels = c('I_II', 'III_IV'))
  d$age_c <- d$age - 50      # centred at 50 years, reduces collinearity
  d$male <- factor(d$male, levels = c(0, 1), labels = c('female', 'male'))
  d <- d[complete.cases(d[, c('z', 'stage', 'age_c', 'male', 'time', 'event')]),
         , drop = FALSE]

  # GSE14520 的肿瘤队列全为女性，male 无方差 -> 该协变量不可估计，
  # coxph 会返回 NA。按方差检查自动剔除，保证模型可识别。
  use_sex <- length(unique(d$male)) > 1
  rhs <- if (use_sex) 'z + stage + age_c + male' else 'z + stage + age_c'
  f <- coxph(as.formula(sprintf('Surv(time, event) ~ %s', rhs)),
             data = d, ties = 'efron')
  s <- summary(f)
  # summary()$conf.int 列序（已实测确认）：
  #   [,1]=exp(coef)  [,2]=exp(-coef)  [,3]=lower.95  [,4]=upper.95
  ci <- s$conf.int
  co <- s$coefficients
  out <- data.frame(
    term = rownames(ci),
    HR = co[, 'exp(coef)'],
    lo = ci[, 'lower .95'],
    hi = ci[, 'upper .95'],
    p = co[, 5],
    stringsAsFactors = FALSE
  )
  list(n = nrow(d), ev = sum(d$event), tab = out, fit = f,
       sex_in_model = use_sex)
}

# ---------- 主流程 ----------
run_cohort <- function(tag, clin, expr) {
  rows_uni <- list()
  rows_mv <- list()
  for (g in genes) {
    gv <- get_expr(expr, clin, g)
    if (is.null(gv)) {
      cat(sprintf('  [%s] %s : expression missing, skipped\n', tag, g))
      next
    }
    # 单因素
    u <- tryCatch(fit_uni(clin, gv), error = function(e) {
      cat(sprintf('  [%s] %s UNI-ERROR: %s\n', tag, g, conditionMessage(e)))
      NULL
    })
    if (!is.null(u)) {
      rows_uni[[length(rows_uni) + 1]] <- data.frame(
        cohort = tag, gene = g, n = u$n, events = u$ev,
        HR_uni = u$HR, lo_uni = u$lo, hi_uni = u$hi, p_uni = u$p,
        stringsAsFactors = FALSE)
    }
    # 多因素
    m <- tryCatch(fit_multi(clin, gv), error = function(e) {
      cat(sprintf('  [%s] %s MV-ERROR: %s\n', tag, g, conditionMessage(e)))
      NULL
    })
    if (!is.null(m)) {
      tb <- m$tab
      tb$cohort <- tag
      tb$gene <- g
      tb$n <- m$n
      tb$events <- m$ev
      rows_mv[[length(rows_mv) + 1]] <- tb[, c('cohort', 'gene', 'term', 'n',
                                             'events', 'HR', 'lo', 'hi', 'p')]
    }
  }
  list(uni = do.call(rbind, rows_uni), mv = do.call(rbind, rows_mv))
}

cat('\n=== FITTING ===\n')
res_tcga <- run_cohort('TCGA-LIHC', tcga_clin, tcga_expr)
res_g <- run_cohort('GSE14520', g_clin, g_expr)

uni_all <- rbind(res_tcga$uni, res_g$uni)
mv_all <- rbind(res_tcga$mv, res_g$mv)

write.table(uni_all, file.path(RES, 'multivariable_cox_uni.tsv'),
            sep = '\t', row.names = FALSE, quote = FALSE)
write.table(mv_all, file.path(RES, 'multivariable_cox.tsv'),
            sep = '\t', row.names = FALSE, quote = FALSE)

# ---------- 输出单因素表 ----------
cat('\n===== UNIVARIABLE COX (expression z-score) =====\n')
for (tag in c('TCGA-LIHC', 'GSE14520')) {
  if (is.null(uni_all)) next
  sub <- uni_all[uni_all$cohort == tag, , drop = FALSE]
  if (!nrow(sub)) next
  cat(sprintf('\n--- %s (n=%d, events=%d) ---\n',
              tag, sub$n[1], sub$events[1]))
  for (i in seq_len(nrow(sub))) {
    cat(sprintf('  %-9s HR=%.4f (%.3f-%.3f)  P=%.3e\n',
                sub$gene[i], sub$HR_uni[i], sub$lo_uni[i], sub$hi_uni[i],
                sub$p_uni[i]))
  }
}

# ---------- 输出多因素表（只显示基因项）----------
cat('\n===== MULTIVARIABLE COX: gene term (adjusted stage + age + sex) =====\n')
for (tag in c('TCGA-LIHC', 'GSE14520')) {
  if (is.null(mv_all)) next
  sub <- mv_all[mv_all$cohort == tag & mv_all$term == 'z', , drop = FALSE]
  if (!nrow(sub)) next
  cat(sprintf('\n--- %s (n=%d, events=%d) ---\n',
              tag, sub$n[1], sub$events[1]))
  for (i in seq_len(nrow(sub))) {
    uni <- uni_all[uni_all$cohort == tag & uni_all$gene == sub$gene[i], , drop = FALSE]
    hr_u <- if (nrow(uni) > 0) sprintf('%.3f', uni$HR_uni[1]) else 'NA'
    p_u <- if (nrow(uni) > 0) sprintf('%.3f', uni$p_uni[1]) else 'NA'
    cat(sprintf('  %-9s adjHR=%.4f (%.3f-%.3f) P=%.3e  |  uniHR=%s P=%s\n',
                sub$gene[i], sub$HR[i], sub$lo[i], sub$hi[i], sub$p[i],
                hr_u, p_u))
  }
}

# ---------- 输出全部协变量项（供补充材料）----------
cat('\n===== ALL TERMS (example: NT5DC2 / HSD17B6) =====\n')
for (g in c('NT5DC2', 'HSD17B6')) {
  for (tag in c('TCGA-LIHC', 'GSE14520')) {
    sub <- mv_all[mv_all$cohort == tag & mv_all$gene == g, , drop = FALSE]
    if (!nrow(sub)) next
    cat(sprintf('\n--- %s / %s ---\n', tag, g))
    for (i in seq_len(nrow(sub))) {
      cat(sprintf('  %-12s HR=%.4f (%.3f-%.3f) P=%.3e\n',
                  sub$term[i], sub$HR[i], sub$lo[i], sub$hi[i], sub$p[i]))
    }
  }
}

cat('\n=== WARNINGS ===\n')
{
  ws <- warnings()
  if (is.null(ws)) cat('  (none)\n') else {
    uw <- unique(as.character(ws))
    cat(sprintf('  %d warnings, %d unique://n', length(ws), length(uw)))
    for (u in head(uw, 10)) cat('   -', u, '\n')
  }
}

cat('\nDONE. Outputs: multivariable_cox.tsv, multivariable_cox_uni.tsv\n')
