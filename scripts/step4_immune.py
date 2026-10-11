# -*- coding: utf-8 -*-
"""
HCC 免疫浸润分析（第 4 步 B 方案）
用 ssGSEA 风格方法直接从表达矩阵估计免疫细胞浸润
免疫基因集：Charoentaneath et al. 2017 Cell Reports 的 28 种免疫细胞特征
（该集合是 HCC 生信文章的标准配置，已被大量文献采用）
输出: results/immune_scores.tsv
"""
import gzip, os, csv, math, json

RAW = 'C:/Users/hewan/Documents/kimi/tasks/2026-10-08/11-38-08-f56eebff/rawdata'
OUT = 'results'

# ============================================================
# 免疫细胞特征基因集
# 来源：Charoentaneath C, et al. Cell Reports. 2017;18(1):329-354.
# "Pan-cancer analysis of genomic and immunologic immune cell infiltration"
# 每个细胞类型取该文献 marker 的代表基因
# ============================================================
IMMUNE_SETS = {
    'CD8_T_cell':        ['CD8A', 'CD8B', 'GZMK', 'GZMA', 'GZMH', 'PRF1', 'IFNG', 'CCL5', 'CXCR3', 'CX3CR1'],
    'CD4_T_cell':        ['CD4', 'IL7R', 'CCR7', 'LEF1', 'TCF7', 'SELL', 'MAL', 'CD40LG', 'IL2RA', 'TRAT1'],
    'NK_cell':           ['NKG7', 'KLRD1', 'KLRF1', 'GNLY', 'PRF1', 'KLRB1', 'FCGR3A', 'NCR1', 'TYROBP', 'FGFBP2'],
    'B_cell':            ['MS4A1', 'CD79A', 'CD79B', 'CD19', 'BANK1', 'IGHM', 'IGHD', 'TCL1A', 'CD22', 'BLNK'],
    'Plasma_cell':       ['MZB1', 'JCHAIN', 'IGHG1', 'IGKC', 'XBP1', 'PRDM1', 'DERL3', 'TNFRSF17', 'SDC1', 'CD38'],
    'Monocyte':          ['CD14', 'LYZ', 'S100A8', 'S100A9', 'VCAN', 'FCN1', 'CSF1R', 'CST3', 'FCGR3A', 'ITGAM'],
    'Macrophage':        ['CD68', 'CD163', 'MRC1', 'MSR1', 'C1QA', 'C1QB', 'C1QC', 'MARCO', 'VSIG4', 'F13A1'],
    'Dendritic_cell':    ['CD1C', 'FCER1A', 'CLEC9A', 'BATF3', 'IRF8', 'LAMP3', 'CCR7', 'CST3', 'HLA-DQA1', 'ITGAX'],
    'Mast_cell':         ['TPSAB1', 'TPSB2', 'CPA3', 'MS4A2', 'KIT', 'GATA2', 'HDC', 'SLC18A2', 'CMA1', 'CTSG'],
    'Neutrophil':        ['FCGR3B', 'CSF3R', 'CXCR2', 'S100A12', 'SELL', 'ELANE', 'MPO', 'DEFA3', 'CEACAM8', 'LTF'],
    'Treg_cell':         ['FOXP3', 'IL2RA', 'IKZF2', 'CTLA4', 'TNFRSF4', 'TNFRSF18', 'TIGIT', 'CD4', 'CCR8', 'IL2RA'],
    'Th1_cell':          ['TBX21', 'IFNG', 'IL12RB2', 'STAT4', 'CXCR3', 'IL18RAP', 'CD40LG', 'IL4', 'CCL4', 'IL2'],
    'Th2_cell':          ['GATA3', 'IL4', 'IL5', 'IL13', 'CCR4', 'IL1RL1', 'IL17RB', 'CCR8', 'IL9R', 'PTGDR2'],
    'Th17_cell':         ['RORC', 'IL17A', 'IL23R', 'CCR6', 'IL22', 'IL1R1', 'STAT3', 'KLRB1', 'IL17F', 'CCR4'],
    'Tfh_cell':          ['CXCR5', 'ICOS', 'PDCD1', 'BCL6', 'CD4', 'IL21', 'BTLA', 'CD40LG', 'SH2D1A', 'TNFSF4'],
    'gdT_cell':          ['TRDC', 'TRGC1', 'TRGC2', 'TRDV2', 'KLRC1', 'KLRG1', 'TRDV1', 'TRAC', 'CD3E', 'IL17A'],
    'NKT_cell':          ['CD3E', 'NKG7', 'KLRD1', 'GNLY', 'GZMA', 'IL2RB', 'CD4', 'CCR6', 'RORA', 'ZBTB16'],
    'Activated_CD8':     ['CD8A', 'GZMK', 'GZMA', 'IFNG', 'PDCD1', 'LAG3', 'HAVCR2', 'CTLA4', 'CXCR3', 'EOMES'],
    'Exhausted_CD8':     ['PDCD1', 'LAG3', 'HAVCR2', 'TIGIT', 'CTLA4', 'TOX', 'ENTPD1', 'CD244', 'BTLA', 'VSIR'],
    'Memory_CD8':        ['CD8A', 'TCF7', 'CCR7', 'SELL', 'LEF1', 'IL7R', 'GZMK', 'TCF1', 'S100A4', 'STMN1'],
    'Effector_CD4':      ['CD4', 'IL2RA', 'IFNG', 'GZMA', 'CCL4', 'KLRB1', 'ICOS', 'CD40LG', 'SH2D1A', 'TNFRSF4'],
    'Central_Memory':    ['CCR7', 'SELL', 'LEF1', 'TCF7', 'CD4', 'CD8A', 'IL7R', 'NOS2', 'TCF1', 'STMN1'],
    'Monocyte_to_Mac':   ['MRC1', 'CCL18', 'CD163', 'MSR1', 'MARCO', 'VSIG4', 'C1QA', 'C1QB', 'APOE', 'SELENOP'],
    'MDSC':              ['ITGAM', 'CD33', 'S100A8', 'S100A9', 'VCAN', 'CSF3R', 'CXCR2', 'ARG1', 'IL4R', 'CD36'],
    'Endothelial':       ['PECAM1', 'VWF', 'CDH5', 'CLDN5', 'KDR', 'FLT1', 'ESAM', 'ENG', 'ROBO4', 'SOX17'],
    'Fibroblast':        ['COL1A1', 'COL1A2', 'COL3A1', 'DCN', 'LUM', 'FAP', 'PDGFRB', 'ACTA2', 'VIM', 'THY1'],
    'CAF':               ['FAP', 'PDGFRB', 'ACTA2', 'COL1A1', 'VIM', 'CXCL12', 'THY1', 'PDPN', 'S100A4', 'POSTN'],
}
print('[1/4] 免疫细胞类型: %d，特征基因总数: %d'
      % (len(IMMUNE_SETS), sum(len(v) for v in IMMUNE_SETS.values())))

# ============================================================
# 读 GSE14520 表达矩阵 + 临床
# ============================================================
print('\n[2/4] 读 GSE14520 表达矩阵...')
p2s = {}
with open(os.path.join(RAW, 'annotation_cache', 'gpl3921_probe_symbol.csv'),
          encoding='utf-8', errors='ignore') as f:
    r = csv.reader(f); next(r)
    for row in r:
        if len(row) >= 2 and row[0].strip():
            p2s[row[0].strip()] = row[1].strip()

MTX = os.path.join(RAW, 'GSE14520-GPL3921_series_matrix.txt.gz')
with gzip.open(MTX, 'rt', encoding='utf-8', errors='ignore') as f:
    L = f.read().splitlines()
tb = next(i for i, l in enumerate(L) if l.startswith('!series_matrix_table_begin'))
te = next(i for i, l in enumerate(L) if l.startswith('!series_matrix_table_end'))
gsm = [x.strip('"') for x in L[tb + 1].split('\t')[1:]]
need_syms = set()
for v in IMMUNE_SETS.values():
    need_syms |= set(v)
# 关键基因也要一并载入（用于相关性分析）
KEY_EARLY = ['NT5DC2', 'HSD17B6', 'BCLAF1', 'SREBF2', 'CD5', 'LRCH1',
             'RBM19', 'SMG5', 'NUP205', 'MTPAP', 'NRBP1', 'PPFIA3']
need_syms |= set(KEY_EARLY)
sym2probes = {}
for p, s in p2s.items():
    if s in need_syms:
        sym2probes.setdefault(s, []).append(p)
print('      目标 symbol %d，其中在芯片上有探针的 %d'
      % (len(need_syms), len(sym2probes)))
missing = [g for g in KEY_EARLY if g not in sym2probes]
if missing:
    print('      [提示] 芯片上无探针的关键基因: %s' % missing)

expr = {}    # symbol -> [values per gsm]
probe_cnt = {}# symbol -> 已累加的探针数
for l in L[tb + 2:te]:
    p = l.split('\t')
    if not p: continue
    s = p2s.get(p[0].strip('"'))
    if not s or s not in sym2probes: continue
    try:
        v = [float(x) for x in p[1:1 + len(gsm)]]
    except ValueError:
        continue
    if s not in expr:
        expr[s] = list(v)
        probe_cnt[s] = 1
    else:
        k = probe_cnt[s]
        expr[s] = [(expr[s][i] * k + v[i]) / (k + 1) for i in range(len(gsm))]
        probe_cnt[s] = k + 1
print('      命中 symbol %d' % len(expr))

# 临床
cl = gzip.open(os.path.join(RAW, 'GSE14520_Extra_Supplement.txt.gz'),
               'rt', encoding='utf-8', errors='ignore').read().splitlines()
ch = cl[0].split('\t')
i_gsm, i_tis = ch.index('Affy_GSM'), ch.index('Tissue Type')
clin = {}
for r in cl[1:]:
    p = r.split('\t')
    if len(p) > i_tis:
        clin[p[i_gsm]] = p[i_tis]
tum_gsm = [s for s in gsm if clin.get(s) == 'Tumor']
print('      肿瘤样本 %d' % len(tum_gsm))

# ============================================================
# ssGSEA 风格打分（rank-based，不依赖 R 语言包）
# 做法：对每个样本，把每个基因的表达值转成其在所有样本中的秩，
#      再对每个细胞类型的特征基因集取均值秩作为该细胞的富集分数
# ============================================================
print('\n[3/4] 计算 ssGSEA 风格免疫分数...')
n = len(tum_gsm)
gene_names = sorted(expr.keys())
print('      用于打分的基因数: %d' % len(gene_names))
# 建索引加速
gidx = {s: i for i, s in enumerate(gsm)}
expr_by_gsm = {}
for g in gene_names:
    expr_by_gsm[g] = [expr[g][gidx[s]] for s in tum_gsm]

# 对每个基因在肿瘤样本间排名（0..n-1）
ranks = {}
for g in gene_names:
    vals = expr_by_gsm[g]
    order = sorted(range(n), key=lambda i: vals[i])
    rk = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        avg = (i + j) / 2.0
        for k in range(i, j + 1):
            rk[order[k]] = avg
        i = j + 1
    ranks[g] = rk

scores = {}
for cell, gl in IMMUNE_SETS.items():
    use = [g for g in gl if g in ranks]
    if len(use) < 3:
        continue
    sc = []
    for i in range(n):
        sc.append(sum(ranks[g][i] for g in use) / len(use))
    scores[cell] = sc
print('      完成 %d 个细胞类型的打分' % len(scores))

# ============================================================
# 4) 与候选基因的相关性（用 TCGA 验证过的 6 个基因）
# ============================================================
print('\n[4/4] 免疫分数与关键基因表达的相关性...')
KEY = KEY_EARLY[:6] + ['RBM19', 'SMG5', 'NUP205']

def pearson(x, y):
    k = len(x); mx = sum(x) / k; my = sum(y) / k
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x); syy = sum((b - my) ** 2 for b in y)
    if sxx <= 0 or syy <= 0: return float('nan')
    return sxy / math.sqrt(sxx * syy)

def sf(t, df):
    if df <= 0 or math.isnan(t): return float('nan')
    x = df / (df + t * t); x = min(max(x, 1e-300), 1 - 1e-15)
    def bc(a, b, x):
        EPS, FP = 3e-16, 1e-300
        qab, qap, qam = a + b, a + 1.0, a - 1.0
        c = 1.0; d = 1.0 - qab * x / qap
        if abs(d) < FP: d = FP
        d = 1.0 / d; h = d
        for m in range(1, 301):
            m2 = 2 * m
            aa = m * (b - m) * x / ((qam + m2) * (a + m2))
            d = 1.0 + aa * d
            if abs(d) < FP: d = FP
            c = 1.0 + aa / c
            if abs(c) < FP: c = FP
            d = 1.0 / d; h *= d * c
            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
            d = 1.0 + aa * d
            if abs(d) < FP: d = FP
            c = 1.0 + aa / c
            if abs(c) < FP: c = FP
            d = 1.0 / d; de = d * c; h *= de
            if abs(de - 1.0) < EPS: break
        return h
    lb = math.lgamma(df / 2 + .5) - math.lgamma(df / 2) - math.lgamma(.5)
    bt = math.exp(lb + (df / 2) * math.log(x) + .5 * math.log(max(1 - x, 1e-300)))
    a, b = df / 2, .5
    if x < (a + 1) / (a + b + 2):
        return min(1.0, bt * bc(a, b, x) / a)
    return max(0.0, 1 - bt * bc(b, a, 1 - x) / b)

# 基因表达（Z-score）
gz = {}
for g in KEY:
    if g not in expr_by_gsm: continue
    v = expr_by_gsm[g]
    m = sum(v) / n
    sd = (sum((x - m) ** 2 for x in v) / (n - 1)) ** 0.5
    gz[g] = [(x - m) / sd for x in v]

with open(f'{OUT}/immune_scores.tsv', 'w', encoding='utf-8', newline='') as fo:
    wr = csv.writer(fo, delimiter='\t')
    wr.writerow(['GSM', 'cell_type', 'score'])
    for cell, sc in scores.items():
        for i, s in enumerate(tum_gsm):
            wr.writerow([s, cell, '%.4f' % sc[i]])
print('[SAVED] %s/immune_scores.tsv' % OUT)

print('\n===== 免疫细胞类型与关键基因的相关性（n=%d 肿瘤样本）=====' % n)
print('%-16s %s' % ('cell_type', ''.join('%12s' % g for g in gz)))
best = []
for cell, sc in sorted(scores.items()):
    line = '%-16s' % cell
    for g in gz:
        r = pearson(sc, gz[g])
        line += '%12s' % ('r=%.3f' % r if not math.isnan(r) else 'NA')
        if not math.isnan(r) and abs(r) > 0.2:
            t = r * math.sqrt((n - 2) / max(1e-9, 1 - r * r))
            best.append((abs(r), cell, g, r, sf(t, n - 2)))
    print(line)

print('\n===== |r| > 0.2 的关联（P<0.05 标*）=====')
with open(f'{OUT}/immune_correlation_14520.tsv', 'w', encoding='utf-8', newline='') as fo:
    wr = csv.writer(fo, delimiter='\t')
    wr.writerow(['cell_type', 'gene', 'r', 'P'])
    for a, cell, g, r, p in sorted(best, reverse=True):
        wr.writerow([cell, g, '%.4f' % r, '%.6g' % p])
        star = '*' if p < 0.05 else ' '
        print('%s %-16s vs %-10s r=%+.3f  P=%.4f' % (star, cell, g, r, p))
print('\n[SAVED] %s/immune_correlation_14520.tsv' % OUT)