# -*- coding: utf-8 -*-
"""
step4-b：免疫浸润在 TCGA-LIHC 上复现（GSE14520 之外第二个队列）
同一套ssGSEA 风格方法 + 同一套免疫基因集，保证跨队列可比
输出: results/immune_scores_TCGA.tsv, results/immune_correlation_2queues.tsv
"""
import gzip, os, csv, math, re

RAW = 'C:/Users/hewan/Documents/kimi/tasks/2026-10-08/11-38-08-f56eebff/rawdata'
OUT = 'results'
GTF = os.path.join(RAW, 'annotation_cache', 'gencode.v22.annotation.gtf.gz')

# 与 step4 完全一致的免疫基因集
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
KEY = ['NT5DC2', 'HSD17B6', 'BCLAF1', 'SREBF2', 'CD5', 'LRCH1', 'RBM19', 'SMG5', 'NUP205']

# ---------------- 读 Ensembl→symbol ----------------
print('[1/4] 读 gencode v22 映射...')
id2sym = {}
with gzip.open(GTF, 'rt', encoding='utf-8', errors='ignore') as f:
    for l in f:
        if l.startswith('#'): continue
        p = l.split('\t')
        if len(p) < 9 or p[2] != 'gene': continue
        mg = re.search(r'gene_id "([^"]+)"', p[8])
        mn = re.search(r'gene_name "([^"]+)"', p[8])
        if mg and mn:
            id2sym[mg.group(1).split('.')[0]] = mn.group(1)
print('      %d 条' % len(id2sym))

need = set()
for v in IMMUNE_SETS.values():
    need |= set(v)
need |= set(KEY)

# ---------------- 读 TCGA-LIHC FPKM ----------------
print('\n[2/4] 读 TCGA-LIHC FPKM（424 样本）...')
FPKM = os.path.join(RAW, 'TCGA-LIHC.star_fpkm.tsv.gz')
samples = []
sy = {}
with gzip.open(FPKM, 'rt', encoding='utf-8', errors='ignore') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    samples = hdr[1:]
    for l in f:
        p = l.rstrip('\n').split('\t')
        gid = p[0].split('.')[0]
        s = id2sym.get(gid)
        if not s or s not in need: continue
        try:
            sy[s] = [float(x) for x in p[1:1 + len(samples)]]
        except ValueError:
            continue
print('      命中 symbol %d / %d' % (len(sy), len(need)))
missing = [g for g in KEY if g not in sy]
if missing: print('      缺失关键基因: %s' % missing)

# 样本分组（barcode 末段 01A/11A）
def stype(b):
    seg = b.split('-')[-1]
    if not seg: return 'other'
    if seg[0] == '0': return 'tumor'
    if seg[0] == '1' and len(seg) > 1 and seg[1].isdigit(): return 'normal'
    return 'other'

tum = [i for i, s in enumerate(samples) if stype(s) == 'tumor']
print('      肿瘤 %d' % len(tum))
n = len(tum)
tum_samples = [samples[i] for i in tum]

# ---------------- ssGSEA ----------------
print('\n[3/4] 计算 ssGSEA 免疫分数（%d 肿瘤样本）...' % n)
expr = {g: [sy[g][i] for i in tum] for g in sy}

gene_names = sorted(expr.keys())
ranks = {}
for g in gene_names:
    vals = expr[g]
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
    if len(use) < 3: continue
    scores[cell] = [sum(ranks[g][i] for g in use) / len(use) for i in range(n)]
print('      完成 %d 个细胞类型' % len(scores))

with open(f'{OUT}/immune_scores_TCGA.tsv', 'w', encoding='utf-8', newline='') as fo:
    wr = csv.writer(fo, delimiter='\t')
    wr.writerow(['sample', 'cell_type', 'score'])
    for cell, sc in scores.items():
        for i, s in enumerate(tum_samples):
            wr.writerow([s, cell, '%.4f' % sc[i]])
print('[SAVED] %s/immune_scores_TCGA.tsv' % OUT)

# ---------------- 相关性 ----------------
print('\n[4/4] 相关性分析...')
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

gz = {}
for g in KEY:
    if g not in expr: continue
    v = expr[g]; m = sum(v) / n
    sd = (sum((x - m) ** 2 for x in v) / (n - 1)) ** 0.5
    if sd > 0: gz[g] = [(x - m) / sd for x in v]

# 读 GSE14520 的结果做对照
prev = {}
if os.path.exists(f'{OUT}/immune_correlation_14520.tsv'):
    with open(f'{OUT}/immune_correlation_14520.tsv', encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            prev[(r['cell_type'], r['gene'])] = float(r['r'])

with open(f'{OUT}/immune_correlation_2queues.tsv', 'w', encoding='utf-8', newline='') as fo:
    wr = csv.writer(fo, delimiter='\t')
    wr.writerow(['cell_type', 'gene', 'r_TCGA', 'P_TCGA', 'r_14520', 'direction_consistent'])
    for cell in sorted(scores):
        sc = scores[cell]
        for g in gz:
            r = pearson(sc, gz[g])
            if math.isnan(r): continue
            t = r * math.sqrt((n - 2) / max(1e-9, 1 - r * r))
            p = sf(t, n - 2)
            pr = prev.get((cell, g))
            cons = '' if pr is None else ('YES' if (r > 0) == (pr > 0) else 'no')
            wr.writerow([cell, g, '%.4f' % r, '%.6g' % p, pr if pr is not None else '', cons])
print('[SAVED] %s/immune_correlation_2queues.tsv' % OUT)

# 汇总
rows = []
with open(f'{OUT}/immune_correlation_2queues.tsv', encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        rows.append(r)
both = [r for r in rows if r['direction_consistent'] == 'YES'
        and r['r_14520'] and abs(float(r['r_14520'])) > 0.2]
cons = [r for r in rows if r['direction_consistent'] == 'YES']
print('\n===== 双队列方向一致性 =====')
print('可比关联数: %d，方向一致: %d (%.0f%%)'
      % (len(rows), len(cons), len(cons) / max(len(rows), 1) * 100))
print('\n===== 两队列 |r|>0.2 且方向一致（P<0.05 标*）=====')
for r in sorted(both, key=lambda x: -abs(float(x['r_TCGA']))):
    p = float(r['P_TCGA'])
    star = '*' if p < 0.05 else ' '
    print('%s %-16s vs %-10s  TCGA r=%+.3f (P=%.4f) | 14520 r=%+.3f'
          % (star, r['cell_type'], r['gene'], float(r['r_TCGA']), p, float(r['r_14520'])))