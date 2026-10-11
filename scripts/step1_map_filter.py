# -*- coding: utf-8 -*-
"""
HCC 全景 circRNA 网络构建 — 第 1 步
1) circRNA 坐标 → 基因映射（gencode v22）
2) 配对差异表达（肿瘤 vs 癌旁，10 对）
3) 候选 circRNA 筛选（未校正 P < 0.05 且 |log2FC| >= 1）
输出: results/circRNA_gene_map.tsv, results/candidate_circRNA.tsv
"""
import gzip, math, re, os, json

RAW = r'C:\Users\hewan\Documents\kimi\tasks\2026-10-08\11-38-08-f56eebff\rawdata'
GTF = os.path.join(RAW, 'annotation_cache', 'gencode.v22.annotation.gtf.gz')
MTX = 'results/GSE122482_circRNA_matrix.tsv'
OUT = 'results'
os.makedirs(OUT, exist_ok=True)

# ============ 1) 读基因区间 ============
print('[1/3] 读取 gencode v22 基因区间...')
genes = []          # (chr, start, end, strand, gene_name, biotype)
with gzip.open(GTF, 'rt', encoding='utf-8', errors='ignore') as f:
    for l in f:
        if l.startswith('#'):
            continue
        p = l.split('\t')
        if len(p) < 9 or p[2] != 'gene':
            continue
        m = re.search(r'gene_name "([^"]+)"', p[8])
        if not m:
            continue
        name = m.group(1)
        if name.startswith('ENSG'):
            continue
        bi = re.search(r'gene_type "([^"]+)"', p[8])
        genes.append((p[0], int(p[3]), int(p[4]), p[6], name,
                      bi.group(1) if bi else ''))
print(f'      基因数 {len(genes)}')

# 按染色体分组以加速
by_chr = {}
for g in genes:
    by_chr.setdefault(g[0], []).append(g)
for c in by_chr:
    by_chr[c].sort(key=lambda x: x[1])
print(f'      染色体数 {len(by_chr)}')

# ============ 2) 读 circRNA 表达矩阵并做配对差异 ============
print('\n[2/3] 读取 circRNA 矩阵并做配对差异分析...')
with open(MTX, encoding='utf-8') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    samples = hdr[1:]
# 前10=肿瘤(GSM3467136-145)，后10=癌旁(GSM3467146-155)
tum_idx = list(range(0, 10))
nor_idx = list(range(10, 20))
print(f'      肿瘤 {[samples[i] for i in tum_idx[:3]]}...')
print(f'      癌旁 {[samples[i] for i in nor_idx[:3]]}...')


def sf(t, df):
    if df <= 0 or math.isnan(t):
        return float('nan')
    x = df / (df + t * t)
    x = min(max(x, 1e-300), 1 - 1e-15)

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


PAT = re.compile(r'^hsa_circ:(chr[\w]+):(\d+)-(\d+)$')
res = []
with open(MTX, encoding='utf-8') as f:
    f.readline()
    for l in f:
        p = l.rstrip('\n').split('\t')
        name = p[0]
        try:
            v = [float(x) for x in p[1:]]
        except ValueError:
            continue
        tv = [v[i] for i in tum_idx]
        nv = [v[i] for i in nor_idx]
        diffs = [tv[k] - nv[k] for k in range(10)]
        md = sum(diffs) / len(diffs)
        var = sum((d - md) ** 2 for d in diffs) / (len(diffs) - 1)
        se = math.sqrt(var / len(diffs))
        t = md / se if se > 0 else 0.0
        pv = sf(t, 9)
        lfc = math.log2((sum(tv) / 10 + 1) / (sum(nv) / 10 + 1))
        res.append((name, sum(nv) / 10, sum(tv) / 10, lfc, t, pv))

print(f'      完成 {len(res)} 个 circRNA 配对检验')

# BH 校正
order = sorted(range(len(res)), key=lambda i: res[i][5])
m_ = len(res); adj = [1.0] * m_; prev = 1.0
for rank in range(m_, 0, -1):
    v = min(prev, res[order[rank - 1]][5] * m_ / rank)
    adj[order[rank - 1]] = v; prev = v
padj = {order[i]: adj[i] for i in range(m_)}

# ============ 3) 映射到基因 ============
print('\n[3/3] circRNA → 基因映射（区间重叠）...')
mapped = []
unmapped = 0
for i, r in enumerate(res):
    m = PAT.match(r[0])
    if not m:
        unmapped += 1
        continue
    ch, s, e = m.group(1), int(m.group(2)), int(m.group(3))
    cands = []
    lst = by_chr.get(ch, [])
    # 二分找第一个 start > e 的位置，再往前回溯
    lo, hi = 0, len(lst)
    while lo < hi:
        mid = (lo + hi) // 2
        if lst[mid][1] <= e:
            lo = mid + 1
        else:
            hi = mid
    j = lo - 1
    while j >= 0 and lst[j][2] >= s:
        gg = lst[j]
        if gg[2] >= s and lst[j][1] <= e:
            # 判定为外显子来源：circRNA 长度 < 50kb 且不完全覆盖整基因
            span = e - s
            overlap = min(e, gg[2]) - max(s, gg[1])
            if overlap > 0 and span < 500000:
                cands.append(gg[4])
        j -= 1
        if lo - j > 200:      # 安全上限，避免极端情况全扫
            break
    if cands:
        mapped.append((r, padj[i], cands[0]))
    else:
        mapped.append((r, padj[i], ''))

print(f'      成功映射 {sum(1 for x in mapped if x[2])} / {len(mapped)}')
print(f'      未映射 {unmapped}')

# ============ 输出 ============
with open(f'{OUT}/circRNA_gene_map.tsv', 'w', encoding='utf-8') as fo:
    fo.write('circRNA\tgene\tmean_normal_RPM\tmean_tumor_RPM\tlog2FC\tt\tP\tP_adj\n')
    for r, pa, g in sorted(mapped, key=lambda x: -abs(x[0][3])):
        fo.write('%s\t%s\t%.3f\t%.3f\t%.4f\t%.3f\t%.6g\t%.6g\n'
                 % (r[0], g, r[1], r[2], r[3], r[4], r[5], pa))
print(f'[SAVED] {OUT}/circRNA_gene_map.tsv')

# 候选：未校正 P<0.05 且 |log2FC|>=1 且 已映射到基因
cand = [x for x in mapped if x[2] and x[0][5] < 0.05 and abs(x[0][3]) >= 1.0]
with open(f'{OUT}/candidate_circRNA.tsv', 'w', encoding='utf-8') as fo:
    fo.write('circRNA\tgene\tmean_normal_RPM\tmean_tumor_RPM\tlog2FC\tP\tP_adj\tregulation\n')
    for r, pa, g in sorted(cand, key=lambda x: x[0][5]):
        fo.write('%s\t%s\t%.3f\t%.3f\t%.4f\t%.6g\t%.6g\t%s\n'
                 % (r[0], g, r[1], r[2], r[3], r[5], pa,
                    'up' if r[3] > 0 else 'down'))
up = sum(1 for x in cand if x[0][3] > 0)
print(f'\n[RESULT] 候选 circRNA（未校正 P<0.05, |log2FC|>=1, 已映射基因）: {len(cand)}')
print(f'         上调 {up}，下调 {len(cand)-up}')
print(f'[SAVED] {OUT}/candidate_circRNA.tsv')

# 涉及基因数
gset = {}
for r, pa, g in cand:
    gset.setdefault(g, []).append(r[0])
print(f'         涉及基因数: {len(gset)}')
print('\n=== 差异最显著的候选 circRNA（前 25）===')
print('%-40s %-14s %8s %8s %8s %10s' % ('circRNA', 'gene', 'normal', 'tumor', 'log2FC', 'P'))
for r, pa, g in sorted(cand, key=lambda x: x[0][5])[:25]:
    print('%-40s %-14s %8.2f %8.2f %8.3f %10.3g'
          % (r[0], g, r[1], r[2], r[3], r[5]))

json.dump({g: v for g, v in gset.items()},
          open(f'{OUT}/candidate_genes.json', 'w', encoding='utf-8'),
          ensure_ascii=False, indent=1)