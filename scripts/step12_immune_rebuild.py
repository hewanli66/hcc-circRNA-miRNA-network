# -*- coding: utf-8 -*-
"""
step12：重建免疫关联表（Major 4 / Major 6 修复）
两个问题：
  (1) 原immune_correlation_2queues.tsv 只筛了「方向一致」，未套方法声明的 |r|>0.15 阈值
      → 导致稿件写"seven genes"、实列6 个、表里其实 9 个，三方不一致
  (2) 稿件称使用 Charoentong 2017 特征基因集，实际是本项目自定义的 27×10 基因集
      → 必须如实报告，且不能挂 Charoentong 的引用

同时输出：
  - CD5 是否在任一特征集内（Major 4 循环论证核查）
  - 若在，则做 leave-one-out 重算（本次核查结果：不在，无需重算）

输出:
  results/immune_association_filtered.tsv   按方法阈值(|r|>0.15 双队列)筛出的关联
  results/immune_association_summary.tsv    每基因汇总
  data/immune_feature_sets_used.tsv        实际使用的特征基因集（供 Supplementary）
"""
import csv, os, json
from collections import defaultdict

OUT = 'results'

# ---------- 1) 实际使用的特征基因集（从脚本源码提取） ----------
import re
src = open('step4b_immune_tcga.py', encoding='utf-8').read()
m = re.search(r'IMMUNE_SETS\s*=\s*\{(.*?)\n\}', src, re.S)
sets = {}
for mm in re.finditer(r"'([\w]+)':\s*\[([^\]]+)\]", m.group(1)):
    sets[mm.group(1)] = [g.strip().strip("'\"") for g in mm.group(2).split(',')]

KEY = ['NT5DC2', 'HSD17B6', 'BCLAF1', 'SREBF2', 'CD5', 'LRCH1', 'RBM19', 'SMG5', 'NUP205']

print('[1/5] 特征基因集核查')
print('  细胞类型 %d 个，每类 %d 个基因' % (len(sets), len(list(sets.values())[0])))
membership = {}
for k in KEY:
    hit = [s for s, gs in sets.items() if k in gs]
    membership[k] = hit
    print('    %-9s → %s' % (k, hit if hit else '不在任何特征集内'))

with open('data/immune_feature_sets_used.tsv', 'w', encoding='utf-8', newline='') as fo:
    w = csv.writer(fo, delimiter='\t')
    w.writerow(['cell_type', 'n_genes', 'genes'])
    for s in sorted(sets):
        w.writerow([s, len(sets[s]), ';'.join(sets[s])])
print('  [SAVED] data/immune_feature_sets_used.tsv')

# ---------- 2) 读相关性表，套方法阈值 ----------
print('\n[2/5] 读相关性并套用 |r|>0.15（双队列）阈值')
ic = list(csv.DictReader(open(os.path.join(OUT, 'immune_correlation_2queues.tsv'),
                            encoding='utf-8'), delimiter='\t'))


def f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


TH = 0.15
kept = []
for r in ic:
    a = f(r['r_TCGA']); b = f(r['r_14520'])
    if a is None or b is None:
        continue
    if r['direction_consistent'] != 'YES':
        continue
    if abs(a) > TH and abs(b) > TH:          # 两队列都超阈值
        kept.append(dict(cell_type=r['cell_type'], gene=r['gene'],
                         r_TCGA=a, r_14520=b, absmean=(abs(a) + abs(b)) / 2))

print('  原表（仅方向一致）: %d 行 | 套阈值后: %d 行' % (len(ic), len(kept)))

with open(os.path.join(OUT, 'immune_association_filtered.tsv'), 'w',
          encoding='utf-8', newline='') as fo:
    w = csv.writer(fo, delimiter='\t')
    w.writerow(['gene', 'cell_type', 'r_TCGA', 'r_14520', 'mean_abs_r'])
    for r in sorted(kept, key=lambda x: -x['absmean']):
        w.writerow([r['gene'], r['cell_type'], '%.4f' % r['r_TCGA'],
                    '%.4f' % r['r_14520'], '%.4f' % r['absmean']])
print('  [SAVED] results/immune_association_filtered.tsv')

# ---------- 3) 汇总 ----------
print('\n[3/5] 每基因汇总')
byg = defaultdict(list)
for r in kept:
    byg[r['gene']].append(r)
rows = []
for g, rs in byg.items():
    rs.sort(key=lambda x: -x['absmean'])
    top = rs[0]
    rows.append(dict(gene=g, n_pairs=len(rs), top_cell=top['cell_type'],
                     top_r_tcga=top['r_TCGA'], top_r_14520=top['r_14520'],
                     mean_abs=sum(x['absmean'] for x in rs) / len(rs)))
rows.sort(key=lambda x: -x['top_r_tcga'] if abs(x['top_r_tcga']) > 0 else 0, reverse=True)
rows.sort(key=lambda x: -x['mean_abs'])

with open(os.path.join(OUT, 'immune_association_summary.tsv'), 'w',
          encoding='utf-8', newline='') as fo:
    w = csv.writer(fo, delimiter='\t')
    w.writerow(['gene', 'n_cell_pairs_passing', 'top_cell_type', 'r_TCGA', 'r_14520',
                'mean_abs_r', 'in_feature_set'])
    for r in rows:
        w.writerow([r['gene'], r['n_pairs'], r['top_cell'], '%.4f' % r['top_r_tcga'],
                    '%.4f' % r['top_r_14520'], '%.4f' % r['mean_abs'],
                    ';'.join(membership.get(r['gene'], [])) or 'no'])

print('  %-9s %-6s %-16s %-8s %-8s %s' % ('gene', 'pairs', 'top cell', 'r_TCGA', 'r_14520', 'mean|r|'))
for r in rows:
    print('  %-9s %-6d %-16s %+8.3f %+8.3f %.3f'
          % (r['gene'], r['n_pairs'], r['top_cell'], r['top_r_tcga'],
             r['top_r_14520'], r['mean_abs']))
print('  [SAVED] results/immune_association_summary.tsv')

# ---------- 4) 未通过阈值的基因（须从"稳定关联"中剔除） ----------
print('\n[4/5] 方向一致但不达阈值的基因（不得列入稳定关联）')
allg = set(r['gene'] for r in ic)
passed = set(byg)
for g in sorted(allg - passed):
    rs = [r for r in ic if r['gene'] == g]
    rs2 = [r for r in rs if f(r['r_TCGA']) is not None and f(r['r_14520']) is not None]
    rs2.sort(key=lambda x: -abs(f(x['r_TCGA'])))
    if rs2:
        t = rs2[0]
        print('  %-9s 最强 %-16s r_TCGA=%+.3f r_14520=%+.3f  → 剔除'
              % (g, t['cell_type'], f(t['r_TCGA']), f(t['r_14520'])))

# ---------- 5) Major 4 结论 ----------
print('\n[5/5] Major 4 循环论证核查')
cd5_hits = membership['CD5']
if cd5_hits:
    print('  CD5 在特征集内: %s → 需 leave-one-out 重算' % cd5_hits)
else:
    print('  ✅ CD5 不在任何免疫特征基因集内')
    print('     → 与 T 细胞浸润评分的相关性不是由共同基因导致的构造性相关')
    print('     → r = 0.877 不存在 Major 4 所指的循环论证')
    print('     → 但仍须在文中说明：bulk 组织中 CD5 信号主要反映免疫细胞组成，')
    print('        而非肿瘤细胞内源功能（审稿人要求的另一半）')
json.dump({'cd5_in_feature_set': bool(cd5_hits),
           'cd5_sets': cd5_hits,
           'feature_set_note': '27 cell types x 10 genes, curated for this study '
                               '(not the full Charoentong et al. 2017 sets)',
           'n_genes_passing_threshold': len(rows),
           'threshold': TH},
          open('results/immune_major4_check.json', 'w', encoding='utf-8'),
          ensure_ascii=False, indent=2)