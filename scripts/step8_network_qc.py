# -*- coding: utf-8 -*-
"""
第8 步：网络质量筛选 + circRNA-miRNA-mRNA 三元网络构建
质量维度（circBank 对人源条目未填 highConfidence/validated/conserved，实测 95 条全空，
故只能用这三项+ 结合自由能）：
  1) match_type：exact > loose > gene_fallback（匹配可信度）
  2) 结合自由能 energy（越负越强，阈值 -20 / -22 / -25 kcal/mol）
  3) circBank score（越高越可信）
分层输出：high / medium / low 三档网络
输出: results/network_edges_qc.tsv, results/trinet_nodes.tsv, results/trinet_summary.tsv
"""
import csv, json, os, math
from collections import defaultdict, Counter

OUT = 'results'

# ---------- 载入 ----------
edges_raw = list(csv.DictReader(open(f'{OUT}/candidate_circRNA_miRNA.tsv', encoding='utf-8'),
                               delimiter='\t'))
mapinfo = {r['circRNA_GSE']: r for r in
           csv.DictReader(open(f'{OUT}/circbank_id_map.tsv', encoding='utf-8'), delimiter='\t')}
genes = {r['gene']: r for r in
         csv.DictReader(open(f'{OUT}/network_genes.tsv', encoding='utf-8'), delimiter='\t')}
imm = defaultdict(list)
p2 = f'{OUT}/immune_correlation_2queues.tsv'
if os.path.exists(p2):
    for r in csv.DictReader(open(p2, encoding='utf-8'), delimiter='\t'):
        if r['direction_consistent'] == 'YES':
            try:
                imm[r['gene']].append((r['cell_type'], float(r['r_TCGA'])))
            except (ValueError, KeyError):
                pass

print('原始边: %d' % len(edges_raw))


def fe(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def fs(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


# ---- 去重：circBank 分页边界会返回重叠记录（实测多余 766 行 / 5%）----
_seen_edge = set()
_dedup = []
_dup_n = 0
for r in edges_raw:
    k = (r['circRNA_GSE'], r['host_gene'], r['miRNA'],
         r['energy_kcal_mol'], r['score'])
    if k in _seen_edge:
        _dup_n += 1
        continue
    _seen_edge.add(k)
    _dedup.append(r)
print('去重: %d → %d 行（移除重复 %d 行，%.1f%%）'
      % (len(edges_raw), len(_dedup), _dup_n, _dup_n / len(edges_raw) * 100))
edges_raw = _dedup


# ---------- 质量打分 ----------
# match_type 权重
MT_W = {'exact': 3.0, 'loose': 2.0, 'gene_fallback': 1.0}
rows = []
for r in edges_raw:
    mt = r['match_type']
    e = fe(r['energy_kcal_mol'])
    s = fs(r['score'])
    w_mt = MT_W.get(mt, 1.0)
    # energy 归一：-15 →0 分，-30 → 1 分
    w_e = 0.0
    if e is not None:
        w_e = max(0.0, min(1.0, (-e - 15.0) / 15.0))
    # score 归一：140 → 0，192 → 1
    w_s = 0.0
    if s is not None:
        w_s = max(0.0, min(1.0, (s - 140.0) / 52.0))
    # 总分 0~5
    q = w_mt + w_e * 1.5 + w_s * 0.5
    rows.append(dict(circ=r['circRNA_GSE'], cb=r['circBankID'], gene=r['host_gene'],
                     mt=mt, miRNA=r['miRNA'], energy=e, score=s,
                     log2fc=float(r['circ_log2FC']), reg=r['circ_reg'],
                     q=round(q, 3)))

print('打分完成，%d 条' % len(rows))

# ---------- 分档 ----------
# HIGH: exact 且 energy<=-20；MEDIUM: (exact|loose) 且 energy<=-20；LOW: 其余
def tier(r):
    if r['mt'] == 'exact' and r['energy'] is not None and r['energy'] <= -20:
        return 'high'
    if r['mt'] in ('exact', 'loose') and r['energy'] is not None and r['energy'] <= -20:
        return 'medium'
    return 'low'


for r in rows:
    r['tier'] = tier(r)

tc = Counter(r['tier'] for r in rows)
print('\n分层结果:')
for k in ('high', 'medium', 'low'):
    sub = [r for r in rows if r['tier'] == k]
    print('  %-7s %6d 边 | %5d miRNA | %3d circRNA | %3d 基因'
          % (k, len(sub), len(set(x['miRNA'] for x in sub)),
             len(set(x['circ'] for x in sub)), len(set(x['gene'] for x in sub))))

with open(f'{OUT}/network_edges_qc.tsv', 'w', encoding='utf-8', newline='') as fo:
    w = csv.writer(fo, delimiter='\t')
    w.writerow(['circRNA', 'circBankID', 'host_gene', 'match_type', 'miRNA',
                'energy_kcal_mol', 'score', 'quality_score', 'tier', 'circ_log2FC', 'circ_reg'])
    for r in sorted(rows, key=lambda x: -x['q']):
        w.writerow([r['circ'], r['cb'], r['gene'], r['mt'], r['miRNA'],
                    r['energy'] if r['energy'] is not None else '',
                    r['score'] if r['score'] is not None else '',
                    '%.3f' % r['q'], r['tier'], '%.4f' % r['log2fc'], r['reg']])
print('\n[SAVED] %s/network_edges_qc.tsv' % OUT)

# ---------- 构建三元网络节点 ----------
# miRNA 节点：连接度 + 平均能量
mir_info = defaultdict(lambda: {'deg': 0, 'energies': [], 'circs': set(), 'genes': set()})
for r in rows:
    m = mir_info[r['miRNA']]
    m['deg'] += 1
    if r['energy'] is not None:
        m['energies'].append(r['energy'])
    m['circs'].add(r['circ'])
    m['genes'].add(r['gene'])

print('\n===== miRNA 节点 hub（degree>=8）=====')
hubs = sorted(mir_info.items(), key=lambda kv: -kv[1]['deg'])
for name, m in hubs[:15]:
    avg_e = sum(m['energies']) / len(m['energies']) if m['energies'] else 0
    print('  %-20s deg=%-4d 平均能量=%6.1f  连接基因=%d'
          % (name, m['deg'], avg_e, len(m['genes'])))

# circRNA 节点：连接度 + 所属基因的预后属性
circ_info = defaultdict(lambda: {'mirnas': set(), 'gene': '', 'log2fc': 0, 'reg': ''})
for r in rows:
    c = circ_info[r['circ']]
    c['mirnas'].add(r['miRNA'])
    c['gene'] = r['gene']
    c['log2fc'] = r['log2fc']
    c['reg'] = r['reg']

# ---------- 输出三元网络（分层） ----------
val = json.load(open(f'{OUT}/validated_genes.json')) if os.path.exists(f'{OUT}/validated_genes.json') else []

with open(f'{OUT}/trinet_nodes.tsv', 'w', encoding='utf-8', newline='') as fo:
    w = csv.writer(fo, delimiter='\t')
    w.writerow(['node_type', 'node_id', 'degree', 'avg_energy', 'host_gene',
                'log2FC', 'regulation', 'HR_TCGA', 'P_cox_TCGA', 'cross_cohort',
                'immune_max_r', 'validated_gene', 'tiers'])
    for name, m in hubs:
        if m['deg'] < 2:
            continue
        avg_e = sum(m['energies']) / len(m['energies']) if m['energies'] else 0
        # 该 miRNA 连接的基因里预后最好的
        best = None
        for g in m['genes']:
            if g in genes:
                if best is None or (genes[g]['n_circRNA'] and g):
                    best = g
        gi = genes.get(best, {}) if best else {}
        imm_max = max((abs(x[1]) for g2 in m['genes'] for x in imm.get(g2, [])), default=0)
        tiers = sorted(set(r['tier'] for r in rows if r['miRNA'] == name))
        w.writerow(['miRNA', name, m['deg'], '%.1f' % avg_e,
                    ';'.join(sorted(m['genes'])), '', '',
                    gi.get('HR_TCGA', ''), gi.get('P_cox_TCGA', ''), '',
                    '%.3f' % imm_max, 'YES' if best in val else '',
                    ';'.join(tiers)])
    for name, c in sorted(circ_info.items(), key=lambda kv: -len(kv[1]['mirnas'])):
        if len(c['mirnas']) < 2:
            continue
        g = c['gene']
        gi = genes.get(g, {})
        tiers = sorted(set(r['tier'] for r in rows if r['circ'] == name))
        w.writerow(['circRNA', name, len(c['mirnas']), '',
                    g, '%.4f' % c['log2fc'], c['reg'],
                    gi.get('HR_TCGA', ''), gi.get('P_cox_TCGA', ''),
                    gi.get('cross_cohort_concordant', ''),
                    gi.get('immune_max_r', ''), 'YES' if g in val else '',
                    ';'.join(tiers)])
print('\n[SAVED] %s/trinet_nodes.tsv' % OUT)

# ---------- 摘要 ----------
print('\n===== 三元网络摘要 =====')
for k in ('high', 'medium', 'low'):
    sub = [r for r in rows if r['tier'] == k]
    hg = [g for g in set(x['gene'] for x in sub) if g in val]
    print('%s 档: %d 边 | miRNA %d | circRNA %d | 基因 %d | 验证基因命中 %d %s'
          % (k, len(sub), len(set(x['miRNA'] for x in sub)),
             len(set(x['circ'] for x in sub)), len(set(x['gene'] for x in sub)),
             len(hg), hg))

# 保存摘要
with open(f'{OUT}/trinet_summary.tsv', 'w', encoding='utf-8', newline='') as fo:
    w = csv.writer(fo, delimiter='\t')
    w.writerow(['tier', 'edges', 'n_miRNA', 'n_circRNA', 'n_gene', 'validated_genes_hit'])
    for k in ('high', 'medium', 'low'):
        sub = [r for r in rows if r['tier'] == k]
        hg = [g for g in set(x['gene'] for x in sub) if g in val]
        w.writerow([k, len(sub), len(set(x['miRNA'] for x in sub)),
                    len(set(x['circ'] for x in sub)), len(set(x['gene'] for x in sub)),
                    ';'.join(hg)])
print('[SAVED] %s/trinet_summary.tsv' % OUT)