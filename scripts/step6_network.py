# -*- coding: utf-8 -*-
"""
第6 步：构建 circRNA–mRNA 两层网络骨架 + 生成图表
1) 网络节点/边整合（circRNA → 基因 → 预后/免疫属性）
2) 核心模块识别
3) 生成投稿级图表（SVG）
输出: results/network_nodes.tsv, network_edges.tsv, network_summary.tsv
      figures/Fig1-Fig5
"""
import csv, os, json, math, gzip, re
from collections import defaultdict

RAW = 'C:/Users/hewan/Documents/kimi/tasks/2026-10-08/11-38-08-f56eebff/rawdata'
OUT = 'results'
FIG = 'figures'
os.makedirs(FIG, exist_ok=True)

# ============ 1) 读候选 circRNA ============
cands = []
with open(f'{OUT}/candidate_circRNA.tsv', encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        cands.append(dict(circ=r['circRNA'], gene=r['gene'],
                          normal=float(r['mean_normal_RPM']), tumor=float(r['mean_tumor_RPM']),
                          log2fc=float(r['log2FC']), P=float(r['P']), padj=float(r['P_adj']),
                          reg=r['regulation']))
print('[1/5] 候选 circRNA: %d（基因 %d）' % (len(cands), len(set(c['gene'] for c in cands))))

# ============ 2) 读 TCGA 预后 ============
tcga = {}
with open(f'{OUT}/tcga_candidate_gene_stats.tsv', encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        tcga[r['gene']] = r
print('[2/5] TCGA 统计: %d 基因' % len(tcga))

# ============ 3) 读 GSE14520 验证 ============
g14520 = {}
p14520 = f'{OUT}/gse14520_validation.tsv'
if os.path.exists(p14520):
    with open(p14520, encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            g14520[r['gene']] = r
print('[3/5] GSE14520 验证: %d 基因' % len(g14520))

# ============ 4) 读免疫关联 ============
imm = []
p2 = f'{OUT}/immune_correlation_2queues.tsv'
if os.path.exists(p2):
    with open(p2, encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            imm.append(r)
print('[4/5] 免疫关联记录: %d' % len(imm))

# ============ 5) 构建网络 ============
# 节点：circRNA（外圈）+ 基因（内圈）
# 属性：预后显著性、跨队列一致性、免疫关联强度
nodes = []
edges = []

# 基因 -> 免疫关联汇总
imm_by_gene = defaultdict(list)
for r in imm:
    if r['direction_consistent'] == 'YES':
        try:
            imm_by_gene[r['gene']].append((r['cell_type'], float(r['r_TCGA'])))
        except (ValueError, KeyError):
            pass

for c in cands:
    g = c['gene']
    t = tcga.get(g, {})
    v = g14520.get(g, {})
    # 跨队列一致性
    hr_t = float(t['HR_per_SD']) if t.get('HR_per_SD') not in (None, '', 'nan') else None
    hr_v = float(v['HR_14520']) if v.get('HR_14520') not in (None, '', 'nan') else None
    if hr_t and hr_v:
        cons = (hr_t > 1) == (hr_v > 1)
    else:
        cons = None
    imms = imm_by_gene.get(g, [])
    max_imm = max((abs(x[1]) for x in imms), default=0.0)

    nodes.append(dict(
        node_id=c['circ'], node_type='circRNA', gene=g,
        log2fc=c['log2fc'], P=c['P'], reg=c['reg'],
        normal=c['normal'], tumor=c['tumor'],
        hr_tcga=hr_t, p_cox_tcga=float(t['P_cox']) if t.get('P_cox') not in (None, '', 'nan') else None,
        p_lr_tcga=float(t['P_logrank']) if t.get('P_logrank') not in (None, '', 'nan') else None,
        hr_g14520=hr_v, p_cox_g14520=float(v['P_cox_14520']) if v.get('P_cox_14520') not in (None, '', 'nan') else None,
        concordant=('YES' if cons else ('no' if cons is not None else '')),
        immune_max=round(max_imm, 4), immune_n=len(imms),
        immune_cells=';'.join(x[0] for x in sorted(imms, key=lambda y: -abs(y[1]))[:3]),
    ))
    edges.append(dict(source=c['circ'], target=g, edge_type='circRNA_encodes_gene',
                      weight=abs(c['log2fc'])))

# 汇总基因层
genes = {}
for n in nodes:
    g = n['gene']
    if g not in genes:
        genes[g] = dict(gene=g, n_circ=0, log2fc_sum=0.0, min_P=1.0,
                        hr_tcga=n['hr_tcga'], p_cox_tcga=n['p_cox_tcga'],
                        p_lr_tcga=n['p_lr_tcga'], hr_g14520=n['hr_g14520'],
                        p_cox_g14520=n['p_cox_g14520'], concordant=n['concordant'],
                        immune_max=n['immune_max'], immune_n=0, immune_cells='')
    G = genes[g]
    G['n_circ'] += 1
    G['log2fc_sum'] += n['log2fc']
    G['min_P'] = min(G['min_P'], n['P'])
    if n['immune_max'] > G['immune_max']:
        G['immune_max'] = n['immune_max']
        G['immune_cells'] = n['immune_cells']
    G['immune_n'] += n['immune_n']

with open(f'{OUT}/network_nodes.tsv', 'w', encoding='utf-8', newline='') as fo:
    w = csv.writer(fo, delimiter='\t')
    w.writerow(['node_id', 'node_type', 'gene', 'log2FC', 'P', 'regulation',
                'HR_TCGA', 'P_cox_TCGA', 'P_logrank_TCGA',
                'HR_14520', 'P_cox_14520', 'cross_cohort_concordant',
                'immune_max_r', 'immune_cell_types'])
    for n in nodes:
        w.writerow([n['node_id'], n['node_type'], n['gene'],
                    '%.4f' % n['log2fc'], '%.6g' % n['P'], n['reg'],
                    n['hr_tcga'] if n['hr_tcga'] is not None else '',
                    '%.6g' % n['p_cox_tcga'] if n['p_cox_tcga'] is not None else '',
                    '%.6g' % n['p_lr_tcga'] if n['p_lr_tcga'] is not None else '',
                    n['hr_g14520'] if n['hr_g14520'] is not None else '',
                    '%.6g' % n['p_cox_g14520'] if n['p_cox_g14520'] is not None else '',
                    n['concordant'], '%.4f' % n['immune_max'], n['immune_cells']])

with open(f'{OUT}/network_genes.tsv', 'w', encoding='utf-8', newline='') as fo:
    w = csv.writer(fo, delimiter='\t')
    w.writerow(['gene', 'n_circRNA', 'sum_log2FC', 'min_P', 'HR_TCGA', 'P_cox_TCGA',
                'P_logrank_TCGA', 'HR_14520', 'P_cox_14520', 'cross_cohort_concordant',
                'immune_max_r', 'immune_cell_types'])
    for g, G in sorted(genes.items(), key=lambda x: -x[1]['n_circ']):
        w.writerow([g, G['n_circ'], '%.3f' % G['log2fc_sum'], '%.6g' % G['min_P'],
                    G['hr_tcga'] if G['hr_tcga'] is not None else '',
                    '%.6g' % G['p_cox_tcga'] if G['p_cox_tcga'] is not None else '',
                    '%.6g' % G['p_lr_tcga'] if G['p_lr_tcga'] is not None else '',
                    G['hr_g14520'] if G['hr_g14520'] is not None else '',
                    '%.6g' % G['p_cox_g14520'] if G['p_cox_g14520'] is not None else '',
                    G['concordant'], '%.4f' % G['immune_max'], G['immune_cells']])

with open(f'{OUT}/network_edges.tsv', 'w', encoding='utf-8', newline='') as fo:
    w = csv.writer(fo, delimiter='\t')
    w.writerow(['source_circRNA', 'target_gene', 'edge_type', 'weight_abs_log2FC'])
    for e in edges:
        w.writerow([e['source'], e['target'], e['edge_type'], '%.4f' % e['weight']])

print('[SAVED] network_nodes.tsv / network_genes.tsv / network_edges.tsv')

# ============ 网络摘要 ============
multi = [g for g, G in genes.items() if G['n_circ'] > 1]
both_sig = [g for g, G in genes.items()
            if G['p_cox_tcga'] is not None and G['p_lr_tcga'] is not None
            and G['p_cox_tcga'] < 0.05 and G['p_lr_tcga'] < 0.05]
validated = [g for g, G in genes.items()
             if G['concordant'] == 'YES' and G['p_cox_g14520'] is not None
             and G['p_cox_g14520'] < 0.05]
print('\n===== 网络摘要 =====')
print('circRNA 节点: %d' % len(nodes))
print('基因节点: %d' % len(genes))
print('关联多个 circRNA 的基因: %d %s' % (len(multi), sorted(multi)))
print('TCGA 双检验显著基因: %d' % len(both_sig))
print('跨队列验证一致且显著: %d %s' % (len(validated), validated))
print('有免疫关联的基因: %d' % sum(1 for G in genes.values() if G['immune_max'] > 0))

with open(f'{OUT}/network_summary.tsv', 'w', encoding='utf-8', newline='') as fo:
    w = csv.writer(fo, delimiter='\t')
    w.writerow(['metric', 'value'])
    for k, v in [('circRNA_nodes', len(nodes)), ('gene_nodes', len(genes)),
                 ('multi_circRNA_genes', ';'.join(sorted(multi))),
                 ('TCGA_dual_significant', len(both_sig)),
                 ('cross_cohort_validated', ';'.join(validated)),
                 ('genes_with_immune_link', sum(1 for G in genes.values() if G['immune_max'] > 0))]:
        w.writerow([k, v])
print('[SAVED] network_summary.tsv')