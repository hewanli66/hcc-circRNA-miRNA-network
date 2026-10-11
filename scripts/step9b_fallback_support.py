# -*- coding: utf-8 -*-
"""
第 9 步（扩展）：把所有 gene_fallback 匹配的基因做「间接支持」分析
背景：这些基因的 circRNA-miRNA 边靠「同基因同染色体兜底」匹配，可信度低。
     但若其连接的 miRNA 同时被多个 exact 匹配基因连接，则构成间接可信证据。
输出: results/fallback_genes_indirect_support.tsv（每基因汇总）
      results/fallback_gene_mirna_support.tsv（每条边明细）
"""
import csv, json, os
from collections import defaultdict, Counter

OUT = 'results'
rows = list(csv.DictReader(open(f'{OUT}/network_edges_qc.tsv', encoding='utf-8'), delimiter='\t'))
genes = {r['gene']: r for r in csv.DictReader(open(f'{OUT}/network_genes.tsv', encoding='utf-8'),
                                             delimiter='\t')}
val = set(json.load(open(f'{OUT}/validated_genes.json'))) if os.path.exists(
    f'{OUT}/validated_genes.json') else set()
imm = defaultdict(list)
if os.path.exists(f'{OUT}/immune_correlation_2queues.tsv'):
    for r in csv.DictReader(open(f'{OUT}/immune_correlation_2queues.tsv', encoding='utf-8'),
                            delimiter='\t'):
        if r['direction_consistent'] == 'YES':
            try:
                imm[r['gene']].append((r['cell_type'], float(r['r_TCGA'])))
            except (ValueError, KeyError):
                pass

# 所有 fallback 基因
fb_genes = sorted(set(r['host_gene'] for r in rows if r['match_type'] == 'gene_fallback'))
print('gene_fallback 基因: %d 个' % len(fb_genes))
print('  %s\n' % ', '.join(fb_genes))


def fe(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


mir_map = defaultdict(list)
for r in rows:
    mir_map[r['miRNA']].append(r)

# ---------- 每条 fallback 边的间接支持 ----------
detail = []
for r in rows:
    if r['match_type'] != 'gene_fallback':
        continue
    m = r['miRNA']
    others = [x for x in mir_map[m] if x['match_type'] == 'exact' and x['host_gene'] != r['host_gene']]
    og = sorted(set(x['host_gene'] for x in others))
    ens = [fe(x['energy_kcal_mol']) for x in others if fe(x['energy_kcal_mol']) is not None]
    best_e = min(ens) if ens else None
    if len(og) >= 5:
        t = 'strong'
    elif len(og) >= 2:
        t = 'moderate'
    elif len(og) == 1:
        t = 'weak'
    else:
        t = 'none'
    vo = [g for g in og if g in val]
    # miRNA 是否连接免疫相关基因
    im = []
    for g in og:
        if g in imm:
            im.append('%s(r=%.2f)' % (g, max(abs(x[1]) for x in imm[g])))
    detail.append(dict(circ=r['circRNA'], cb=r['circBankID'], gene=r['host_gene'],
                       miRNA=m, energy=fe(r['energy_kcal_mol']), score=r['score'],
                       n_og=len(og), og=';'.join(og[:12]), n_vo=len(vo),
                       vo=';'.join(vo), best_e=best_e, tier=t,
                       immune=';'.join(sorted(set(im))[:4])))

with open(f'{OUT}/fallback_gene_mirna_support.tsv', 'w', encoding='utf-8', newline='') as fo:
    w = csv.writer(fo, delimiter='\t')
    w.writerow(['circRNA', 'circBankID', 'fallback_gene', 'miRNA', 'energy',
                'score', 'n_exact_genes_supporting', 'exact_genes', 'n_validated_among_them',
                'validated_genes', 'best_other_energy', 'indirect_support_tier',
                'immune_linked_genes'])
    for r in sorted(detail, key=lambda x: (-x['n_og'], -(x['best_e'] or 0))):
        w.writerow([r['circ'], r['cb'], r['gene'], r['miRNA'],
                    r['energy'] if r['energy'] is not None else '',
                    r['score'], r['n_og'], r['og'], r['n_vo'], r['vo'],
                    r['best_e'] if r['best_e'] is not None else '', r['tier'], r['immune']])
print('[SAVED] %s/fallback_gene_mirna_support.tsv' % OUT)

# ---------- 每个基因汇总 ----------
summ = []
for g in fb_genes:
    sub = [d for d in detail if d['gene'] == g]
    if not sub:
        continue
    tc = Counter(d['tier'] for d in sub)
    gi = genes.get(g, {})
    # 该基因的最佳 miRNA（strong 优先，再看能量）
    strong = [d for d in sub if d['tier'] in ('strong', 'moderate')]
    best = sorted(strong or sub, key=lambda d: (-d['n_og'], -(d['energy'] or 0)))[:5]
    # 免疫交叉：多少 miRNA 连到了免疫相关基因
    ims = [d for d in sub if d['immune']]
    summ.append(dict(
        gene=g, n_edges=len(sub),
        n_strong=tc['strong'], n_moderate=tc['moderate'],
        n_weak=tc['weak'], n_none=tc['none'],
        pct_strong=round((tc['strong'] + tc['moderate']) / len(sub) * 100, 1),
        n_mirna=len(set(d['miRNA'] for d in sub)),
        HR_TCGA=gi.get('HR_TCGA', ''), P_cox=gi.get('P_cox_TCGA', ''),
        HR_14520=gi.get('HR_14520', ''), P_cox_14520=gi.get('P_cox_14520', ''),
        cross=gi.get('cross_cohort_concordant', ''),
        immune_max=gi.get('immune_max_r', ''), immune_cells=gi.get('immune_cell_types', ''),
        validated='YES' if g in val else '',
        n_immune_linked=len(ims),
        top_mirna=';'.join(d['miRNA'] for d in best),
        top_mirna_support=';'.join('%d个exact基因' % d['n_og'] for d in best)))

summ.sort(key=lambda x: (-x['n_strong'] - x['n_moderate'], -x['n_edges']))
with open(f'{OUT}/fallback_genes_indirect_support.tsv', 'w', encoding='utf-8', newline='') as fo:
    w = csv.writer(fo, delimiter='\t')
    w.writerow(['gene', 'n_edges', 'n_mirna', 'n_strong', 'n_moderate', 'n_weak', 'n_none',
                'pct_supported', 'HR_TCGA', 'P_cox_TCGA', 'HR_14520', 'P_cox_14520',
                'cross_cohort', 'immune_max_r', 'immune_cell_types', 'validated',
                'n_mirna_immune_linked', 'top_miRNA', 'top_miRNA_support'])
    for r in summ:
        w.writerow([r['gene'], r['n_edges'], r['n_mirna'], r['n_strong'], r['n_moderate'],
                    r['n_weak'], r['n_none'], r['pct_strong'], r['HR_TCGA'], r['P_cox'],
                    r['HR_14520'], r['P_cox_14520'], r['cross'], r['immune_max'],
                    r['immune_cells'], r['validated'], r['n_immune_linked'],
                    r['top_mirna'], r['top_mirna_support']])
print('[SAVED] %s/fallback_genes_indirect_support.tsv' % OUT)

print('\n===== 支持度最高的 fallback 基因（可作为机制候选）=====')
print('%-14s %-6s %-6s %-7s %-6s %-8s %s' % ('gene', '边', 'miRNA', 'strong', '支持率', '验证', 'top miRNA'))
for r in summ[:14]:
    print('%-14s %-6d %-6d %-7d %-6.1f%% %-8s %s'
          % (r['gene'], r['n_edges'], r['n_mirna'], r['n_strong'], r['pct_strong'],
             r['validated'] or '-', r['top_mirna'][:44]))

print('\n===== fallback 基因的双队列验证情况 =====')
for r in summ:
    if r['validated'] == 'YES':
        print('  %-12s HR(TCGA)=%-7s P=%-9s HR(14520)=%-7s P=%-9s 支持率 %.1f%%'
              % (r['gene'], r['HR_TCGA'], r['P_cox'], r['HR_14520'], r['P_cox_14520'],
                 r['pct_strong']))

#免疫交叉汇总
print('\n===== fallback 基因的免疫交叉证据（仅 |r|>=0.15）=====')
found = 0
for r in summ:
    try:
        mx = float(r['immune_max'])
    except (TypeError, ValueError):
        mx = 0.0
    if mx >= 0.15 and r['immune_cells']:
        print('  %-14s %s (r=%.3f) | %d 个 miRNA 有免疫交叉'
              % (r['gene'], r['immune_cells'][:44], mx, r['n_immune_linked']))
        found += 1
if found == 0:
    print('  （无：所有 fallback 基因的免疫关联 |r| 均 < 0.15，'
          '不构成免疫机制证据）')
else:
    print('\n  说明：全局仅 7 个基因有双队列一致的免疫关联，'
          '其中仅 CD5 相关性较强（r=0.836–0.877）')

print('\n===== 全局免疫关联基因一览（共 %d 个）=====' % len(imm))
for g, v in sorted(imm.items(), key=lambda kv: -max(abs(x[1]) for x in kv[1])):
    top = max(v, key=lambda x: abs(x[1]))
    print('  %-10s %-18s r=%+.3f' % (g, top[0], top[1]))