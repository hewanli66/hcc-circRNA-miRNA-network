# -*- coding: utf-8 -*-
"""
NT5DC2 的 miRNA 可信度间接验证（不依赖 TargetScan/ENCORI）
思路：NT5DC2 因 circBank 坐标无法精确匹配，其 circRNA-miRNA 边只能靠
     gene_fallback 兜底 → 可信度低。
     但若某个 miRNA 同时连接了**多个 exact 坐标匹配的基因**，
     说明该 miRNA 的结合预测在多基因上被 circBank 一致认可，
     可作为 NT5DC2 侧 miRNA 的「间接可信证据」。
输出: results/NT5DC2_mirna_indirect_support.tsv
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

TARGET = 'NT5DC2'
mir_map = defaultdict(list)
for r in rows:
    mir_map[r['miRNA']].append(r)


def fe(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


out = []
for r in rows:
    if r['host_gene'] != TARGET:
        continue
    m = r['miRNA']
    # 同一 miRNA 在 exact 匹配基因上的其他边
    others = [x for x in mir_map[m]
              if x['match_type'] == 'exact' and x['host_gene'] != TARGET]
    og = sorted(set(x['host_gene'] for x in others))
    # 最强结合能
    ens = [fe(x['energy_kcal_mol']) for x in others if fe(x['energy_kcal_mol']) is not None]
    best_e = min(ens) if ens else None
    # 间接支持等级
    if len(og) >= 5:
        tier_ind = 'strong'
    elif len(og) >= 2:
        tier_ind = 'moderate'
    elif len(og) == 1:
        tier_ind = 'weak'
    else:
        tier_ind = 'none'
    # 那些 exact 基因里有几个是双队列验证基因
    val_others = [g for g in og if g in val]
    out.append(dict(miRNA=m, energy=fe(r['energy_kcal_mol']), score=r['score'],
                    n_exact_genes=len(og), exact_genes=';'.join(og[:12]),
                    n_validated_others=len(val_others), validated_others=';'.join(val_others),
                    best_other_energy=best_e, indirect_tier=tier_ind))

out.sort(key=lambda x: (-x['n_exact_genes'], -(x['best_other_energy'] or 0)))

with open(f'{OUT}/{TARGET}_mirna_indirect_support.tsv', 'w', encoding='utf-8', newline='') as fo:
    w = csv.writer(fo, delimiter='\t')
    w.writerow(['miRNA', 'energy_NT5DC2', 'best_other_energy', 'n_exact_genes',
                'exact_genes_supporting', 'n_validated_genes_among_them',
                'validated_genes', 'indirect_support_tier'])
    for r in out:
        w.writerow([r['miRNA'], r['energy'] if r['energy'] is not None else '',
                    r['best_other_energy'] if r['best_other_energy'] is not None else '',
                    r['n_exact_genes'], r['exact_genes'], r['n_validated_others'],
                    r['validated_others'], r['indirect_tier']])
print('[SAVED] %s/%s_mirna_indirect_support.tsv' % (OUT, TARGET))

print('\n%s 的 miRNA 间接支持分布:' % TARGET)
for k in ('strong', 'moderate', 'weak', 'none'):
    sub = [r for r in out if r['indirect_tier'] == k]
    print('  %-9s %4d 个' % (k, len(sub)))

print('\n=== 强支持 + 结合最强的 12 个（优先机制候选）===')
print('%-20s %-9s %-10s %-6s %s' % ('miRNA', 'E(NT5DC2)', 'E(其他最佳)', 'exact基因', '含验证基因'))
for r in out:
    if r['indirect_tier'] not in ('strong', 'moderate'):
        continue
    print('%-20s %-9s %-10s %-6d %s' % (
        r['miRNA'],
        '%.1f' % r['energy'] if r['energy'] is not None else '-',
        '%.1f' % r['best_other_energy'] if r['best_other_energy'] is not None else '-',
        r['n_exact_genes'], r['validated_others'] or '-'))
    if sum(1 for x in out[:20] if x['indirect_tier'] in ('strong', 'moderate')) > 12:
        pass
    if r['indirect_tier'] == 'moderate':
        break

#免疫关联：这些 miRNA 连接的 exact 基因是否与免疫浸润相关
print('\n=== 免疫浸润相关的支撑证据 ===')
cnt = 0
for r in out:
    if r['indirect_tier'] not in ('strong', 'moderate'):
        continue
    im = []
    for g in r['exact_genes'].split(';'):
        if g in imm:
            im.append('%s(r=%.2f)' % (g, max(abs(x[1]) for x in imm[g])))
    if im:
        print('  %-20s ← %s' % (r['miRNA'], ', '.join(im[:3])))
        cnt += 1
    if cnt >= 8:
        break