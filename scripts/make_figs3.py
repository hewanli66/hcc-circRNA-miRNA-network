# -*- coding: utf-8 -*-
"""
三元网络图：circRNA - miRNA - mRNA
分层展示 high/medium/low，画 hub miRNA 与核心基因
输出: figures/Fig6_trinetwork.svg/png, figures/Fig7_qc_tiers.svg/png
"""
import csv, os, math, json
from collections import defaultdict, Counter

OUT = 'results'
FIG = 'figures'
os.makedirs(FIG, exist_ok=True)


def esc(s):
    return str(s).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def txt(x, y, s, size=11, anchor='start', fill='#222', weight='normal'):
    return ('<text x="%.1f" y="%.1f" font-size="%d" text-anchor="%s" fill="%s" '
            'font-weight="%s" font-family="Helvetica,Arial,sans-serif">%s</text>'
            % (x, y, size, anchor, fill, weight, esc(s)))


rows = list(csv.DictReader(open(f'{OUT}/network_edges_qc.tsv', encoding='utf-8'), delimiter='\t'))
genes = {r['gene']: r for r in csv.DictReader(open(f'{OUT}/network_genes.tsv', encoding='utf-8'),
                                             delimiter='\t')}
val = set(json.load(open(f'{OUT}/validated_genes.json'))) if os.path.exists(
    f'{OUT}/validated_genes.json') else set()

print('边 %d' % len(rows))

TIERS = [('high', '#c0392b', 'High tier (exact + energy <= -20)'),
         ('medium', '#e08214', 'Medium tier (exact/loose + energy <= -20)'),
         ('low', '#8a93a0', 'Low tier (fallback or weak binding)')]


def fe(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


# ==================== Fig6: 三元网络（high+medium）===================
use = [r for r in rows if r['tier'] in ('high', 'medium')]
print('绘制用边: %d（high+medium）' % len(use))

# hub miRNA（按 degree）
mir = defaultdict(lambda: {'circs': set(), 'genes': set(), 'es': []})
for r in use:
    m = mir[r['miRNA']]
    m['circs'].add(r['circRNA'])
    m['genes'].add(r['host_gene'])
    e = fe(r['energy_kcal_mol'])
    if e is not None:
        m['es'].append(e)

HUB_N = 16
hubs = [k for k, v in sorted(mir.items(), key=lambda kv: -len(kv[1]['circs']))[:HUB_N]]
hubset = set(hubs)
# 只画与 hub 相连的 circRNA / gene
hub_circs, hub_genes = set(), set()
for h in hubs:
    hub_circs |= mir[h]['circs']
    hub_genes |= mir[h]['genes']
print('  hub miRNA %d，涉及 circRNA %d，基因 %d'
      % (len(hubs), len(hub_circs), len(hub_genes)))

W, H = 1300, 900
cx, cy = W / 2, H / 2 + 20
R1, R2 = 355, 215      # circRNA 环/ 基因环
svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">' % (W, H, W, H),
       '<rect width="%d" height="%d" fill="#fff"/>' % (W, H),
       txt(24, 30, 'Figure 6  circRNA-miRNA-mRNA tri-network (high and medium confidence tiers)', 17, weight='bold'),
       txt(24, 50, 'Middle: %d hub miRNAs   Outer left: %d circRNAs   Outer right: %d host genes' % (len(hubs), len(hub_circs), len(hub_genes)), 11, fill='#666')]

# 节点角度
circs_sorted = sorted(hub_circs)
genes_sorted = sorted(hub_genes)
nc, ng = len(circs_sorted), len(genes_sorted)

circ_ang = {}
for i, c in enumerate(circs_sorted):
    circ_ang[c] = 2 * math.pi * i / max(nc, 1) - math.pi / 2
gene_ang = {}
for i, g in enumerate(genes_sorted):
    gene_ang[g] = 2 * math.pi * i / max(ng, 1) - math.pi / 2

# 边：miRNA(中心) -> circRNA / gene
TCOL = {'high': '#c0392b', 'medium': '#e08214'}
for r in use:
    h = r['miRNA']
    if h not in hubset:
        continue
    col = TCOL[r['tier']]
    op = 0.10 if r['tier'] == 'medium' else 0.16
    mx, my = cx, cy
    a = circ_ang.get(r['circRNA'])
    if a is None:
        continue
    px = cx + R1 * math.cos(a); py = cy + R1 * math.sin(a)
    svg.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-opacity="%.2f" stroke-width="0.8"/>'
               % (mx, my, px, py, col, op))
    # gene 点
    a2 = gene_ang.get(r['host_gene'])
    if a2 is not None:
        gx = cx + R2 * math.cos(a2); gy = cy + R2 * math.sin(a2)
        svg.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-opacity="%.2f" stroke-width="0.8"/>'
                   % (mx, my, gx, gy, col, op))

# circRNA 节点
for c in circs_sorted:
    a = circ_ang[c]
    x = cx + R1 * math.cos(a); y = cy + R1 * math.sin(a)
    tier = 'high' if any(r['tier'] == 'high' for r in use if r['circRNA'] == c) else 'medium'
    col = TCOL[tier]
    rr = 4 if tier == 'medium' else 5.5
    svg.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s" fill-opacity="0.75"/>' % (x, y, rr, col))

# gene 节点（预后显著者加大加粗）
for g in genes_sorted:
    a = gene_ang[g]
    x = cx + R2 * math.cos(a); y = cy + R2 * math.sin(a)
    gi = genes.get(g, {})
    isv = g in val
    col = '#8b1a1a' if isv else ('#4a6fa5' if gi.get('P_cox_TCGA') and
            gi['P_cox_TCGA'] not in ('', 'nan') and float(gi['P_cox_TCGA']) < 0.05 else '#aab0b8')
    rr = 8 if isv else 5.5
    svg.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s" stroke="#fff" stroke-width="1.3"/>'
               % (x, y, rr, col))

# miRNA 中心节点（统一放在中心，正对小圆环）
for h in hubs:
    n = len(mir[h]['circs'])
    rr = 5 + min(13, n * 0.55)
    avg_e = sum(mir[h]['es']) / len(mir[h]['es']) if mir[h]['es'] else 0
    shade = max(0, min(1, (-avg_e - 18) / 16))
    col = '#%02x%02x%02x' % (int(120 + 100 * shade), int(30 + 60 * shade), int(90 + 60 * shade))
    svg.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s" stroke="#fff" stroke-width="1.5"/>'
               % (cx, cy, rr, col))

# 标签：miRNA 沿垂直方向排在中心右侧，按 y 排序后拉开最小间距
hl = []
for h in hubs:
    n = len(mir[h]['circs'])
    rr = 5 + min(13, n * 0.55)
    a = circ_ang.get(sorted(mir[h]['circs'])[0], 0)
    x = cx + (55 + rr + 6) * math.cos(a)
    y = cy + (55 + rr + 6) * math.sin(a)
    hl.append([h.replace('hsa-', ''), x, y])
# 中心区域：强制放到中右侧的标签带，按 y 拉开
band_x = cx + 78
band_cy = []
for name, x, y in hl:
    band_cy.append([name, min(max(x, band_x), cx + R1 - 70), y])
band_cy.sort(key=lambda t: t[2])
GAP = 13
for i in range(1, len(band_cy)):
    if band_cy[i][2] - band_cy[i - 1][2] < GAP:
        band_cy[i][2] = band_cy[i - 1][2] + GAP
# 垂直居中整体
if band_cy:
    over = band_cy[-1][2] - (cy + R1 - 20)
    if over > 0:
        for t in band_cy:
            t[2] -= over
    under = (cy - R1 + 30) - band_cy[0][2]
    if under > 0:
        for t in band_cy:
            t[2] += under
for name, x, ty in band_cy:
    svg.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#d8c0cc" stroke-width="0.7"/>'
               % (cx + 18, cy, x - 3, ty - 3))
    svg.append(txt(x, ty, name, 10, 'start', '#5a1a3a', 'bold'))

# 基因标签（只标验证基因 + 预后显著者）
gl = []
for g in genes_sorted:
    a = gene_ang[g]
    gi = genes.get(g, {})
    sig = gi.get('P_cox_TCGA') not in ('', None, 'nan') and gi.get('P_cox_TCGA') is not None \
        and float(gi['P_cox_TCGA']) < 0.05
    if g not in val and not sig:
        continue
    x = cx + (R2 + 14) * math.cos(a); y = cy + (R2 + 14) * math.sin(a)
    gl.append([g, x, y, a])
gl.sort(key=lambda t: (t[2]))
last_y = -999
for name, x, y, a in gl:
    ty = max(y, last_y + 12)
    last_y = ty
    anc = 'start' if math.cos(a) > 0.2 else ('end' if math.cos(a) < -0.2 else 'middle')
    col = '#8b1a1a' if name in val else '#2a4a7a'
    svg.append(txt(x, ty, name, 10, anc, col, 'bold'))

# 图例
lgx, lgy = 24, H - 110
svg.append(txt(lgx, lgy - 12, 'Legend', 12, weight='bold'))
for i, (k, col, desc) in enumerate(TIERS):
    svg.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" stroke-width="2"/>'
               % (lgx, lgy + i * 16, lgx + 22, lgy + i * 16, col))
    svg.append(txt(lgx + 28, lgy + i * 16 + 4, desc, 10, fill='#444'))
svg.append('<circle cx="%d" cy="%d" r="7" fill="#8b1a1a"/>' % (lgx + 8, lgy + 60))
svg.append(txt(lgx + 28, lgy + 64, 'gene validated in two cohorts', 10, fill='#444'))
svg.append('<circle cx="%d" cy="%d" r="5" fill="#4a6fa5"/>' % (lgx + 8, lgy + 78))
svg.append(txt(lgx + 28, lgy + 82, 'gene P_cox < 0.05 (TCGA)', 10, fill='#444'))
svg.append('<circle cx="%d" cy="%d" r="5" fill="#aab0b8"/>' % (lgx + 8, lgy + 96))
svg.append(txt(lgx + 28, lgy + 100, 'gene not significant', 10, fill='#444'))
svg.append('</svg>')
open(f'{FIG}/Fig6_trinetwork.svg', 'w', encoding='utf-8').write('\n'.join(svg))
print('[SAVED] Fig6_trinetwork.svg')


# ==================== Fig7: 质量分层对比 ====================
W2, H2 = 900, 460
svg2 = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">' % (W2, H2, W2, H2),
        '<rect width="%d" height="%d" fill="#fff"/>' % (W2, H2),
        txt(24, 30, 'Figure 7  Quality stratification of circRNA-miRNA edges', 17, weight='bold'),
        txt(24, 50, 'Tier assignment: match confidence (exact > loose > gene fallback) + binding free energy + circBank score', 11, fill='#666')]

# 三个面板：边数 / miRNA数 / 验证基因命中
panels = [('Edges', 'n_edges'), ('Distinct miRNAs', 'n_miRNA'), ('Validated genes covered', 'n_val')]
bx, bw, bh = 90, 210, 220
by = H2 - 150
mx = 0
for name, _ in panels:
    for k, _, _ in TIERS:
        sub = [r for r in rows if r['tier'] == k]
        if name == 'Edges':
            v = len(sub)
        elif name == 'Distinct miRNAs':
            v = len(set(r['miRNA'] for r in sub))
        else:
            v = len(set(r['host_gene'] for r in sub) & val)
        mx = max(mx, v)

for pi, (name, _) in enumerate(panels):
    x0 = bx + pi * (bw + 60)
    svg2.append('<rect x="%.1f" y="%d" width="%d" height="%d" fill="#fafafa" stroke="#ddd"/>'
                % (x0, by - bh, bw, bh))
    for gi, (k, col, desc) in enumerate(TIERS):
        sub = [r for r in rows if r['tier'] == k]
        if name == 'Edges':
            v = len(sub)
        elif name == 'Distinct miRNAs':
            v = len(set(r['miRNA'] for r in sub))
        else:
            v = len(set(r['host_gene'] for r in sub) & val)
        hh = v / mx * (bh - 30)
        xx = x0 + 18 + gi * 62
        svg2.append('<rect x="%.1f" y="%.1f" width="46" height="%.1f" fill="%s" fill-opacity="0.85"/>'
                    % (xx, by - hh, hh, col))
        svg2.append(txt(xx + 23, by - hh - 6, str(v), 11, 'middle', '#222', 'bold'))
        svg2.append(txt(xx + 23, by + 14, k, 10, 'middle', '#666'))
    svg2.append(txt(x0 + bw / 2, by - bh - 12, name, 12, 'middle', '#222', 'bold'))
svg2.append(txt(24, H2 - 30, 'Red = high, orange = medium, grey = low. High tier requires exact coordinate match and binding energy <= -20 kcal/mol.', 10, fill='#666'))
svg2.append('</svg>')
open(f'{FIG}/Fig7_qc_tiers.svg', 'w', encoding='utf-8').write('\n'.join(svg2))
print('[SAVED] Fig7_qc_tiers.svg')
print('\n完成')