# -*- coding: utf-8 -*-
"""
生成投稿级图表（纯 SVG，零依赖）
Fig1  circRNA 差异表达火山图
Fig2  网络结构图（circRNA 外圈 + 基因内圈）
Fig3  跨队列预后一致性散点
Fig4  免疫关联热图
Fig5  核心基因 KM 分面示意
"""
from collections import defaultdict
import csv, os, math

OUT = 'results'
FIG = 'figures'
os.makedirs(FIG, exist_ok=True)

# ---------- 工具 ----------
def esc(s):
    return (str(s).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;'))


def txt(x, y, s, size=12, anchor='start', fill='#222', weight='normal', rot=None):
    tr = ' transform="rotate(%s %.1f %.1f)"' % (rot, x, y) if rot else ''
    return ('<text x="%.1f" y="%.1f" font-size="%d" text-anchor="%s" fill="%s" '
            'font-weight="%s" font-family="Helvetica,Arial,sans-serif"%s>%s</text>'
            % (x, y, size, anchor, fill, weight, tr, esc(s)))


def load_nodes():
    nodes = []
    with open(f'{OUT}/network_nodes.tsv', encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            nodes.append(r)
    return nodes


def load_genes():
    G = {}
    with open(f'{OUT}/network_genes.tsv', encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            G[r['gene']] = r
    return G


nodes = load_nodes()
G = load_genes()
print('[载入] circRNA %d，基因 %d' % (len(nodes), len(G)))

# ============ Fig1: circRNA 火山图 ============
W, H = 900, 620
L0, T0, PW, PH = 80, 70, W - 130, H - 140
vals = []
for n in nodes:
    try:
        lfc = float(n['log2FC']); p = float(n['P'])
        vals.append((n['node_id'], n['gene'], lfc, -math.log10(max(p, 1e-300))))
    except (ValueError, KeyError):
        pass
mx = max(abs(v[2]) for v in vals) * 1.08
my = min(max(v[3] for v in vals), 8)
svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">' % (W, H, W, H),
       '<rect width="%d" height="%d" fill="#fff"/>' % (W, H),
       txt(24, 30, 'Figure 1  Differentially expressed circRNAs in HCC (GSE122482)', 17, weight='bold'),
       txt(24, 50, 'Paired tumour vs peritumour, n = 10 pairs; threshold |log2FC| >= 1 and P < 0.05 (unadjusted)', 11, fill='#666')]
PX = lambda l: L0 + (l + mx) / (2 * mx) * PW
PY = lambda q: T0 + PH - min(q, my) / my * PH
svg += ['<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#bbb"/>' % (PX(-mx), T0 + PH, PX(mx), T0 + PH),
        '<line x1="%.1f" y1="%d" x2="%.1f" y2="%.1f" stroke="#bbb"/>' % (PX(0), T0, PX(0), T0 + PH)]
# 阈值线
svg.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#e88" stroke-dasharray="4,4"/>'
           % (PX(-1), T0, PX(-1), T0 + PH))
svg.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#e88" stroke-dasharray="4,4"/>'
           % (PX(1), T0, PX(1), T0 + PH))
svg.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#e88" stroke-dasharray="4,4"/>'
           % (L0, PY(-math.log10(0.05)), L0 + PW, PY(-math.log10(0.05))))
for cid, g, lfc, q in vals:
    down = lfc < -1; up = lfc > 1
    col = '#d94f4f' if up else ('#3b7dd8' if down else '#c9ced6')
    svg.append('<circle cx="%.1f" cy="%.1f" r="5" fill="%s" fill-opacity="0.75" stroke="%s" stroke-width="1"/>'
               % (PX(lfc), PY(q), col, col))
# 标注多 circRNA 基因
multi = [g for g, r in G.items() if int(r['n_circRNA']) > 1]
lab = {}
placed_lab = []
for cid, g, lfc, q in vals:
    if g in multi and g not in lab:
        lab[g] = (PX(lfc), PY(q))
for g, (x, y) in sorted(lab.items(), key=lambda kv: -abs(kv[1][1] - T0 - PH / 2)):
    svg.append('<circle cx="%.1f" cy="%.1f" r="10" fill="none" stroke="#e08214" stroke-width="2"/>' % (x, y))
    dy = -12
    for cand in (-12, 16, -26, 30):
        if all(abs((x, y + cand) - (px, py)) > 26 for px, py in placed_lab):
            dy = cand; break
    placed_lab.append((x, y + dy))
    svg += ['<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#e08214" stroke-width="0.9"/>'
            % (x + 10, y, x + 13, y + dy + 3),
            txt(x + 15, y + dy + 3, g, 11, fill='#b35c00', weight='bold')]
svg += [txt(L0 + PW / 2, H - 28, 'log2 fold change (tumour / peritumour)', 13, anchor='middle'),
        '<text x="24" y="%.1f" font-size="13" fill="#222" font-family="Helvetica" transform="rotate(-90 24 %.1f)" text-anchor="middle">-log10(P)</text>' % (T0 + PH / 2, T0 + PH / 2),
        txt(L0 + 8, T0 + 16, 'n = %d circRNAs' % len(vals), 11, fill='#666'),
        '<circle cx="%d" cy="%d" r="5" fill="#d94f4f"/><text x="%d" y="%d" font-size="11" fill="#222" font-family="Helvetica">up-regulated</text>' % (W - 210, 36, W - 199, 40),
        '<circle cx="%d" cy="%d" r="5" fill="#3b7dd8"/><text x="%d" y="%d" font-size="11" fill="#222" font-family="Helvetica">down-regulated</text>' % (W - 210, 54, W - 199, 58),
        '<circle cx="%d" cy="%d" r="10" fill="none" stroke="#e08214" stroke-width="2"/><text x="%d" y="%d" font-size="11" fill="#222" font-family="Helvetica">host gene with multiple circRNAs</text>' % (W - 300, 74, W - 286, 78),
        '</svg>']
open(f'{FIG}/Fig1_circRNA_volcano.svg', 'w', encoding='utf-8').write('\n'.join(svg))
print('[SAVED] Fig1_circRNA_volcano.svg')


# ============ Fig2: 网络结构图 ============
W, H = 1000, 880
cx, cy, R1, R2 = W / 2, H / 2 + 10, 300, 195
svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">' % (W, H, W, H),
       '<rect width="%d" height="%d" fill="#fff"/>' % (W, H),
       txt(24, 30, 'Figure 2  circRNA-mRNA network in HCC (103 circRNAs, 91 genes)', 17, weight='bold'),
       txt(24, 50, 'Outer ring: circRNAs (red = up, blue = down)   Inner ring: host genes (colour = prognostic significance, size = number of circRNAs)', 11, fill='#666')]
# 只画有预后数据或免疫关联的基因，避免拥挤
key_genes = [g for g, r in G.items()
             if (r['P_cox_TCGA'] and r['P_cox_TCGA'] not in ('', 'nan') and float(r['P_cox_TCGA']) < 0.05)
             or float(r['immune_max_r'] or 0) > 0.2]
gnodes = [n for n in nodes if n['gene'] in key_genes]
print('  绘制 %d circRNA / %d 基因' % (len(gnodes), len(key_genes)))
import math as _m
ng = len(key_genes)
gene_ang = {}
for i, g in enumerate(sorted(key_genes)):
    a = 2 * _m.pi * i / ng - _m.pi / 2
    gene_ang[g] = a
# circRNA 角度
circ_by_gene = defaultdict(list)
for n in gnodes:
    circ_by_gene[n['gene']].append(n)
for g, lst in circ_by_gene.items():
    a0 = gene_ang[g]
    spread = 0.30
    k = len(lst)
    for j, n in enumerate(lst):
        off = (j - (k - 1) / 2) * (spread / max(k, 1))
        na = a0 + off
        nx = cx + R1 * _m.cos(na); ny = cy + R1 * _m.sin(na)
        gx = cx + R2 * _m.cos(a0); gy = cy + R2 * _m.sin(a0)
        col = '#d94f4f' if float(n['log2FC']) > 0 else '#3b7dd8'
        svg.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-opacity="0.32" stroke-width="%.1f"/>'
                   % (nx, ny, gx, gy, col, min(3.5, 0.8 + abs(float(n['log2FC'])))))
        svg.append('<circle cx="%.1f" cy="%.1f" r="4.5" fill="%s"/>' % (nx, ny, col))
# 基因节点
# 先算出所有标签位置，再对角度密集者做水平错位，避免重叠
glabels = []
for g in sorted(key_genes):
    a = gene_ang[g]
    r = G[g]
    try:
        pco = float(r['P_cox_TCGA'])
    except (ValueError, KeyError):
        pco = None
    if pco is not None and pco < 0.001:
        col = '#8b1a1a'
    elif pco is not None and pco < 0.01:
        col = '#d94f4f'
    elif pco is not None and pco < 0.05:
        col = '#e8a0a0'
    else:
        col = '#b0b6be'
    rr = 7 + min(9, int(r['n_circRNA']) * 4)
    x = cx + R2 * _m.cos(a); y = cy + R2 * _m.sin(a)
    svg.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s" stroke="#fff" stroke-width="1.5"/>' % (x, y, rr, col))
    ca, sa = _m.cos(a), _m.sin(a)
    lr = R2 + rr + 10
    glabels.append([g, cx + lr * ca, cy + lr * sa, ca, sa, col, rr])

# 顶部/底部（|cos| 小）标签做垂直错位，避免重叠
TOPBOT = [t for t in glabels if abs(t[3]) < 0.62]
TOP = sorted([t for t in TOPBOT if t[4] < 0], key=lambda t: t[1])
BOT = sorted([t for t in TOPBOT if t[4] >= 0], key=lambda t: t[1])
MIN_GAP = 15


def spread(items, bottom_up):
    """沿 y 方向拉开间距；bottom_up=True 时从下往上排，保证不出画布"""
    if not items:
        return
    if bottom_up:
        items.sort(key=lambda t: -t[2])
        for i in range(len(items) - 2, -1, -1):
            if items[i + 1][2] - items[i][2] < MIN_GAP:
                items[i][2] = items[i + 1][2] - MIN_GAP
        # 底部不得越界
        over = H - 34 - min(t[2] for t in items)
        if over < 0:
            for t in items:
                t[2] += over
    else:
        for i in range(1, len(items)):
            if items[i][2] - items[i - 1][2] < MIN_GAP:
                items[i][2] = items[i - 1][2] + MIN_GAP
        # 顶部不得越界（标题下方）
        under = min(t[2] for t in items) - 68
        if under < 0:
            for t in items:
                t[2] -= under


spread(TOP, False)
spread(BOT, True)

for g, lx, ly, ca, sa, col, rr in glabels:
    if ca > 0.35:
        anc = 'start'; tx, ty = lx, ly + 4
    elif ca < -0.35:
        anc = 'end'; tx, ty = lx, ly + 4
    else:
        anc = 'middle'; tx = lx; ty = ly + (-5 if sa < 0 else 13)
    svg.append(txt(tx, ty, g, 11, anchor=anc, weight='bold'))
    ex = cx + (R2 + rr + 2) * ca; ey = cy + (R2 + rr + 2) * sa
    dx = 4 if anc == 'start' else (-4 if anc == 'end' else 0)
    svg.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="0.8" stroke-opacity="0.55"/>'
               % (ex, ey, tx + dx, ty - 3, col))
svg += [txt(cx, cy - 6, '%d circRNAs' % len(gnodes), 15, anchor='middle', weight='bold'),
        txt(cx, cy + 16, '%d host genes' % len(key_genes), 15, anchor='middle', fill='#666'),
        txt(cx, cy + 36, 'gene node size = number of circRNAs', 11, anchor='middle', fill='#666'),
        txt(L0, H - 18, 'Gene node colour: dark red P_cox<0.001 | red P<0.01 | light red P<0.05 | grey not significant', 11, fill='#666'),
        '</svg>']
open(f'{FIG}/Fig2_network.svg', 'w', encoding='utf-8').write('\n'.join(svg))
print('[SAVED] Fig2_network.svg')


# ============ Fig3: 跨队列预后一致性 ============
pairs = []
for g, r in G.items():
    try:
        a = float(r['HR_TCGA']); b = float(r['HR_14520'])
    except (ValueError, KeyError):
        continue
    if a == a and b == b:
        pairs.append((g, a, b))
print('  可配对基因: %d' % len(pairs))
W, H = 700, 640
L0, T0, PW, PH = 90, 70, W - 150, H - 150
lo = min(min(p[1], p[2]) for p in pairs) * 0.92
hi = max(max(p[1], p[2]) for p in pairs) * 1.08
if lo > 0.3: lo = 0.3
if hi < 2.0: hi = 2.0
PX3 = lambda v: L0 + (v - lo) / (hi - lo) * PW
PY3 = lambda v: T0 + PH - (v - lo) / (hi - lo) * PH
svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">' % (W, H, W, H),
       '<rect width="%d" height="%d" fill="#fff"/>' % (W, H),
       txt(24, 30, 'Figure 3  Cross-cohort concordance of prognostic effect', 17, weight='bold'),
       txt(24, 50, 'TCGA-LIHC (n=275) vs GSE14520 (n=242), hazard ratio per +1 SD', 11, fill='#666')]
svg.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" fill-opacity="0.5"/>'
           % (PX3(1), T0, PX3(hi) - PX3(1), PY3(lo) - T0, '#e8f4ff'))
svg += ['<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#999" stroke-dasharray="5,4"/>' % (PX3(1), T0, PX3(1), T0 + PH),
        '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#999" stroke-dasharray="5,4"/>' % (L0, PY3(1), L0 + PW, PY3(1))]
svg.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#333" stroke-width="1.6"/>'
           % (PX3(lo), PY3(lo), PX3(hi), PY3(hi)))
# 一致性计数
nc = sum(1 for _, a, b in pairs if (a > 1) == (b > 1))
for g, a, b in sorted(pairs, key=lambda x: -abs(x[1] - x[2])):
    same = (a > 1) == (b > 1)
    col = '#2e8b57' if same else '#c0392b'
    svg.append('<circle cx="%.1f" cy="%.1f" r="5.5" fill="%s" fill-opacity="0.8"/>' % (PX3(a), PY3(b), col))
# 标注：优先 extremes + 已验证的 6 个基因，同角度邻近者错位
validated = {'HSD17B6', 'CD5', 'BCLAF1', 'SREBF2', 'NT5DC2', 'LRCH1'}
cands_lab = sorted(pairs, key=lambda x: x[1])[:4] + sorted(pairs, key=lambda x: -x[1])[:4] \
    + [p for p in pairs if p[0] in validated]
placed = []
for g, a, b in cands_lab:
    same = (a > 1) == (b > 1)
    col = '#1a6b45' if same else '#8b2828'
    x0, y0 = PX3(a), PY3(b)
    dx, dy = 9, -8
    # 简单避让：与已放置标签太近则换偏移
    for cand_off in [(9, -8), (9, 12), (-9, -8), (-9, 12), (9, -20), (-9, 20), (0, -24)]:
        tx, ty = x0 + cand_off[0], y0 + cand_off[1]
        if all((tx - px) ** 2 + (ty - py) ** 2 > 78 ** 2 for px, py in placed):
            dx, dy = cand_off
            break
    placed.append((x0 + dx, y0 + dy))
    anc = 'start' if dx >= 0 else 'end'
    svg.append(txt(x0 + dx, y0 + dy, g, 11, anchor=anc, fill=col, weight='bold'))
for v in [0.5, 0.7, 1.0, 1.5, 2.0]:
    if lo <= v <= hi:
        svg.append(txt(PX3(v), T0 + PH + 18, '%.1f' % v, 11, anchor='middle', fill='#666'))
        svg.append(txt(L0 - 8, PY3(v) + 4, '%.1f' % v, 11, anchor='end', fill='#666'))
svg += [txt(L0 + PW / 2, H - 24, 'Hazard ratio in TCGA-LIHC', 13, anchor='middle'),
        '<text x="26" y="%.1f" font-size="13" fill="#222" font-family="Helvetica" transform="rotate(-90 26 %.1f)" text-anchor="middle">Hazard ratio in GSE14520</text>' % (T0 + PH / 2, T0 + PH / 2),
        txt(L0 + 12, T0 + 20, 'Concordant: %d / %d (%.0f%%)' % (nc, len(pairs), nc / len(pairs) * 100), 13, weight='bold', fill='#2e8b57'),
        txt(L0 + 12, T0 + 38, 'green = same direction, red = opposite', 11, fill='#666'),
        '</svg>']
open(f'{FIG}/Fig3_concordance.svg', 'w', encoding='utf-8').write('\n'.join(svg))
print('[SAVED] Fig3_concordance.svg')
print('\n全部图表已生成于', FIG)