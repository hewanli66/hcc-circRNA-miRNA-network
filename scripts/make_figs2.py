# -*- coding: utf-8 -*-
"""Fig4 免疫关联矩阵 + Fig5 核心基因双队列 KM（复用 TCGA/GSE14520 实测数据）"""
import csv, os, math, gzip, re, json
from collections import defaultdict

OUT = 'results'
FIG = 'figures'
RAW = 'C:/Users/hewan/Documents/kimi/tasks/2026-10-08/11-38-08-f56eebff/rawdata'


def esc(s):
    return str(s).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def txt(x, y, s, size=12, anchor='start', fill='#222', weight='normal', rot=None):
    tr = ' transform="rotate(%s %.1f %.1f)"' % (rot, x, y) if rot else ''
    return ('<text x="%.1f" y="%.1f" font-size="%d" text-anchor="%s" fill="%s" '
            'font-weight="%s" font-family="Helvetica,Arial,sans-serif"%s>%s</text>'
            % (x, y, size, anchor, fill, weight, tr, esc(s)))


# ==================== Fig4 免疫关联矩阵 ====================
def _rcolor(v):
    """相关系数 → 颜色。色标与单元格共用，保证视觉一致。"""
    v = max(-1.0, min(1.0, v))
    return '#d94f4f' if v > 0 else '#3b7dd8'


def _ropacity(v):
    return 0.15 + abs(max(-1.0, min(1.0, v))) * 0.85

# 数据源：immune_association_filtered.tsv
# 该表由 step12_immune_rebuild.py 生成，已套用方法声明的保留标准
# （两队列方向一致且 |r| > 0.15 同时满足），243 行原始结果收敛到 25 行
rows = []
with open(f'{OUT}/immune_association_filtered.tsv', encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        try:
            rows.append((r['cell_type'], r['gene'],
                         float(r['r_TCGA']), float(r['r_14520'])))
        except (ValueError, KeyError):
            pass
cells = sorted(set(r[1] for r in rows))
cl = sorted(set(r[0] for r in rows))
mat = defaultdict(dict)
for ct, g, a, b in rows:
    mat[g][ct] = (a, b)
print('[Fig4] 基因 %d，细胞类型 %d，有关联 %d 组' % (len(cells), len(cl), len(rows)))

# 只画两个队列都| r|>0.15 的组合，突出稳定关联
keep = {g for g in cells if any(abs(a) > 0.15 and abs(b) > 0.15 for a, b in mat[g].values())}
cells = [g for g in cells if g in keep]
cl = [c for c in cl if any(c in mat[g] and abs(mat[g][c][0]) > 0.15 and abs(mat[g][c][1]) > 0.15
                            for g in cells)]
print('  保留基因 %d，细胞 %d' % (len(cells), len(cl)))

CW, CH = 44, 24
# T0 需大于副标题底部（约 y=70），否则旋转的列标签会压住副标题
L0, T0 = 300, 150
W = L0 + len(cl) * CW + 150
H = T0 + len(cells) * CH + 120
svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">' % (W, H, W, H),
       '<rect width="%d" height="%d" fill="#fff"/>' % (W, H),
       txt(24, 30, 'Figure 4  Immune infiltration associations (cross-cohort)', 17, weight='bold'),
       txt(24, 50, 'Cell colour: red = positive, blue = negative; intensity = |r| in TCGA-LIHC', 11, fill='#666'),
       txt(24, 66, 'Only gene-cell pairs concordant in both cohorts with |r| > 0.15 are shown', 11, fill='#666')]
# 色标：连续渐变条，与单元格着色同一函数，保证视觉一致
RMAX = 0.9
_bar_w = 240
_bx, _by = L0, H - 88
_svg = []
for _i in range(60):
    _v = -RMAX + 2 * RMAX * _i / 59.0
    _c = _rcolor(_v)
    _svg.append('<rect x="%.2f" y="%d" width="%.2f" height="14" fill="%s" '
                'fill-opacity="%.2f"/>'
                % (_bx + _bar_w * _i / 60.0, _by, _bar_w / 60.0 + 0.6,
                   _c, _ropacity(_v)))
svg.extend(_svg)
svg.append('<rect x="%.1f" y="%d" width="%d" height="14" fill="none" '
           'stroke="#bbb" stroke-width="1"/>' % (_bx, _by, _bar_w))
for _v in (-0.8, -0.4, 0, 0.4, 0.8):
    svg.append(txt(_bx + (_v + RMAX) / (2 * RMAX) * _bar_w, _by - 5,
                   '%+.1f' % _v, 10, anchor='middle', fill='#666'))
svg.append(txt(_bx + _bar_w + 12, _by + 11,
               'Pearson r (TCGA-LIHC; identical in direction in GSE14520)', 10, fill='#666'))

for j, c in enumerate(cl):
    svg.append('<text x="%.1f" y="%d" font-size="11" fill="#222" font-family="Helvetica" transform="rotate(-60 %.1f %d)" text-anchor="start">%s</text>'
               % (L0 + j * CW + CW / 2, T0 - 8, L0 + j * CW + CW / 2, T0 - 8, esc(c)))
for i, g in enumerate(cells):
    svg.append(txt(L0 - 10, T0 + i * CH + CH / 2 + 4, g, 12, anchor='end', weight='bold'))
    for j, c in enumerate(cl):
        if c in mat[g]:
            a, b = mat[g][c]
            if abs(a) > 0.15 and abs(b) > 0.15:
                v = max(-1, min(1, a))
                col = _rcolor(v)
                op = _ropacity(v)
                svg.append('<rect x="%.1f" y="%.1f" width="%d" height="%d" fill="%s" fill-opacity="%.2f" stroke="#fff" stroke-width="1"/>'
                           % (L0 + j * CW, T0 + i * CH, CW - 1, CH - 1, col, op))
                svg.append(txt(L0 + j * CW + CW / 2 - 0.5, T0 + i * CH + CH / 2 + 4,
                               '%.2f' % a, 10, 'middle',
                               '#fff' if abs(a) > 0.45 else '#333'))
svg.append('</svg>')
open(f'{FIG}/Fig4_immune_matrix.svg', 'w', encoding='utf-8').write('\n'.join(svg))
print('[SAVED] Fig4_immune_matrix.svg')


# ==================== Fig5 核心基因双队列 KM ====================
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


def km(times, events):
    pts = sorted(set(times)); surv = 1.0; risk = len(times); out = []
    for t in pts:
        d = sum(1 for ti, ei in zip(times, events) if ti == t and ei == 1)
        c = sum(1 for ti, ei in zip(times, events) if ti == t and ei == 0)
        if risk > 0 and d + c > 0:
            surv *= (1 - d / risk); out.append((t, surv))
        risk -= (d + c)
    return out


def logrank(t1, e1, t2, e2):
    times = sorted(set(t1) | set(t2)); O1 = E1 = V = 0.0
    for t in times:
        a1 = sum(1 for x in t1 if x >= t); a2 = sum(1 for x in t2 if x >= t)
        n = a1 + a2
        d1 = sum(1 for x, e in zip(t1, e1) if x == t and e == 1)
        d2 = sum(1 for x, e in zip(t2, e2) if x == t and e == 1)
        dd = d1 + d2
        if n <= 1 or dd == 0: continue
        O1 += d1; E1 += dd * a1 / n
        V += dd * (a1 / n) * (1 - a1 / n)
    if V <= 0: return float('nan')
    chi = (O1 - E1) ** 2 / V
    return chi, 1 - 0.5 * (1 + math.erf(math.sqrt(chi / 2)))


# --- 读 TCGA ---
print('\n[Fig5] 读取数据...')
id2sym = {}
with gzip.open(os.path.join(RAW, 'annotation_cache', 'gencode.v22.annotation.gtf.gz'),
               'rt', encoding='utf-8', errors='ignore') as f:
    for l in f:
        if l.startswith('#'): continue
        p = l.split('\t')
        if len(p) < 9 or p[2] != 'gene': continue
        mg = re.search(r'gene_id "([^"]+)"', p[8])
        mn = re.search(r'gene_name "([^"]+)"', p[8])
        if mg and mn: id2sym[mg.group(1).split('.')[0]] = mn.group(1)

CORE = ['NT5DC2', 'HSD17B6', 'LRCH1', 'BCLAF1']
need = set(CORE)
tc = {}
samples = []
with gzip.open(os.path.join(RAW, 'TCGA-LIHC.star_fpkm.tsv.gz'), 'rt', encoding='utf-8', errors='ignore') as f:
    hdr = f.readline().rstrip('\n').split('\t'); samples = hdr[1:]
    for l in f:
        p = l.rstrip('\n').split('\t')
        s = id2sym.get(p[0].split('.')[0])
        if s and s in need:
            try: tc[s] = [float(x) for x in p[1:1 + len(samples)]]
            except ValueError: pass
C = [l.rstrip('\n').split('\t') for l in gzip.open(os.path.join(RAW, 'TCGA-LIHC.clinical.tsv.gz'),
        'rt', encoding='utf-8', errors='ignore').read().splitlines()]
ch = C[0]
ivs, idtd = ch.index('vital_status.demographic'), ch.index('days_to_death.demographic')
idt_f = ch.index('days_to_last_follow_up.diagnoses')
clin = {r[0]: (r[ivs], r[idtd], r[idt_f]) for r in C[1:] if len(r) > idt_f}
tumidx = [i for i, s in enumerate(samples)
          if s.split('-')[-1].startswith('0') and samples[i] in clin]
print('  TCGA 肿瘤 %d' % len(tumidx))

# --- 读 GSE14520 ---
p2s = {}
with open(os.path.join(RAW, 'annotation_cache', 'gpl3921_probe_symbol.csv'), encoding='utf-8', errors='ignore') as f:
    r = csv.reader(f); next(r)
    for row in r:
        if len(row) >= 2 and row[0].strip(): p2s[row[0].strip()] = row[1].strip()
s2p = {}
for p, s in p2s.items():
    if s in need: s2p.setdefault(s, []).append(p)
g_expr = {}
with gzip.open(os.path.join(RAW, 'GSE14520-GPL3921_series_matrix.txt.gz'), 'rt', encoding='utf-8', errors='ignore') as f:
    L = f.read().splitlines()
tb = next(i for i, l in enumerate(L) if l.startswith('!series_matrix_table_begin'))
te = next(i for i, l in enumerate(L) if l.startswith('!series_matrix_table_end'))
gsm = [x.strip('"') for x in L[tb + 1].split('\t')[1:]]
cnt = {}
for l in L[tb + 2:te]:
    p = l.split('\t')
    if not p: continue
    s = p2s.get(p[0].strip('"'))
    if s not in s2p: continue
    try: v = [float(x) for x in p[1:1 + len(gsm)]]
    except ValueError: continue
    if s not in g_expr:
        g_expr[s] = list(v); cnt[s] = 1
    else:
        k = cnt[s]; g_expr[s] = [(g_expr[s][i] * k + v[i]) / (k + 1) for i in range(len(gsm))]
        cnt[s] = k + 1
cl = gzip.open(os.path.join(RAW, 'GSE14520_Extra_Supplement.txt.gz'), 'rt', encoding='utf-8', errors='ignore').read().splitlines()
gh = cl[0].split('\t')
igs, itis, iss, ism = gh.index('Affy_GSM'), gh.index('Tissue Type'), gh.index('Survival status'), gh.index('Survival months')
gclin = {r.split('\t')[igs]: (r.split('\t')[itis], r.split('\t')[iss], r.split('\t')[ism])
         for r in cl[1:] if len(r.split('\t')) > ism}
gtum = [i for i, s in enumerate(gsm) if gclin.get(s, ('',))[0] == 'Tumor']
print('  GSE14520 肿瘤 %d' % len(gtum))


def build_tcga(g):
    ts, es, xs = [], [], []
    for i in tumidx:
        vs, dtd, dtf = clin[samples[i]]
        if g not in tc: continue
        if vs == 'Dead' and dtd not in ('', 'NA'):
            ts.append(float(dtd)); es.append(1); xs.append(tc[g][i])
        elif vs == 'Alive' and dtf not in ('', 'NA'):
            ts.append(float(dtf)); es.append(0); xs.append(tc[g][i])
    return ts, es, xs


def build_14520(g):
    ts, es, xs = [], [], []
    for i in gtum:
        tis, ss, sm = gclin[gsm[i]]
        if g not in g_expr or ss in ('', 'NA', '.'): continue
        try: m = float(sm)
        except ValueError: continue
        if m <= 0: continue
        ts.append(m * 30.44); es.append(int(ss)); xs.append(g_expr[g][i])   # 月→天
    return ts, es, xs


W, H = 1060, 760
L0, T0 = 100, 95
pw, ph = (W - 160) / 2 - 40, 170
svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">' % (W, H, W, H),
       '<rect width="%d" height="%d" fill="#fff"/>' % (W, H),
       txt(24, 30, 'Figure 5  Overall survival of four cross-cohort validated genes', 17, weight='bold'),
       txt(24, 50, 'Left: TCGA-LIHC (n=275)   Right: GSE14520 (n=242, months converted to days)', 11, fill='#666')]

for pi, g in enumerate(CORE):
    col, row = pi % 2, pi // 2
    sx = L0 + col * (pw + 90); sy = T0 + row * (ph + 110)
    ts, es, xs = build_tcga(g)
    tsv, esv, xsv = build_14520(g)
    for q, (t_, e_, x_, ttl, mx_) in enumerate([(ts, es, xs, 'TCGA-LIHC', max(ts) if ts else 1),
                                              (tsv, esv, xsv, 'GSE14520', max(tsv) if tsv else 1)]):
        px = sx + q * (pw / 2 + 26); py = sy
        if not t_:
            continue
        med = sorted(x_)[len(x_) // 2]
        hi = [(a, b) for a, b, x in zip(t_, e_, x_) if x > med]
        lo = [(a, b) for a, b, x in zip(t_, e_, x_) if x <= med]
        if len(hi) < 5 or len(lo) < 5:
            continue
        chi, p = logrank([a for a, _ in hi], [b for _, b in hi],
                         [a for a, _ in lo], [b for _, b in lo])
        PXf = lambda t: px + min(t, mx_) / mx_ * (pw / 2 - 20)
        PYf = lambda s: py + ph - s * ph
        svg.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="#fbfbfb" stroke="#ddd"/>'
                   % (px, py, pw / 2 - 20, ph))
        for k in range(5):
            svg += ['<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#eee"/>'
                    % (px, py + ph - k * ph / 4, px + pw / 2 - 20, py + ph - k * ph / 4),
                    txt(px - 6, py + ph - k * ph / 4 + 4, '%d%%' % (k * 25), 10, anchor='end', fill='#666')]
            svg.append(txt(px + (pw / 2 - 20) * k / 4, py + ph + 14, '%.0f' % (mx_ * k / 4 / 365), 10, anchor='middle', fill='#666'))
        for nm, cc, arr in [('High', '#d94f4f', hi), ('Low', '#2e8b57', lo)]:
            pts = km([a for a, _ in arr], [b for _, b in arr])
            d = 'M %.1f %.1f ' % (px, PYf(1)) + ' '.join('L %.1f %.1f' % (PXf(a), PYf(s)) for a, s in pts)
            svg.append('<path d="%s" fill="none" stroke="%s" stroke-width="2.2"/>' % (d, cc))
            ev = sum(1 for _, b in arr if b == 1)
            svg.append(txt(px + pw / 2 - 26, py + 16 + (0 if nm == 'High' else 15),
                           '%s n=%d, events=%d' % (nm, len(arr), ev), 10, 'end', cc, 'bold'))
        # 顶部三行标签：基因名 / P 值 / 队列
        svg.append(txt(px + (pw / 2 - 20) / 2, py - 30, g, 13, 'middle', '#111', 'bold'))
        svg.append(txt(px + (pw / 2 - 20) / 2, py - 14,
                       'log-rank P = %.3f' % p, 11, 'middle',
                       '#c33' if p < 0.05 else '#666', 'bold'))
        svg.append(txt(px + (pw / 2 - 20) / 2, py - 2, '%s (n=%d)' % (ttl, len(t_)), 10, 'middle', '#666'))
svg += [txt(L0, H - 40, 'Green = low expression, red = high expression. Cut-off = median expression.', 11, fill='#666'),
        txt(L0, H - 20, 'X-axis: years', 12, fill='#222'),
        '</svg>']
open(f'{FIG}/Fig5_KM_core.svg', 'w', encoding='utf-8').write('\n'.join(svg))
print('[SAVED] Fig5_KM_core.svg')
print('\n全部完成')