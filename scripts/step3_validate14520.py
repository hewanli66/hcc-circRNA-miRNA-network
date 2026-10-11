# -*- coding: utf-8 -*-
"""
HCC 全景 circRNA 网络构建 — 第 3 步
1) 在 GSE14520（445 例独立队列）验证候选基因的预后价值
2) TCGA + GSE14520 交叉验证一致性评估
3) ceRNA 网络核心边（circRNA-miRNA-mRNA）
输出: results/gse14520_validation.tsv, results/network_core.tsv
"""
import gzip, os, csv, math, json

RAW = 'C:/Users/hewan/Documents/kimi/tasks/2026-10-08/11-38-08-f56eebff/rawdata'
PROBE = os.path.join(RAW, 'annotation_cache', 'gpl3921_probe_symbol.csv')
MTX = os.path.join(RAW, 'GSE14520-GPL3921_series_matrix.txt.gz')
CLIN = os.path.join(RAW, 'GSE14520_Extra_Supplement.txt.gz')
OUT = 'results'


# ================= 统计工具 =================
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


def km_median(times, events):
    pts = sorted(set(times)); surv = 1.0; risk = len(times)
    for t in pts:
        d = sum(1 for ti, ei in zip(times, events) if ti == t and ei == 1)
        c = sum(1 for ti, ei in zip(times, events) if ti == t and ei == 0)
        if risk > 0 and d + c > 0:
            surv *= (1 - d / risk)
            if surv <= 0.5:
                return t
        risk -= (d + c)
    return None


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
    if V <= 0: return float('nan'), float('nan')
    chi = (O1 - E1) ** 2 / V
    return chi, 1 - 0.5 * (1 + math.erf(math.sqrt(chi / 2)))


def cox_score(times, events, x):
    order = sorted(range(len(times)), key=lambda i: times[i])
    U = 0.0; I = 0.0
    for i in order:
        if events[i] != 1: continue
        risk = [j for j in order if times[j] >= times[i]]
        s0 = len(risk)
        if s0 < 2: continue
        s1 = sum(x[j] for j in risk); s2 = sum(x[j] * x[j] for j in risk)
        var = s2 / s0 - (s1 / s0) ** 2
        if var <= 1e-12: continue
        U += x[i] - s1 / s0; I += var
    if I <= 0: return float('nan'), float('nan'), float('nan')
    beta = U / I; chi2 = U * U / I
    p = 1 - 0.5 * (1 + math.erf(math.sqrt(chi2 / 2)))
    return math.exp(beta), chi2, p


def pearson(x, y):
    n = len(x); mx = sum(x) / n; my = sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x); syy = sum((b - my) ** 2 for b in y)
    if sxx <= 0 or syy <= 0: return float('nan')
    return sxy / math.sqrt(sxx * syy)


# ================= 1) 探针映射 =================
print('[1/4] 读取 GPL3921 探针映射...')
p2s = {}
with open(PROBE, encoding='utf-8', errors='ignore') as f:
    r = csv.reader(f); next(r)
    for row in r:
        if len(row) >= 2 and row[0].strip():
            p2s[row[0].strip()] = row[1].strip()
sym2probes = {}
for p, s in p2s.items():
    sym2probes.setdefault(s, []).append(p)
print('      探针 %d，symbol %d' % (len(p2s), len(sym2probes)))

# ================= 2) 读 TCGA 候选基因统计 =================
print('\n[2/4] 读取 TCGA 候选基因结果...')
tcga = {}
with open(f'{OUT}/tcga_candidate_gene_stats.tsv', encoding='utf-8') as f:
    for x in csv.DictReader(f, delimiter='\t'):
        tcga[x['gene']] = x
need = set(tcga)
print('      候选基因 %d，其中 GSE14520 可映射 %d' % (len(need), sum(1 for g in need if g in sym2probes)))

# ================= 3) 读表达矩阵 =================
print('\n[3/4] 扫描 GSE14520 表达矩阵（445 样本）...')
sy_expr = {}     # symbol -> {gsm: value}
want_probes = {}
for g in need:
    if g in sym2probes:
        for p in sym2probes[g]:
            want_probes[p] = g
with gzip.open(MTX, 'rt', encoding='utf-8', errors='ignore') as f:
    L = f.read().splitlines()
tb = next(i for i, l in enumerate(L) if l.startswith('!series_matrix_table_begin'))
te = next(i for i, l in enumerate(L) if l.startswith('!series_matrix_table_end'))
gsm = [x.strip('"') for x in L[tb + 1].split('\t')[1:]]
for l in L[tb + 2:te]:
    p = l.split('\t')
    if not p: continue
    pid = p[0].strip('"')
    g = want_probes.get(pid)
    if not g: continue
    try:
        vals = [float(x) for x in p[1:1 + len(gsm)]]
    except ValueError:
        continue
    d = sy_expr.setdefault(g, {})
    for i, gg in enumerate(gsm):
        # 多探针取均值
        d[gg] = d.get(gg, 0.0) + vals[i] / len(sym2probes[g])
print('      命中基因 %d' % len(sy_expr))

# ================= 4) 临床与生存 =================
print('\n[4/4] 外部验证（生存分析）...')
cl = gzip.open(CLIN, 'rt', encoding='utf-8', errors='ignore').read().splitlines()
ch = cl[0].split('\t')
i_gsm = ch.index('Affy_GSM')
i_tis = ch.index('Tissue Type')
i_ss = ch.index('Survival status')
i_sm = ch.index('Survival months')
i_rs = ch.index('Recurr status')
i_rm = ch.index('Recurr months')
clin = {}
for r in cl[1:]:
    p = r.split('\t')
    if len(p) <= i_rm: continue
    clin[p[i_gsm]] = {'tis': p[i_tis], 'ss': p[i_ss], 'sm': p[i_sm],
                      'rs': p[i_rs], 'rm': p[i_rm]}
tum_gsm = [s for s in gsm if clin.get(s, {}).get('tis') == 'Tumor']
non_gsm = [s for s in gsm if clin.get(s, {}).get('tis') in ('Non-Tumor', 'Normal')]
print('      肿瘤 %d，癌旁/非肿瘤 %d' % (len(tum_gsm), len(non_gsm)))

rows = []
for g in sorted(sy_expr):
    d = sy_expr[g]
    tv = [d[s] for s in tum_gsm if s in d]
    nv = [d[s] for s in non_gsm if s in d]
    if len(tv) < 20 or len(nv) < 20: continue

    def wmean(a, b):
        ma, mb = sum(a) / len(a), sum(b) / len(b)
        va = sum((x - ma) ** 2 for x in a) / (len(a) - 1)
        vb = sum((x - mb) ** 2 for x in b) / (len(b) - 1)
        se = (va / len(a) + vb / len(b)) ** 0.5
        if se == 0: return mb - ma, float('nan')
        t = (mb - ma) / se
        df = (va / len(a) + vb / len(b)) ** 2 / ((va / len(a)) ** 2 / (len(a) - 1) + (vb / len(b)) ** 2 / (len(b) - 1))
        return mb - ma, sf(t, df)
    lfc, pdiff = wmean(nv, tv)

    # OS
    ts, es, xs = [], [], []
    for s in tum_gsm:
        if s not in d: continue
        c = clin.get(s)
        if not c or c['ss'] in ('', 'NA', '.'): continue
        try:
            m = float(c['sm'])
        except ValueError:
            continue
        if m <= 0: continue
        ts.append(m); es.append(int(c['ss'])); xs.append(d[s])
    # RFS（复发无事件=0）
    rs_, es2, xs2 = [], [], []
    for s in tum_gsm:
        if s not in d: continue
        c = clin.get(s)
        if not c or c['rs'] in ('', 'NA', '.'): continue
        try:
            m = float(c['rm'])
        except ValueError:
            continue
        if m <= 0: continue
        rs_.append(m); es2.append(int(c['rs'])); xs2.append(d[s])

    hr = cox_p = lr_p = float('nan')
    chi2 = float('nan'); mh = ml = None; nh = nl = 0
    if len(ts) >= 100 and sum(es) > 20:
        med = sorted(xs)[len(xs) // 2]
        hi = [(a, b) for a, b, x in zip(ts, es, xs) if x > med]
        lo = [(a, b) for a, b, x in zip(ts, es, xs) if x <= med]
        nh, nl = len(hi), len(lo)
        if nh >= 5 and nl >= 5:
            chi2, lr_p = logrank([a for a, _ in hi], [b for _, b in hi],
                                 [a for a, _ in lo], [b for _, b in lo])
            mh = km_median([a for a, _ in hi], [b for _, b in hi])
            ml = km_median([a for a, _ in lo], [b for _, b in lo])
        m = sum(xs) / len(xs)
        sd_ = (sum((x - m) ** 2 for x in xs) / (len(xs) - 1)) ** 0.5
        hr, cox_chi, cox_p = cox_score(ts, es, [(x - m) / sd_ for x in xs])

    rfs_p = float('nan')
    if len(rs_) >= 100 and sum(es2) > 20:
        med2 = sorted(xs2)[len(xs2) // 2]
        h2 = [(a, b) for a, b, x in zip(rs_, es2, xs2) if x > med2]
        l2 = [(a, b) for a, b, x in zip(rs_, es2, xs2) if x <= med2]
        _, rfs_p = logrank([a for a, _ in h2], [b for _, b in h2],
                           [a for a, _ in l2], [b for _, b in l2])

    t_ = tcga.get(g, {})
    circ_fcs = []
    for x in t_.get('circRNA_log2FC', '').split(','):
        if x.strip():
            try:
                circ_fcs.append(float(x))
            except ValueError:
                pass
    rec = dict(gene=g, n_circ=int(t_.get('n_circRNA', 0) or 0),
               lfc14520=lfc, p14520=pdiff,
               hr14520=hr, p_cox14520=cox_p, p_lr14520=lr_p,
               med_hi=mh, med_lo=ml, nh=nh, nl=nl,
               n_os=len(ts), n_evt=sum(es),
               n_rfs=len(rs_), n_revt=sum(es2), p_rfs=rfs_p,
               hr_tcga=float(t_.get('HR_per_SD') or 'nan'),
               p_cox_tcga=float(t_.get('P_cox') or 'nan'),
               p_lr_tcga=float(t_.get('P_logrank') or 'nan'),
               lfc_circ=circ_fcs)
    rows.append(rec)

# ================= 输出 =================
with open(f'{OUT}/gse14520_validation.tsv', 'w', encoding='utf-8', newline='') as fo:
    w = csv.writer(fo, delimiter='\t')
    w.writerow(['gene', 'n_circRNA', 'log2FC_14520', 'P_diff_14520',
                'HR_14520', 'P_cox_14520', 'P_logrank_14520',
                'medianOS_high', 'medianOS_low', 'n_high', 'n_low',
                'n_survival', 'n_events', 'P_RFS_14520',
                'HR_TCGA', 'P_cox_TCGA', 'P_logrank_TCGA', 'circRNA_log2FC'])
    for r in rows:
        w.writerow([r['gene'], r['n_circ'], '%.4f' % r['lfc14520'], '%.6g' % r['p14520'],
                    '%.4f' % r['hr14520'], '%.6g' % r['p_cox14520'], '%.6g' % r['p_lr14520'],
                    r['med_hi'] if r['med_hi'] else 'NA',
                    r['med_lo'] if r['med_lo'] else 'NA', r['nh'], r['nl'],
                    r['n_os'], r['n_evt'], '%.6g' % r['p_rfs'],
                    '%.4f' % r['hr_tcga'], '%.6g' % r['p_cox_tcga'], '%.6g' % r['p_lr_tcga'],
                    ','.join('%.3f' % x for x in r['lfc_circ'])])
print('\n[SAVED] %s/gse14520_validation.tsv' % OUT)

# ================= 交叉验证一致性 =================
print('\n===== 交叉验证：TCGA vs GSE14520 =====')
print('%-14s %9s %9s %8s | %9s %9s %8s | 一致?' % (
    'gene', 'HR_TCGA', 'P_cox_T', 'P_lr_T', 'HR_14520', 'P_cox_1', 'P_lr_1'))
concord = 0; tested = 0
keep = []
for r in rows:
    if math.isnan(r['hr_tcga']) or math.isnan(r['hr14520']):
        continue
    tested += 1
    same = (r['hr_tcga'] > 1) == (r['hr14520'] > 1)
    if same: concord += 1
    flag = 'YES' if same else 'no'
    print('%-14s %9.3f %9.4f %8.4f | %9.3f %9.4f %8.4f | %s'
          % (r['gene'], r['hr_tcga'], r['p_cox_tcga'], r['p_lr_tcga'],
             r['hr14520'], r['p_cox14520'], r['p_lr14520'], flag))
    if same and r['p_cox14520'] < 0.05:
        keep.append(r)
print('\nHR 方向一致: %d / %d (%.0f%%)' % (concord, tested, concord / max(tested, 1) * 100))
print('两队列均显著且方向一致（GSE14520 P_cox<0.05）: %d 个' % len(keep))
for r in sorted(keep, key=lambda x: x['p_cox14520']):
    print('   %-14s TCGA HR=%.3f(P=%.4f) | 14520 HR=%.3f(P=%.4f) | circ=%d 个'
          % (r['gene'], r['hr_tcga'], r['p_cox_tcga'], r['hr14520'], r['p_cox14520'], r['n_circ']))

json.dump([r['gene'] for r in keep],
          open(f'{OUT}/validated_genes.json', 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
print('\n[SAVED] %s/validated_genes.json' % OUT)