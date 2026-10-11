# -*- coding: utf-8 -*-
"""
HCC 全景 circRNA 网络构建 — 第 2 步
候选基因在 TCGA-LIHC 中的：表达差异 / 生存预后 / 分期关联
输出: results/tcga_candidate_gene_stats.tsv
"""
import gzip, json, math, os

RAW = r'C:\Users\hewan\Documents\kimi\tasks\2026-10-08\11-38-08-f56eebff\rawdata'
FPKM = os.path.join(RAW, 'TCGA-LIHC.star_fpkm.tsv.gz')
CLIN = os.path.join(RAW, 'TCGA-LIHC.clinical.tsv.gz')
SURV = os.path.join(RAW, 'TCGA-LIHC.survival.tsv.gz')
OUT = 'results'


# ---------- 统计工具 ----------
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


def welch(a, b):
    na, nb = len(a), len(b)
    ma, mb = sum(a) / na, sum(b) / nb
    va = sum((x - ma) ** 2 for x in a) / (na - 1)
    vb = sum((x - mb) ** 2 for x in b) / (nb - 1)
    se = (va / na + vb / nb) ** 0.5
    if se == 0: return 0.0, float('nan')
    t = (mb - ma) / se
    df = (va / na + vb / nb) ** 2 / ((va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1))
    return t, df


def km(times, events):
    pts = sorted(set(times)); surv = 1.0; risk = len(times); out = []
    for t in pts:
        d = sum(1 for ti, ei in zip(times, events) if ti == t and ei == 1)
        c = sum(1 for ti, ei in zip(times, events) if ti == t and ei == 0)
        if risk > 0 and d + c > 0:
            surv *= (1 - d / risk); out.append((t, surv))
        risk -= (d + c)
    return out


def km_median(times, events):
    for t, s in km(times, events):
        if s <= 0.5:
            return t
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


# ---------- 读候选基因 ----------
cand = {}
with open(f'{OUT}/candidate_circRNA.tsv', encoding='utf-8') as f:
    f.readline()
    for l in f:
        p = l.rstrip('\n').split('\t')
        g = p[1]
        cand.setdefault(g, []).append((p[0], float(p[4]), float(p[5])))
print(f'[1/4] 候选基因 {len(cand)} 个，关联 circRNA {sum(len(v) for v in cand.values())} 条')

# ---------- 读表达矩阵（只提候选基因）----------
print('[2/4] 扫描 TCGA-LIHC FPKM（只提取候选基因）...')
need = set(cand)
expr = {}
with gzip.open(FPKM, 'rt', encoding='utf-8', errors='ignore') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    samples = hdr[1:]
    for l in f:
        p = l.rstrip('\n').split('\t')
        gid = p[0].split('.')[0]
        try:
            expr[gid] = [float(x) for x in p[1:]]
        except ValueError:
            continue
print(f'      基因 {len(expr)} × 样本 {len(samples)}')

# Ensembl ID → symbol：用 GTF 映射（只查需要的）
import re
GTF = os.path.join(RAW, 'annotation_cache', 'gencode.v22.annotation.gtf.gz')
id2sym = {}
with gzip.open(GTF, 'rt', encoding='utf-8', errors='ignore') as f:
    for l in f:
        if l.startswith('#'): continue
        p = l.split('\t')
        if len(p) < 9 or p[2] != 'gene': continue
        mg = re.search(r'gene_id "([^"]+)"', p[8])
        mn = re.search(r'gene_name "([^"]+)"', p[8])
        if mg and mn:
            id2sym[mg.group(1).split('.')[0]] = mn.group(1)
print(f'      Ensembl→symbol 映射 {len(id2sym)} 条')

# 转换表达矩阵为 symbol
sy_expr = {}
for gid, v in expr.items():
    s = id2sym.get(gid)
    if s and s in need:
        sy_expr.setdefault(s, v)
print(f'      命中候选基因 {len(sy_expr)} / {len(need)}')

# ---------- 样本分组（barcode 末段为 01A / 11A 格式）----------
# 🔴 barcode 第三、四位编码样本类型：01 原发肿瘤、02 复发、11 正常组织。
#原判定 `seg[0] == '0'` 会把 02 复发样本一并算成肿瘤，使本步的单因素
#  Cox / log-rank 分析比校正模型多纳入 3 例复发样本，也与多因素模型口径不一致。
# 现改为精确判定 01；复发样本与多因素模型（step14_export_for_R.py）一并排除。
def stype(b):
    seg = b.split('-')[-1]
    if not seg: return 'other'
    if seg[0] == '0' and seg[1:2] == '1': return 'tumor'      # 01 = Primary Tumor
    if seg[0] == '0': return 'recurrent'                     # 02 / 03 / 04
    if seg[0] == '1' and len(seg) > 1 and seg[1].isdigit(): return 'normal'
    return 'other'

tum = [i for i, s in enumerate(samples) if stype(s) == 'tumor']
nor = [i for i, s in enumerate(samples) if stype(s) == 'normal']
rec = [i for i, s in enumerate(samples) if stype(s) == 'recurrent']
print(f'      原发肿瘤 {len(tum)}，癌旁 {len(nor)}，复发（已排除）{len(rec)}')

# ---------- 临床 ----------
C = [l.rstrip('\n').split('\t') for l in gzip.open(CLIN, 'rt', encoding='utf-8', errors='ignore').read().splitlines()]
ch = C[0]
i_vs = ch.index('vital_status.demographic')
i_dtd = ch.index('days_to_death.demographic')
i_dtf = ch.index('days_to_last_follow_up.diagnoses')
i_stage = ch.index('ajcc_pathologic_stage.diagnoses')
i_gender = ch.index('gender.demographic')
i_age = ch.index('age_at_index.demographic')
clin = {}
for r in C[1:]:
    if len(r) <= max(i_vs, i_dtd, i_dtf, i_stage):
        continue
    clin[r[0]] = {
        'vs': r[i_vs], 'dtd': r[i_dtd], 'dtf': r[i_dtf],
        'stage': r[i_stage], 'gender': r[i_gender], 'age': r[i_age],
    }

S = [l.rstrip('\n').split('\t') for l in gzip.open(SURV, 'rt', encoding='utf-8', errors='ignore').read().splitlines()]
surv = {r[0]: (float(r[1]), int(r[2])) for r in S[1:] if len(r) >= 3 and r[1] not in ('', 'NA')}

tw = {}
matched = 0
for i in tum:
    s = samples[i]
    if s in clin:
        tw[i] = (clin[s], surv.get(s))
        matched += 1

print(f'      肿瘤样本 {len(tum)}，精确匹配临床 {matched}')
evt = sum(1 for i in tw if tw[i][1] and tw[i][1][1] == 1)
print(f'      死亡事件 {evt}')

# ---------- 分析 ----------
print('\n[3/4] 计算差异表达与生存关联...')
rows = []
for g, cl in sorted(cand.items()):
    if g not in sy_expr:
        continue
    v = sy_expr[g]
    a = [v[i] for i in nor]; b = [v[i] for i in tum]
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    t, df = welch(b, a); p_diff = sf(t, df)
    lfc = mb - ma

    # 生存
    ts, es, xs = [], [], []
    for i, (c, sv) in tw.items():
        if not sv or c['vs'] == 'Not Reported': continue
        if c['vs'] == 'Dead':
            if not c['dtd'] or c['dtd'] in ('', 'NA'): continue
            ttime, e = float(c['dtd']), 1
        elif c['vs'] == 'Alive':
            if not c['dtf'] or c['dtf'] in ('', 'NA'): continue
            ttime, e = float(c['dtf']), 0
        else:
            continue
        if ttime <= 0: continue
        ts.append(ttime); es.append(e); xs.append(v[i])
    n_surv = len(ts)
    lr_p = cox_p = float('nan'); hr = float('nan')
    chi2 = float('nan'); med_hi = med_lo = None; nh = nl = 0
    if n_surv >= 80 and sum(es) > 20:
        med = sorted(xs)[len(xs) // 2]
        hi = [(tt, ee) for tt, ee, xx in zip(ts, es, xs) if xx > med]
        lo = [(tt, ee) for tt, ee, xx in zip(ts, es, xs) if xx <= med]
        nh, nl = len(hi), len(lo)
        if nh >= 5 and nl >= 5:
            chi2, lr_p = logrank([x[0] for x in hi], [x[1] for x in hi],
                                 [x[0] for x in lo], [x[1] for x in lo])
            med_hi = km_median([x[0] for x in hi], [x[1] for x in hi])
            med_lo = km_median([x[0] for x in lo], [x[1] for x in lo])
        m = sum(xs) / len(xs)
        sd_ = (sum((x - m) ** 2 for x in xs) / (len(xs) - 1)) ** 0.5
        z = [(x - m) / sd_ for x in xs]
        hr, cox_chi, cox_p = cox_score(ts, es, z)
    # 分期
    e1 = [v[i] for i, (c, sv) in tw.items() if c['stage'] == 'Stage I']
    l4 = [v[i] for i, (c, sv) in tw.items() if c['stage'] and c['stage'].startswith('Stage IV')]
    ts4 = tv4 = pn4 = float('nan')
    if len(e1) >= 5 and len(l4) >= 5:
        ts4, tv4 = sum(e1) / len(e1), sum(l4) / len(l4)
        t4, df4 = welch(l4, e1); pn4 = sf(t4, df4)
    rows.append(dict(gene=g, ncirc=len(cl), ma=ma, mb=mb, lfc=lfc,
                     p_diff=p_diff, n_surv=n_surv, n_evt=sum(es),
                     hr=hr, cox_p=cox_p, lr_chi=chi2, lr_p=lr_p,
                     nh=nh, nl=nl, med_hi=med_hi, med_lo=med_lo,
                     s1=ts4, s4=tv4, p_stage=pn4,
                     lfc_circ=[x[1] for x in cl]))
    print('  %-14s ncirc=%d log2FC=%+.3f P_diff=%.3g | HR=%.3f P_cox=%.4f P_lr=%.4f'
          % (g, len(cl), lfc, p_diff, hr, cox_p, lr_p))

with open(f'{OUT}/tcga_candidate_gene_stats.tsv', 'w', encoding='utf-8') as fo:
    fo.write('gene\tn_circRNA\tmean_normal\tmean_tumor\tlog2FC\tP_diff\t'
             'n_survival\tn_events\tHR_per_SD\tP_cox\tlogrank_chi2\tP_logrank\t'
             'n_high\tn_low\tmedianOS_high\tmedianOS_low\tstageI_mean\tstageIV_mean\tP_stage\t'
             'circRNA_log2FC\n')
    for r in rows:
        fo.write('%s\t%d\t%.4f\t%.4f\t%.4f\t%.6g\t%d\t%d\t%.4f\t%.6g\t%.4f\t%.6g\t'
                 '%d\t%d\t%s\t%s\t%.4f\t%.4f\t%.6g\t%s\n'
                 % (r['gene'], r['ncirc'], r['ma'], r['mb'], r['lfc'], r['p_diff'],
                    r['n_surv'], r['n_evt'], r['hr'], r['cox_p'], r['lr_chi'], r['lr_p'],
                    r['nh'], r['nl'],
                    r['med_hi'] if r['med_hi'] else 'NA',
                    r['med_lo'] if r['med_lo'] else 'NA',
                    r['s1'], r['s4'], r['p_stage'],
                    ','.join('%.3f' % x for x in r['lfc_circ'])))
print(f'\n[SAVED] {OUT}/tcga_candidate_gene_stats.tsv')

# ---------- 汇总排序 ----------
print('\n[4/4] 候选基因排序（按预后 P 值）===')
val = [r for r in rows if not math.isnan(r['cox_p'])]
val.sort(key=lambda r: r['cox_p'])
print('%-14s %8s %8s %9s %9s %10s' % ('gene', 'log2FC', 'HR', 'P_cox', 'P_logrank', 'n_circ'))
for r in val[:20]:
    print('%-14s %+8.3f %8.3f %9.4f %9.4f %10d'
          % (r['gene'], r['lfc'], r['hr'], r['cox_p'], r['lr_p'], r['ncirc']))
sig = [r for r in val if r['cox_p'] < 0.05]
print(f'\n预后显著（P_cox<0.05）的候选基因：{len(sig)}')
for r in sig:
    print('   %-14s HR=%.3f P_cox=%.4f P_logrank=%.4f  circRNA数=%d'
          % (r['gene'], r['hr'], r['cox_p'], r['lr_p'], r['ncirc']))
json.dump([r['gene'] for r in sig], open(f'{OUT}/prognostic_genes.json', 'w', encoding='utf-8'),
          ensure_ascii=False, indent=1)