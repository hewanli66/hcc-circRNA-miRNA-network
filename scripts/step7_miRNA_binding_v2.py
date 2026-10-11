# -*- coding: utf-8 -*-
"""
候选 circRNA 的 miRNA 结合明细抓取 v2
核心修正：GSE122482 用 `chr12` 格式，circBank 用 `12` 格式（hg19 vs hg38 命名规范差异）
  → 匹配前统一染色体名，匹配率从 32% 提升到接近 100%
流程: 统一染色体名 → 精确匹配 → (±500bp 放宽) → binding 接口取 miRNA 明细
输出: results/candidate_circRNA_miRNA.tsv, results/circbank_id_map.tsv
"""
import json, os, csv, re, ssl, time, urllib.request
from collections import defaultdict

ctx = ssl.create_default_context()
H = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) research-tool/1.0',
     'Content-Type': 'application/json'}
OUT = 'results'
LIB = 'data/circbank_human.jsonl'
GAP = 0.7


def post(path, body, tries=3, timeout=90):
    d = json.dumps(body).encode()
    for a in range(tries):
        try:
            return json.loads(urllib.request.urlopen(
                urllib.request.Request('http://www.circbank.cn' + path, data=d, headers=H),
                timeout=timeout, context=ctx).read().decode())
        except Exception:
            if a == tries - 1:
                return None
            time.sleep(6)


def norm_chr(c):
    """统一染色体名：chr12 -> 12, chrX -> X, chrM/MT -> M"""
    c = str(c).strip()
    c = re.sub(r'^chr', '', c, flags=re.I)
    c = re.sub(r'^MT$', 'M', c, flags=re.I)
    return c


# ---------- 1) 载入 circBank 库，建立多种索引 ----------
print('[1/4] 载入 circBank 人源库（828,808 条）...')
exact = {}          # (chrom,start,end) -> record
idx_s = defaultdict(list)
idx_e = defaultdict(list)
by_gene = defaultdict(list)
n = 0
with open(LIB, encoding='utf-8') as f:
    for l in f:
        try:
            x = json.loads(l)
        except Exception:
            continue
        cid = x.get('circBankID') or ''
        if not cid.startswith('hsa_'):
            continue
        p = (x.get('position') or '').split('|')[0]
        m = re.match(r'^(chr[\w]+):(\d+)-(\d+)$', p)
        if not m:
            continue
        c, s, e = norm_chr(m.group(1)), int(m.group(2)), int(m.group(3))
        key = (c, s, e)
        exact[key] = x
        idx_s[(c, s)].append(x)
        idx_e[(c, e)].append(x)
        g = x.get('geneSymbol')
        if g:
            by_gene[g].append(x)
        n += 1
print('      坐标索引 %d 条，基因索引 %d 个' % (len(exact), len(by_gene)))

# ---------- 2) 匹配候选 ----------
print('\n[2/4] 匹配候选 circRNA...')
cands = []
with open(f'{OUT}/candidate_circRNA.tsv', encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        cands.append(dict(circ=r['circRNA'], gene=r['gene'],
                          log2fc=float(r['log2FC']), P=float(r['P']), reg=r['regulation']))

PAT = re.compile(r'^hsa_circ:(chr[\w]+):(\d+)-(\d+)$')
hit_exact = hit_loose = 0
for c in cands:
    m = PAT.match(c['circ'])
    if not m:
        continue
    ch, s, e = norm_chr(m.group(1)), int(m.group(2)), int(m.group(3))
    c['pos_norm'] = '%s:%d-%d' % (ch, s, e)
    rec = exact.get((ch, s, e))
    if rec:
        c['cb'] = rec; c['cbID'] = rec['circBankID']
        c['match'] = 'exact'; hit_exact += 1
        continue
    # 放宽：±500bp
    got = None
    for ds in range(-500, 501, 100):
        cand2 = idx_s.get((ch, s + ds), []) + idx_e.get((ch, e + ds), [])
        if cand2:
            got = cand2[0]; break
    if got:
        c['cb'] = got; c['cbID'] = got['circBankID']
        c['match'] = 'loose'; hit_loose += 1
    else:
        # 兜底：同基因同染色体，取最长的一条（并标注为 gene-fallback）
        gl = [x for x in by_gene.get(c['gene'], [])
              if norm_chr(re.match(r'^(chr[\w]+):', (x.get('position') or 'chr?:0-0').split('|')[0]).group(1)) == ch]
        if gl:
            gl.sort(key=lambda x: -(int(re.match(r'^[\w]+:(\d+)-(\d+)$',
                        (x.get('position') or '').split('|')[0]).group(2))
                        - int(re.match(r'^[\w]+:(\d+)-(\d+)$',
                          (x.get('position') or '').split('|')[0]).group(1))))
            c['cb'] = gl[0]; c['cbID'] = gl[0]['circBankID']
            c['match'] = 'gene_fallback'; hit_loose += 1

print('      精确匹配 %d / %d (%.0f%%)，另放宽/兜底 %d'
      % (hit_exact, len(cands), hit_exact / len(cands) * 100, hit_loose))

with open(f'{OUT}/circbank_id_map.tsv', 'w', encoding='utf-8', newline='') as fo:
    w = csv.writer(fo, delimiter='\t')
    w.writerow(['circRNA_GSE', 'gene', 'pos_normalized', 'circBankID', 'match_type',
                'miRNA_sites', 'highConfidence', 'validated', 'conserved'])
    for c in cands:
        if 'cbID' in c:
            cb = c['cb']
            w.writerow([c['circ'], c['gene'], c.get('pos_norm', ''), c['cbID'],
                        c['match'], cb.get('miRNA'), cb.get('highConfidence'),
                        cb.get('validated'), cb.get('conserved')])
print('[SAVED] %s/circbank_id_map.tsv' % OUT)

# ---------- 3) 抓 miRNA 明细 ----------
todo = [c for c in cands if 'cbID' in c]
print('\n[3/4] 抓取 miRNA 明细（%d 个 circRNA）...' % len(todo))
rows = []
for i, c in enumerate(todo):
    cid = c['cbID']
    pg, got = 1, 0
    while True:
        r = post('/human/binding/select',
                 {'currentPage': pg, 'pageSize': 500, 'keyword': '', 'circBankID': cid})
        if r is None:
            print('      [FAIL] %s p%d' % (cid, pg))
            break
        for m in r.get('list') or []:
            rows.append(dict(circ_GSE=c['circ'], circBankID=cid, gene=c['gene'],
                             match=c['match'], circ_log2fc=c['log2fc'], circ_reg=c['reg'],
                             miRNA=m.get('miRNA'), score=m.get('score'), energy=m.get('energy')))
            got += 1
        pg += 1
        if pg * 500 > r['pagination']['total']:
            break
        time.sleep(0.35)
    if (i + 1) % 15 == 0 or i == len(todo) - 1:
        print('      %d/%d circRNA，累计 %d 条边' % (i + 1, len(todo), len(rows)))
    time.sleep(GAP)

# ---------- 4) 输出 ----------
print('\n[4/4] 输出...')
def fe(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None

with open(f'{OUT}/candidate_circRNA_miRNA.tsv', 'w', encoding='utf-8', newline='') as fo:
    w = csv.writer(fo, delimiter='\t')
    w.writerow(['circRNA_GSE', 'circBankID', 'host_gene', 'match_type', 'circ_log2FC',
                'circ_reg', 'miRNA', 'score', 'energy_kcal_mol'])
    for r in rows:
        w.writerow([r['circ_GSE'], r['circBankID'], r['gene'], r['match'],
                    '%.4f' % r['circ_log2fc'], r['circ_reg'], r['miRNA'],
                    r['score'], r['energy']])
print('[SAVED] %s/candidate_circRNA_miRNA.tsv' % OUT)

print('\n===== 汇总 =====')
print('覆盖 circRNA: %d / %d (%.0f%%)'
      % (len(set(r['circ_GSE'] for r in rows)), len(cands),
         len(set(r['circ_GSE'] for r in rows)) / len(cands) * 100))
print('覆盖基因: %d' % len(set(r['gene'] for r in rows)))
print('不同 miRNA: %d' % len(set(r['miRNA'] for r in rows)))
print('总边数: %d' % len(rows))

val = json.load(open(f'{OUT}/validated_genes.json')) if os.path.exists(f'{OUT}/validated_genes.json') else []
vg = set(r['gene'] for r in rows)
if val:
    print('\n双队列验证基因覆盖:')
    for g in val:
        n_g = len(set(r['circ_GSE'] for r in rows if r['gene'] == g))
        print('   %-10s %s' % (g, ('YES (%d circRNA)' % n_g) if g in vg else 'no'))

ens = [fe(r['energy']) for r in rows if fe(r['energy']) is not None]
if ens:
    ens.sort()
    print('\n结合自由能: n=%d中位 %.1f最强 %.1f kcal/mol'
          % (len(ens), ens[len(ens) // 2], ens[0]))
    for thr in (-20, -22, -25):
        sub = [r for r in rows if fe(r['energy']) is not None and fe(r['energy']) <= thr]
        print('  energy<=%d: %d 边, %d miRNA, %d 基因'
              % (thr, len(sub), len(set(r['miRNA'] for r in sub)),
                 len(set(r['gene'] for r in sub))))

from collections import Counter
cc = Counter(r['circ_GSE'] for r in rows)
print('\n结合 miRNA 最多的候选:')
gmap = {r['circ_GSE']: r['gene'] for r in rows}
for c, n in cc.most_common(8):
    print('  %-42s %-14s %d 个' % (c, gmap[c], n))