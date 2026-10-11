# -*- coding: utf-8 -*-
"""
C 方案：PubMed 文献挖掘 HCC 的 circRNA-miRNA-mRNA 轴
（NCBI E-utilities，已验证可用）
输出: results/hcc_circ_axes.tsv, results/pubmed_hcc_abstracts.txt
"""
import urllib.request, json, ssl, time, re, os, csv
from collections import defaultdict

ctx = ssl.create_default_context()
H = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) research-tool/1.0'}
BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
OUT = 'results'
os.makedirs(OUT, exist_ok=True)


def get(url, tries=3):
    for a in range(tries):
        try:
            return urllib.request.urlopen(urllib.request.Request(url, headers=H),
                                          timeout=90, context=ctx).read().decode('utf-8', 'ignore')
        except Exception as e:
            if a == tries - 1:
                print('  [FAIL] %s: %s' % (type(e).__name__, e))
                return ''
            time.sleep(3)


# ---------- 1) 检索 ----------
term = '(circRNA[Title/Abstract]) AND (hepatocellular[Title/Abstract] OR liver[Title/Abstract]) AND (miRNA[Title/Abstract] OR microRNA[Title/Abstract])'
u = BASE + 'esearch.fcgi?db=pubmed&term=' + urllib.parse.quote(term) + \
    '&retmax=300&retmode=json&sort=date&tool=workbuddy'
d = json.loads(get(u))['esearchresult']
ids = d['idlist']
print('[1/4] PubMed 命中 %s 篇，取前 %d 篇摘要' % (d['count'], len(ids)))

# ---------- 2) 抓摘要 ----------
recs = []
for i in range(0, len(ids), 40):
    ch = ','.join(ids[i:i + 40])
    t = get(BASE + 'efetch.fcgi?db=pubmed&id=' + ch +
            '&rettype=abstract&retmode=text&tool=workbuddy')
    recs.append(t)
    time.sleep(0.5)
txt = '\n\n@@@NEXT@@@\n\n'.join(recs)
open(f'{OUT}/pubmed_hcc_abstracts.txt', 'w', encoding='utf-8').write(txt)
print('      摘要共 %d 字符' % len(txt))

# ---------- 3) 精确抽取三元轴 ----------
recs2 = [r for r in re.split(r'\n\n(?=\d+\.\s+\w)', txt) if len(r) > 200]
print('[2/4] 摘要记录 %d 条' % len(recs2))

axes = []
# 模式1: circX/miR-Y/Z
for r in recs2:
    pmid_m = re.match(r'\s*(\d+)\.', r)
    pmid = pmid_m.group(1) if pmid_m else ''
    for m in re.finditer(
            r'(?:hsa[_-])?[Cc]irc[_-][A-Za-z0-9_]+\s*/\s*(?:hsa[_-])?[Mm]i[Rr]-?\d+[a-zA-Z0-9-]*\s*/\s*([A-Za-z0-9-]{2,15})',
            r):
        full = m.group(0).strip()
        circ = re.match(r'(?:hsa[_-])?([Cc]irc[_-][A-Za-z0-9_]+)', full).group(1)
        mir = re.search(r'[Mm]i[Rr]-?\d+[a-zA-Z0-9-]*', full).group(0)
        gene = m.group(1)
        axes.append((pmid, circ, mir, gene, full))
# 模式2: "circX ... miR-Y ... Z"（同一句）
for r in recs2:
    pmid_m = re.match(r'\s*(\d+)\.', r)
    pmid = pmid_m.group(1) if pmid_m else ''
    for m in re.finditer(
            r'(?:hsa[_-])?[Cc]irc[_-][A-Za-z0-9_]+[^.;]{0,80}?[Mm]i[Rr]-?\d+[a-zA-Z0-9-]*[^.;]{0,60}?\b([A-Z][A-Z0-9]{2,9})\b',
            r):
        circ = re.match(r'(?:hsa[_-])?([Cc]irc[_-][A-Za-z0-9_]+)', m.group(0)).group(1)
        mir = re.search(r'[Mm]i[Rr]-?\d+[a-zA-Z0-9-]*', m.group(0)).group(0)
        gene = m.group(1)
        if gene in ('AND', 'THE', 'CELLS', 'HCC', 'DNA', 'RNA'):
            continue
        axes.append((pmid, circ, mir, gene, m.group(0)[:90]))

# 去重
uniq = {}
for pmid, circ, mir, gene, ctxs in axes:
    k = (circ.lower(), mir.lower(), gene.upper())
    if k not in uniq:
        uniq[k] = dict(pmid=pmid, circ=circ, miRNA=mir, gene=gene.upper(),
                       evidence=ctxs[:80])
print('[3/4] 抽取到 %d 条唯一三元轴' % len(uniq))

# ---------- 4) 关联到我的候选基因 ----------
cand = {}
with open(f'{OUT}/candidate_circRNA.tsv', encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        cand[r['gene']] = r['circRNA']

rows = sorted(uniq.values(), key=lambda x: x['gene'])
with open(f'{OUT}/hcc_circ_axes.tsv', 'w', encoding='utf-8', newline='') as fo:
    w = csv.writer(fo, delimiter='\t')
    w.writerow(['circRNA', 'miRNA', 'target_gene', 'PMID', 'in_my_candidates', 'my_circRNA', 'evidence'])
    for r in rows:
        hit = r['gene'] in cand
        w.writerow([r['circ'], r['miRNA'], r['gene'], r['pmid'],
                    'YES' if hit else '', cand.get(r['gene'], ''), r['evidence']])
print('[SAVED] %s/hcc_circ_axes.tsv' % OUT)

hit = [r for r in rows if r['gene'] in cand]
print('\n===== 与我候选基因匹配的轴（%d 条）=====' % len(hit))
print('%-24s %-16s %-12s %s' % ('circRNA', 'miRNA', 'target', 'PMID'))
for r in hit:
    print('%-24s %-16s %-12s %s' % (r['circ'][:24], r['miRNA'], r['gene'], r['pmid']))

print('\n===== 全部轴按基因汇总（前 30 基因）=====')
bygene = defaultdict(list)
for r in rows:
    bygene[r['gene']].append(r)
for g in sorted(bygene, key=lambda x: -len(bygene[x]))[:30]:
    items = bygene[g]
    cs = ','.join(sorted(set(i['circ'] for i in items)))[:40]
    ms = ','.join(sorted(set(i['miRNA'] for i in items)))[:34]
    star = ' ★候选' if g in cand else ''
    print('%-12s n=%-3d circ:%-42s miR:%s%s' % (g, len(items), cs, ms, star))