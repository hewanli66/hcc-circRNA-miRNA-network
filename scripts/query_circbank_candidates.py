# -*- coding: utf-8 -*-
"""
GSE122482 候选 circRNA × circBank 定向匹配（按 position 前缀，单线程限速）

背景：
  全库翻页抓取（4068 页 × 平均 25s ≈ 28 小时）代价过高，且深翻页越来越慢。
  实测 /search 的 keyword 字段对 position 做前缀匹配：
    kw="chr19:58858900"    -> total=1  ['chr19:58858900-58861793|+']
    kw="chr19:5885890"     -> total=2  ['chr19:58858900-...', 'chr19:58858907-...']
  ⟹ 对每条候选直接用 "chrN:start" 做前缀查询即可拿到该起点下的全部记录，
     再在本地按"完整坐标 + 链方向"做严格判定。这与全库扫描等价，但只需 ~103 次请求。

安全阀：
  · 单线程，请求间隔 GAP 秒
  · 断点续跑（已查过的候选记入 data/cb_query_done.json）
  · 只查 103 条候选，不做全库遍历，服务器压力远低于翻页抓取

输入: results/candidate_circRNA.tsv
输出: results/circbank_targeted_match.tsv   逐候选匹配明细
      results/circbank_targeted_summary.json 汇总
"""
import urllib.request, json, ssl, time, os, re, sys

ctx = ssl.create_default_context()
H = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) research-tool/1.0',
     'Content-Type': 'application/json'}
URL = 'http://www.circbank.cn/search'
CAND = 'results/candidate_circRNA.tsv'
DONE = 'data/cb_query_done.json'
OUT = 'results/circbank_targeted_match.tsv'
SUM = 'results/circbank_targeted_summary.json'
GAP = 1.5          # 请求间隔（秒），保持礼貌
TIMEOUT = 60


def post(body, tries=3):
    d = json.dumps(body).encode()
    for a in range(tries):
        try:
            return json.loads(urllib.request.urlopen(
                urllib.request.Request(URL, data=d, headers=H),
                timeout=TIMEOUT, context=ctx).read().decode())
        except Exception as e:
            if a == tries - 1:
                print('    [FAIL] %s' % e, flush=True)
                return None
            time.sleep(4 * (a + 1))


def bare(cid):
    return cid.split(':', 1)[1] if cid.startswith('hsa_circ:') else cid


def norm(p):
    m = re.match(r'^(chr[0-9XYM]+):(\d+)-(\d+)\|([+-])$', (p or '').strip())
    return (m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)) if m else None


# ---------- 读候选 ----------
cands = []
with open(CAND, encoding='utf-8') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    for line in f:
        line = line.rstrip('\n')
        if line.strip():
            cands.append(dict(zip(hdr, line.split('\t'))))
print('[CAND] %d 条候选' % len(cands), flush=True)

# ---------- 断点 ----------
done = {}
if os.path.exists(DONE):
    try:
        done = json.load(open(DONE, encoding='utf-8'))
    except Exception:
        done = {}
print('[RESUME] 已完成 %d 条' % len(done), flush=True)

rows = []
t0 = time.time()
nq = 0

for i, cd in enumerate(cands, 1):
    cid = cd['circRNA']
    coord = bare(cid)
    chrom, st, en = coord.split(':')[0], int(coord.split(':')[1].split('-')[0]), int(coord.split('-')[1])

    if cid in done:
        hits = done[cid]
    else:
        # 前缀查询：该染色体+起点的全部记录
        r = post({'currentPage': 1, 'pageSize': 50,
                  'keyword': '%s:%d' % (chrom, st), 'circBankID': ''})
        nq += 1
        if r is None:
            print('  [%d/%d] %s 查询失败，跳过' % (i, len(cands), cid), flush=True)
            time.sleep(GAP)
            continue
        hits = [x for x in (r.get('list') or [])
                if (norm(x.get('position')) or ('', -1, -1, ''))[0] == chrom
                and (norm(x.get('position')) or ('', -1, -1, ''))[1] == st]
        hits.sort(key=lambda x: int((norm(x.get('position')) or ('', 0, 0, ''))[2]))
        done[cid] = hits
        json.dump(done, open(DONE, 'w', encoding='utf-8'), ensure_ascii=False)
        if i % 20 == 0:
            el = time.time() - t0
            print('  [%d/%d] 新查询 %d 次 | %.0f 秒' % (i, len(cands), nq, el), flush=True)
        time.sleep(GAP)

    # 本地严格判定：完整坐标 + 链方向必须完全一致
    exact, variant = [], []
    for x in hits:
        pp = norm(x.get('position'))
        if not pp:
            continue
        c, s, e, st_ = pp
        if s == st and e == en:
            exact.append(x)
        elif s == st:
            variant.append(x)

    for x in exact + variant:
        c, s, e, st_ = norm(x['position'])
        rows.append(dict(
            circRNA=cid, gene=cd.get('gene', ''),
            log2FC=cd.get('log2FC', ''), P_adj=cd.get('P_adj', ''),
            regulation=cd.get('regulation', ''),
            match='exact' if (s, e) == (st, en) else 'same_start_diff_end',
            circBankID=x.get('circBankID', ''), cb_position=x.get('position', ''),
            cb_gene=x.get('geneSymbol', ''), cb_geneID=x.get('geneID', ''),
            circType=x.get('circType', ''), splicedLength=x.get('splicedLength', ''),
            samples=x.get('samples', ''), miRNA=x.get('miRNA', ''),
            m6A=x.get('m6A', ''), ires=x.get('ires', ''), orf=x.get('orf', ''),
            conserved=x.get('conserved', ''), validated=x.get('validated', ''),
            source=x.get('source', ''),
        ))

# ---------- 输出 ----------
cols = ['circRNA', 'gene', 'log2FC', 'P_adj', 'regulation', 'match',
        'circBankID', 'cb_position', 'cb_gene', 'cb_geneID', 'circType',
        'splicedLength', 'samples', 'miRNA', 'm6A', 'ires', 'orf',
        'conserved', 'validated', 'source']
os.makedirs('results', exist_ok=True)
# 注意：同一候选可能命中多条 circBank 记录，必须全部落盘，
# 不可用 {circRNA: row} 折叠（那会丢掉同起点异终点的多条记录）
byc = {}
for r in rows:
    byc.setdefault(r['circRNA'], []).append(r)
with open(OUT, 'w', encoding='utf-8', newline='') as f:
    f.write('\t'.join(cols) + '\n')
    for cd in cands:                       # 按候选原顺序输出，命中与未命中都留痕
        cid = cd['circRNA']
        hits = byc.get(cid)
        if hits:
            for r in hits:
                f.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')
        else:
            f.write('\t'.join([cid, cd.get('gene', ''), cd.get('log2FC', ''),
                               cd.get('P_adj', ''), cd.get('regulation', ''),
                               'no_hit'] + [''] * 14) + '\n')
n_rows = sum(len(v) for v in byc.values()) + len(cands) - len(byc)

n_exact = len({r['circRNA'] for r in rows if r['match'] == 'exact'})
n_var = len({r['circRNA'] for r in rows if r['match'] == 'same_start_diff_end'})
n_any = len({r['circRNA'] for r in rows})
gene_ok = len({r['circRNA'] for r in rows
               if r['match'] == 'exact' and r['cb_gene'] and r['gene'] and
               r['cb_gene'].upper() == r['gene'].upper()})

summary = dict(
    candidates=len(cands),
    queries_issued=nq,
    matched_any=n_any,
    matched_exact=n_exact,
    matched_same_start_diff_end=n_var,
    no_hit=len(cands) - n_any,
    match_rate_any=round(n_any / len(cands) * 100, 2),
    match_rate_exact=round(n_exact / len(cands) * 100, 2),
    exact_and_gene_symbol_same=gene_ok,
    elapsed_sec=round(time.time() - t0, 1),
    method='circBank /search keyword=chrN:start (position prefix), local strict verify',
)
json.dump(summary, open(SUM, 'w', encoding='utf-8'), ensure_ascii=False, indent=2)

print('\n=====匹配结果=====')
for k, v in summary.items():
    print('  %-32s %s' % (k, v))
print('\n[SAVED] %s (%d 行)' % (OUT, n_rows + 1))
print('[SAVED] %s' % SUM)