# -*- coding: utf-8 -*-
"""
补抓 5 个漏抓 circRNA 的 miRNA binding 明细
背景：step7_miRNA_binding_v2.py 的抓取在 1925/2448 页（79%）中断，
     95 个匹配中 21 个未进网络。其中 14 个 circBank 本身无预测（正常），
     5 个属抓取遗漏（本脚本处理对象）。

漏抓清单（circBank 声明 miRNA 数）：
  hsa_HSD17B6_0000300     297  → 实际接口返回 337
  hsa_ADAMTS13_0000800    128
  hsa_ADAMTS13_0001400    178
  hsa_MOGAT2_0000300       92
  hsa_SLC25A47_0000300    127

输出：
  results/backfill_miRNA.tsv          补抓到的边（与原抓取同字段，追加到主文件用）
  data/backfill_log.jsonl             抓取日志（每条 circRNA 一行）
行为：只追加，不覆盖 results/candidate_circRNA_miRNA.tsv
"""
import json, os, csv, ssl, time, urllib.request

ctx = ssl.create_default_context()
H = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) research-tool/1.0',
     'Content-Type': 'application/json'}
OUT = 'results'
LOG = 'data/backfill_log.jsonl'

# 漏抓的 5 个：circBankID -> (GSE坐标, 宿主基因名, match_type, log2FC, regulation)
TARGETS = [
    ('hsa_HSD17B6_0000300',  'hsa_circ:chr12:56754955-56755968',    'HSD17B6',  'exact',          None, None),
    ('hsa_ADAMTS13_0000800', 'hsa_circ:chr9:136302869-136303486',   'MGC50722', 'exact',          None, None),
    ('hsa_ADAMTS13_0001400', 'hsa_circ:chr9:136304487-136310940',   'MGC50722', 'exact',          None, None),
    ('hsa_MOGAT2_0000300',   'hsa_circ:chr11:75438480-75440034',   'MOGAT2',   'exact',          None, None),
    ('hsa_SLC25A47_0000300', 'hsa_circ:chr14:100792125-100792565', 'SLC25A47', 'exact',          None, None),
]


def post(path, body, tries=4, timeout=90):
    d = json.dumps(body).encode()
    last = None
    for a in range(tries):
        try:
            return json.loads(urllib.request.urlopen(
                urllib.request.Request('http://www.circbank.cn' + path, data=d, headers=H),
                timeout=timeout, context=ctx).read().decode())
        except Exception as ex:
            last = str(ex)[:120]
            time.sleep(5 * (a + 1))
    print('      [FAIL] %s' % last)
    return None


# 补 log2FC / regulation（从候选表取）
cand = {}
with open(os.path.join(OUT, 'candidate_circRNA.tsv'), encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        cand[r['circRNA']] = (r['log2FC'], r['regulation'])

rows = []
logs = []
print('=== 补抓 %d 个 circRNA ===' % len(TARGETS))
for cid, coord, gene, mt, _a, _b in TARGETS:
    lfc, reg = cand.get(coord, ('', ''))
    print('\n[%s] %s  (%s)' % (cid, coord, gene))
    pg, got, declared = 1, 0, None
    while True:
        r = post('/human/binding/select',
                 {'currentPage': pg, 'pageSize': 500, 'keyword': '', 'circBankID': cid})
        if r is None:
            print('      抓取中断于第 %d 页，已取 %d 条' % (pg, got))
            break
        if declared is None:
            declared = r.get('pagination', {}).get('total')
        lst = r.get('list') or []
        for m in lst:
            rows.append(dict(circRNA_GSE=coord, circBankID=cid, host_gene=gene,
                             match_type=mt, circ_log2FC=lfc, circ_reg=reg,
                             miRNA=m.get('miRNA'), score=m.get('score'),
                             energy=m.get('energy')))
            got += 1
        pg += 1
        if not lst:
            break
        if r.get('pagination', {}).get('total') and pg * 500 > r['pagination']['total']:
            break
        if pg > 40:          # 安全上限
            print('      达到页数上限')
            break
        time.sleep(0.4)
    print('      声明 %s 条，实取 %d 条' % (declared, got))
    logs.append(dict(circBankID=cid, coord=coord, gene=gene, declared=declared,
                     fetched=got, pages=pg - 1))
    time.sleep(0.8)

# ---------- 保存 ----------
bf = os.path.join(OUT, 'backfill_miRNA.tsv')
with open(bf, 'w', encoding='utf-8', newline='') as fo:
    w = csv.writer(fo, delimiter='\t')
    w.writerow(['circRNA_GSE', 'circBankID', 'host_gene', 'match_type',
                'circ_log2FC', 'circ_reg', 'miRNA', 'score', 'energy'])
    for r in rows:
        w.writerow([r['circRNA_GSE'], r['circBankID'], r['host_gene'], r['match_type'],
                    r['circ_log2FC'], r['circ_reg'], r['miRNA'], r['score'], r['energy']])
with open(LOG, 'w', encoding='utf-8') as fo:
    for l in logs:
        fo.write(json.dumps(l, ensure_ascii=False) + '\n')

print('\n=== 汇总 ===')
print('补抓边数: %d' % len(rows))
print('涉及 circRNA: %d' % len(set(r['circRNA_GSE'] for r in rows)))
print('涉及 miRNA: %d' % len(set(r['miRNA'] for r in rows)))
print('[SAVED] %s' % bf)
print('[SAVED] %s' % LOG)
for l in logs:
    flag='' if l['declared']==l['fetched'] else'  ← 声明与实取不符'
    print('  %-24s 声明 %-5s 实取 %-5d%s' % (l['circBankID'], l['declared'], l['fetched'], flag))