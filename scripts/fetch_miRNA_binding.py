# -*- coding: utf-8 -*-
"""
按需抓取候选 circRNA 的 miRNA 结合明细（circBank /human/binding/select）
输入: results/candidate_circRNA.tsv（103 条候选，含基因与坐标）
流程: 坐标 → circBankID → miRNA 明细（含结合自由能）
输出: results/candidate_circRNA_miRNA.tsv, results/circbank_id_map.tsv
"""
import csv, json, os, re, ssl, time, urllib.request

ctx = ssl.create_default_context()
H = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) research-tool/1.0',
     'Content-Type': 'application/json'}
BASE = 'http://www.circbank.cn'
OUT = 'results'


def post(path, body, tries=4):
    d = json.dumps(body).encode()
    req = urllib.request.Request(BASE + path, data=d, headers=H)
    for a in range(tries):
        try:
            return json.loads(urllib.request.urlopen(req, timeout=120, context=ctx).read().decode())
        except Exception as e:
            if a == tries - 1:
                return {'code': -1, 'msg': '%s: %s' % (type(e).__name__, e),
                        'list': [], 'pagination': {'total': 0}}
            time.sleep(3)


# ---------- 1) 读候选 ----------
cands = []
with open(f'{OUT}/candidate_circRNA.tsv', encoding='utf-8') as f:
    f.readline()
    for l in f:
        p = l.rstrip('\n').split('\t')
        cands.append(dict(circ=p[0], gene=p[1], normal=float(p[2]),
                          tumor=float(p[3]), log2fc=float(p[4]),
                          P=float(p[5]), padj=float(p[6]), reg=p[7]))
print('[1/4] 候选 circRNA: %d（涉及基因 %d）'
      % (len(cands), len(set(c['gene'] for c in cands))))

# ---------- 2) 匹配 circBankID ----------
print('\n[2/4] 匹配 circBankID ...')
pos2id = {}
f_l = 'data/circbank_human.jsonl'
if os.path.exists(f_l):
    n = 0
    with open(f_l, encoding='utf-8') as f:
        for l in f:
            try:
                x = json.loads(l)
            except Exception:
                continue
            cid = x.get('circBankID') or ''
            if not cid.startswith('hsa_'):
                continue
            p = x.get('position') or ''
            # position 形如 chr14:100792125-100792565|+，前半段即坐标
            core = p.split('|')[0]
            pos2id[core] = cid
            n += 1
    print('      本地 circBank 库已载入 %d 条人源记录' % n)
else:
    print('      [警告] 未找到本地 circBank 库')

PAT = re.compile(r'^hsa_circ:(chr[\w]+):(\d+)-(\d+)$')
hit = 0
for c in cands:
    m = PAT.match(c['circ'])
    if not m:
        continue
    key = '%s:%s-%s' % (m.group(1), m.group(2), m.group(3))
    c['cb_pos'] = key
    if key in pos2id:
        c['circBankID'] = pos2id[key]
        hit += 1
print('      坐标精确匹配: %d / %d' % (hit, len(cands)))

with open(f'{OUT}/circbank_id_map.tsv', 'w', encoding='utf-8', newline='') as fo:
    w = csv.writer(fo, delimiter='\t')
    w.writerow(['circRNA', 'position', 'circBankID', 'gene'])
    for c in cands:
        w.writerow([c['circ'], c.get('cb_pos', ''), c.get('circBankID', ''), c['gene']])
print('[SAVED] %s/circbank_id_map.tsv' % OUT)

# ---------- 3) 按基因批量取 miRNA 明细 ----------
# 策略：同一基因的多个候选 circRNA 只查一次（ID 形如 hsa_<GENE>_xxxxxxx）
gene2cids = {}
for c in cands:
    cid = c.get('circBankID', '')
    if cid:
        gene2cids.setdefault(c['gene'], []).append(cid)
print('\n[3/4] 待查基因 %d 个，circBankID %d 个' % (len(gene2cids), sum(len(v) for v in gene2cids.values())))

# 先看单个基因的返回规模
test_cid = next(iter(next(iter(gene2cids.values()))), None)
if test_cid:
    r = post('/human/binding/select',
             {'currentPage': 1, 'pageSize': 5, 'keyword': '', 'circBankID': test_cid})
    print('      试查 %s → total=%d msg=%s' % (test_cid, r['pagination']['total'], r.get('msg')))
    if r['list']:
        print('      字段: %s' % list(r['list'][0].keys()))