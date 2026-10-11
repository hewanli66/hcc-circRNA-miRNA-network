# -*- coding: utf-8 -*-
"""
circBank 人源 circRNA 温和抓取（单线程 + 限速 + 长退避）
针对此前并发压垮服务器的问题：单线程、每次请求间隔、失败指数退避
输出: data/circbank_human.jsonl（追加）
"""
import urllib.request, json, ssl, time, os, sys

ctx = ssl.create_default_context()
H = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) research-tool/1.0',
     'Content-Type': 'application/json'}
OUT = 'data/circbank_human.jsonl'
PS = 500                # 小页，降低单次查询压力
TOTAL = 2033961
MAXPAGE = TOTAL // PS + 1
GAP = 1.2# 每次请求间隔（秒），礼貌抓取


def post(page, tries=2, timeout=60):
    d = json.dumps({'currentPage': page, 'pageSize': PS,
                    'keyword': '', 'circBankID': ''}).encode()
    for a in range(tries):
        try:
            return json.loads(urllib.request.urlopen(
                urllib.request.Request('http://www.circbank.cn/search', data=d, headers=H),
                timeout=timeout, context=ctx).read().decode())
        except Exception:
            if a == tries - 1:
                return None
            time.sleep(5)


seen = set()
if os.path.exists(OUT):
    with open(OUT, encoding='utf-8') as f:
        for l in f:
            try:
                seen.add(json.loads(l)['circBankID'])
            except Exception:
                pass
print('[INFO] 已有 %d 条' % len(seen), flush=True)

# 先探测服务器是否可用
probe = post(1)
if probe is None:
    print('[WAIT] 服务器仍不可用，进入等待重试循环', flush=True)
    wait = 60
    while probe is None and wait <= 3600:
        time.sleep(wait)
        probe = post(1)
        wait *= 1.5
        print('[WAIT] 重试中... 累计等待 %.0f 秒' % (wait,), flush=True)
if probe is None:
    print('[GIVEUP] 服务器长时间不可用，退出', flush=True)
    sys.exit(1)

print('[OK] 服务器恢复，开始抓取', flush=True)

page = 1
hsa = 0
batch = []
t0 = time.time()
stall = 0

while page <= MAXPAGE and hsa < 200000:
    r = post(page)
    if r is None:
        stall += 1
        print('[STALL] 第 %d 页失败，连续 %d 次，休眠 30 秒' % (page, stall), flush=True)
        time.sleep(30)
        if stall >= 20:
            print('[GIVEUP] 连续失败过多，退出（已保存 %d 条）' % len(seen), flush=True)
            break
        continue
    stall = 0
    for x in r.get('list') or []:
        cid = x.get('circBankID') or ''
        if not cid.startswith('hsa_') or cid in seen:
            continue
        seen.add(cid)
        batch.append(json.dumps(x, ensure_ascii=False))
        hsa += 1
    if len(batch) >= 500:
        with open(OUT, 'a', encoding='utf-8') as f:
            f.write('\n'.join(batch) + '\n')
        batch = []
    page += 1
    time.sleep(GAP)
    if page % 50 == 0:
        el = time.time() - t0
        print('  第%d页 | 人源 %d 条 | %.0f 秒 | %.1f 页/秒'
              % (page, len(seen), el, page / max(el, 1)), flush=True)

if batch:
    with open(OUT, 'a', encoding='utf-8') as f:
        f.write('\n'.join(batch) + '\n')

print('[DONE] 共 %d 条，耗时 %.0f 秒' % (len(seen), time.time() - t0), flush=True)