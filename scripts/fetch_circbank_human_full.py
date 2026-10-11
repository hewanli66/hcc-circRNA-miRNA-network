# -*- coding: utf-8 -*-
"""
circBank 全量人源 circRNA 温和抓取（单线程 + 限速 + 断点续传）
策略要点（针对此前压垮服务器的问题）：
  1) 单线程，绝不并发
  2) pageSize=500，间隔 1.0 秒
  3) 失败指数退避，最多重试 40 次后放弃（可安全重跑，续传）
  4) 只保留 circBankID 以 hsa_ 开头的条目，其余丢弃（不写盘）

输出: data/circbank_human.jsonl（增量追加，约 14 万条）
      data/circbank_meta.json（抓取元信息）
"""
import urllib.request, json, ssl, time, os, sys

ctx = ssl.create_default_context()
H = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) research-tool/1.0',
     'Content-Type': 'application/json'}
OUT = 'data/circbank_human.jsonl'
META = 'data/circbank_meta.json'
PS = 500
TOTAL = 2033961
MAXPAGE = TOTAL // PS + 1
GAP = 1.5
FLUSH_EVERY = 300          # 每累积 300 条写盘一次


def post(page, tries=2, timeout=90):
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
            time.sleep(6)


# ---------- 断点续传：读已有 ----------
seen = set()
if os.path.exists(OUT):
    with open(OUT, encoding='utf-8') as f:
        for l in f:
            try:
                seen.add(json.loads(l)['circBankID'])
            except Exception:
                pass
print('[INFO] 本地已有 %d 条' % len(seen), flush=True)

# ---------- 启动探测 ----------
probe = post(1)
if probe is None:
    print('[ERROR] 服务器仍不可用，退出（稍后重跑即可续传）', flush=True)
    sys.exit(2)
print('[OK] 服务器可用，开始抓取', flush=True)

batch = []
fail = 0
page = 1
t0 = time.time()
hsa_new = 0

while page <= MAXPAGE:
    r = post(page)
    if r is None:
        fail += 1
        print('[STALL] 第 %d 页失败（连续 %d 次），退避 20 秒' % (page, fail), flush=True)
        time.sleep(20)
        if fail >= 40:
            print('[GIVEUP] 连续失败过多，安全退出（可重跑续传）', flush=True)
            break
        continue
    fail = 0
    for x in r.get('list') or []:
        cid = x.get('circBankID') or ''
        if not cid.startswith('hsa_') or cid in seen:
            continue
        seen.add(cid)
        batch.append(json.dumps(x, ensure_ascii=False))
        hsa_new += 1
        if len(batch) >= FLUSH_EVERY:
            with open(OUT, 'a', encoding='utf-8') as f:
                f.write('\n'.join(batch) + '\n')
            batch = []
    page += 1
    time.sleep(GAP)
    if page % 100 == 0:
        el = time.time() - t0
        rate = (page - 1) / max(el, 1)
        eta = (MAXPAGE - page) / max(rate, 0.001)
        print('  %d/%d 页 | 人源 %d 条(累计%d) | %.2f 页/秒 | 剩余 ~%.0f 分钟'
              % (page, MAXPAGE, hsa_new, len(seen), rate, eta / 60), flush=True)

if batch:
    with open(OUT, 'a', encoding='utf-8') as f:
        f.write('\n'.join(batch) + '\n')

el = time.time() - t0
meta = dict(total_site=TOTAL, page_size=PS, pages_fetched=page - 1,
            human_records=len(seen), new_this_run=hsa_new,
            elapsed_sec=round(el, 1), finished=page > MAXPAGE,
            source='http://www.circbank.cn/search')
json.dump(meta, open(META, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)

print('[DONE] 人源共 %d 条，本次新增 %d，耗时 %.0f 秒，完成度 %d/%d 页'
      % (len(seen), hsa_new, el, page - 1, MAXPAGE), flush=True)
print('[SAVED] %s' % OUT, flush=True)