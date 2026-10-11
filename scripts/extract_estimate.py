# -*- coding: utf-8 -*-
"""提取 ESTIMATE 评分并检查 HCC 覆盖"""
import zipfile, re, os, csv
from collections import Counter

RAW = 'C:/Users/hewan/Documents/kimi/tasks/2026-10-08/11-38-08-f56eebff/rawdata'
OUT = 'results'


def col_to_idx(ref):
    m = re.match(r'([A-Z]+)', ref)
    if not m: return -1
    n = 0
    for ch in m.group(1):
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def read_all(path, sheet_idx):
    z = zipfile.ZipFile(path)
    shared = []
    if 'xl/sharedStrings.xml' in z.namelist():
        s = z.read('xl/sharedStrings.xml').decode('utf-8', 'ignore')
        for si in re.findall(r'<si>(.*?)</si>', s, re.S):
            shared.append(''.join(re.findall(r'<t[^>]*>(.*?)</t>', si, re.S)))
    names = sorted(n for n in z.namelist() if n.startswith('xl/worksheets/sheet'))
    fn = names[sheet_idx]
    sh = z.read(fn).decode('utf-8', 'ignore')
    rows = re.findall(r'<row[^>]*>(.*?)</row>', sh, re.S)
    out = []
    for r in rows:
        cells = {}
        for m in re.finditer(r'<c([^>]*)>(.*?)</c>', r, re.S):
            attrs, body = m.group(1), m.group(2)
            t = re.search(r'\st="([^"]+)"', attrs)
            t = t.group(1) if t else 'n'
            ref = re.search(r'r="([A-Z]+\d+)"', attrs)
            ci = col_to_idx(ref.group(1)) if ref else len(cells)
            v = re.search(r'<v>(.*?)</v>', body, re.S)
            isv = re.search(r'<is>(.*?)</is>', body, re.S)
            if t == 's' and v:
                try: val = shared[int(v.group(1))]
                except (ValueError, IndexError): val = ''
            elif isv:
                val = ''.join(re.findall(r'<t[^>]*>(.*?)</t>', isv.group(1), re.S))
            else:
                val = v.group(1) if v else ''
            cells[ci] = val
        out.append([cells.get(i, '') for i in range(max(cells) + 1 if cells else 0)])
    return out


p = os.path.join(RAW, 'ESTIMATE_moesm489.xlsx')
allrows = []
for si in range(3):
    rows = read_all(p, si)
    if not rows: continue
    hdr_idx = None
    for i, r in enumerate(rows[:5]):
        joined = ' '.join(r[:8])
        if 'Stromal score' in joined and 'Immune score' in joined:
            hdr_idx = i
            break
    if hdr_idx is None:
        print('sheet%d 未找到表头' % (si + 1))
        continue
    hdr = rows[hdr_idx]
    # 用表头名称定位列索引，避免固定索引错位
    def find(name):
        for j, h in enumerate(hdr):
            if name in str(h):
                return j
        return None
    iID = find('ID')
    iDis = find('Disease')
    iPlat = find('Platform')
    iStr = find('Stromal')
    iImm = find('Immune')
    iEst = find('ESTIMATE')
    iPur = find('purity')
    print('sheet%d 列定位: ID=%s Disease=%s Platform=%s Stromal=%s Immune=%s EST=%s purity=%s'
          % (si + 1, iID, iDis, iPlat, iStr, iImm, iEst, iPur))
    body = rows[hdr_idx + 1:]
    need = max(x for x in (iID, iDis, iPlat, iStr, iImm, iEst, iPur) if x is not None)
    good = []
    for r in body:
        if len(r) <= need: continue
        if not r[iID].strip() or not r[iStr].strip(): continue
        good.append([r[iID].strip(), r[iDis].strip() if iDis is not None else '',
                     r[iPlat].strip() if iPlat is not None else '',
                     r[iStr].strip(), r[iImm].strip() if iImm is not None else '',
                     r[iEst].strip() if iEst is not None else '',
                     r[iPur].strip() if iPur is not None else ''])
    print('  有效行 %d，样例: %s' % (len(good), good[0] if good else None))
    allrows += good

print('\n合并总行: %d' % len(allrows))
dis = Counter(r[1] for r in allrows)
print('\n疾病分布:')
for k, v in dis.most_common(20):
    print('  %-30s %d' % (k, v))

liver = [r for r in allrows
         if 'liver' in r[1].lower() or 'LIHC' in r[0].upper() or 'hepat' in r[1].lower()]
print('\n肝相关条目: %d' % len(liver))
for r in liver[:6]:
    print('    ID=%s | Disease=%s | Platform=%s | Stromal=%s Immune=%s EST=%s purity=%s'
          % (r[0], r[1], r[2], r[3][:8], r[4][:8], r[5][:8], r[6]))

os.makedirs(OUT, exist_ok=True)
with open(f'{OUT}/estimate_all.tsv', 'w', encoding='utf-8', newline='') as fo:
    w = csv.writer(fo, delimiter='\t')
    w.writerow(['ID', 'Disease', 'Platform', 'Stromal', 'Immune', 'ESTIMATE', 'purity'])
    for r in allrows:
        w.writerow(r)
print('\n[SAVED] %s/estimate_all.tsv' % OUT)