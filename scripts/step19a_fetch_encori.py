# -*- coding: utf-8 -*-
"""step19a_fetch_encori.py — 从 ENCORI API 下载指定基因的 miRNA 靶点预测

ENCORI(starBase v3 后继) 整合 7 个独立算法（PITA/RNA22/miRmap/microT/
miRanda/PicTar/TargetScan）并整合 CLIP-seq / degradome-seq 实验支持。
每个基因分别用全部 7 个算法查询，避免 programNum 截断。

用法：python step19a_fetch_encori.py NT5DC2 HSD17B6 ...
"""
import json
import os
import ssl
import sys
import time
import urllib.parse
import urllib.request

HCC = 'C:/Users/hewan/WorkBuddy/2026-10-09-22-01-46/hcc'
OUT = os.path.join(HCC, 'data', 'encori')
os.makedirs(OUT, exist_ok=True)

API = 'https://rnasysu.com/encori/api/miRNATarget/'
PROGRAMS = ['PITA', 'RNA22', 'miRmap', 'microT', 'miRanda',
            'PicTar', 'TargetScan']

CTX = ssl.create_default_context()
HDR = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) research-tool/1.0'}


def fetch_one(gene, program, tries=3, timeout=120):
    """ENCORI 的 program 参数一次只接受一个算法（实测：多算法会报
    'program parameter haven't been set correctly'）。故逐算法查询。"""
    params = {
        'assembly': 'hg38',
        'geneType': 'mRNA',
        'target': gene,
        'clipExpNum': '0',
        'degraExpNum': '0',
        'pancancerNum': '0',
        'programNum': '1',
        'program': program,
    }
    url = API + '?' + urllib.parse.urlencode(params)
    for a in range(tries):
        try:
            req = urllib.request.Request(url, headers=HDR)
            txt = urllib.request.urlopen(req, timeout=timeout,
                                         context=CTX).read().decode()
            lines = [l for l in txt.split('\n') if l and not l.startswith('#')]
            if len(lines) < 2:
                return []
            return [l.rstrip('\r') for l in lines]
        except Exception as ex:
            print('    %s attempt %d failed: %s' % (program, a + 1, str(ex)[:70]))
            time.sleep(3)
    return None


def fetch(gene):
    """合并 7 个算法的结果，按 miRNA 去重（同一 miRNA 记它在哪些算法被预测）"""
    header = None
    rows = {}
    for prog in PROGRAMS:
        lines = fetch_one(gene, prog)
        if lines is None:
            print('    %-10s %s : failed' % (gene, prog))
            continue
        if not lines:
            print('    %-10s %-10s : 0' % (gene, prog))
            continue
        if header is None:
            header = lines[0].split('\t')
            hidx = {c: i for i, c in enumerate(header)}
        n_ok = 0
        for l in lines[1:]:
            c = l.split('\t')
            if len(c) < len(header):
                continue
            mi = c[1]
            if mi not in rows:
                rows[mi] = dict((p, '0') for p in PROGRAMS)
                rows[mi]['miRNAname'] = mi
                rows[mi]['geneName'] = c[3]
                rows[mi]['geneID'] = c[2]
                rows[mi]['chromosome'] = c[5]
                rows[mi]['narrowStart'] = c[6]
                rows[mi]['narrowEnd'] = c[7]
                rows[mi]['clipExpNum'] = c[11]
                rows[mi]['degraExpNum'] = c[12]
            try:
                rows[mi][prog] = c[hidx[prog]]
            except (KeyError, IndexError):
                pass
            n_ok += 1
        print('    %-10s %-10s : %d' % (gene, prog, n_ok))
    if header is None or not rows:
        return 0
    cols = ['miRNAname', 'geneID', 'geneName', 'chromosome',
            'narrowStart', 'narrowEnd', 'clipExpNum', 'degraExpNum'] + PROGRAMS
    path = os.path.join(OUT, gene + '.tsv')
    with open(path, 'w', encoding='utf-8', newline='') as f:
        f.write('\t'.join(cols) + '\n')
        for mi in sorted(rows):
            f.write('\t'.join(str(rows[mi].get(c, '')) for c in cols) + '\n')
    return len(rows)


if __name__ == '__main__':
    genes = sys.argv[1:] or ['NT5DC2', 'HSD17B6', 'BCLAF1', 'SREBF2',
                             'CD5', 'LRCH1', 'ADAMTS13']
    print('ENCORI API: %s' % API)
    print('algorithms: %s' % ', '.join(PROGRAMS))
    print()
    summary = {}
    for g in genes:
        n = fetch(g)
        if not n:
            print('  %-10s -> 0 miRNA' % g)
            summary[g] = 0
        else:
            print('  %-10s -> %d miRNA (7 algorithms merged)' % (g, n))
            summary[g] = n
    with open(os.path.join(OUT, '_fetch_summary.json'), 'w',
              encoding='utf-8') as f:
        json.dump(summary, f, indent=2)
    print('\ndone: %s' % os.path.join(OUT, '_fetch_summary.json'))
