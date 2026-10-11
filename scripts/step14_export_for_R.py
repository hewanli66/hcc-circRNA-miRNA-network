#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""step14_export_for_R.py

把 TCGA-LIHC 与 GSE14520 的表达、临床协变量、生存时间导出为 R 可读 TSV，
供 step15_multivar_cox.R 做多因素 Cox 与数值基准验证。

输出（默认 results/r_input/）：
  tcga_expr.tsv     基因 × 样本 的表达矩阵（Ensembl 原始行名保留，附 symbol 列）
  tcga_clin.tsv     sample, stage, age, male, time, event
  gse14520_expr.tsv 同上（Affy probe → symbol 已映射）
  gse14520_clin.tsv 同上
  gene_list.tsv     本次分析的宿主基因清单
"""
import csv
import gzip
import math
import os
import re
import sys

RAW = 'C:/Users/hewan/Documents/kimi/tasks/2026-10-08/11-38-08-f56eebff/rawdata'
HCC = 'C:/Users/hewan/WorkBuddy/2026-10-09-22-01-46/hcc'
OUT = os.path.join(HCC, 'results', 'r_input')
os.makedirs(OUT, exist_ok=True)

# 待分析的 6 个跨队列验证基因 + 若干对照
GENES = ['NT5DC2', 'HSD17B6', 'BCLAF1', 'SREBF2', 'CD5', 'LRCH1',
         'RBM19', 'SMG5', 'NUP205']


# ==================== TCGA ====================
def load_ensembl_symbol():
    """gencode.v22 注释：Ensembl ID → gene symbol"""
    m = {}
    p = os.path.join(RAW, 'annotation_cache', 'gencode.v22.annotation.gtf.gz')
    with gzip.open(p, 'rt', encoding='utf-8', errors='ignore') as f:
        for l in f:
            if l.startswith('#'):
                continue
            mg = re.search(r'gene_id "([^"]+)"', l)
            mn = re.search(r'gene_name "([^"]+)"', l)
            if mg and mn:
                m.setdefault(mg.group(1).split('.')[0], mn.group(1))
    return m


def export_tcga():
    id2sym = load_ensembl_symbol()
    expr_path = os.path.join(OUT, 'tcga_expr.tsv')
    want = set(GENES)
    n_found = 0
    with gzip.open(os.path.join(RAW, 'TCGA-LIHC.star_fpkm.tsv.gz'),
                   'rt', encoding='utf-8', errors='ignore') as f, \
         open(expr_path, 'w', encoding='utf-8', newline='') as out:
        hdr = f.readline().rstrip('\n').split('\t')
        samples = hdr[1:]
        w = csv.writer(out, delimiter='\t', lineterminator='\n')
        w.writerow(['gene_symbol'] + samples)
        for l in f:
            p = l.rstrip('\n').split('\t')
            sym = id2sym.get(p[0].split('.')[0])
            if sym in want:
                w.writerow([sym] + p[1:])
                n_found += 1

    # 临床 + 生存
    with gzip.open(os.path.join(RAW, 'TCGA-LIHC.clinical.tsv.gz'),
                   'rt', encoding='utf-8', errors='ignore') as f:
        clin = {r['sample']: r for r in csv.DictReader(f, delimiter='\t')}
    with gzip.open(os.path.join(RAW, 'TCGA-LIHC.survival.tsv.gz'),
                   'rt', encoding='utf-8', errors='ignore') as f:
        surv = {r['sample']: (float(r['OS.time']), int(r['OS']))
                for r in csv.DictReader(f, delimiter='\t')}

    def stage_of(s):
        v = (clin[s].get('ajcc_pathologic_stage.diagnoses') or '').strip().upper()
        if v in ('STAGE I', 'STAGE IA', 'STAGE IB'):
            return 0
        if v.startswith(('STAGE II', 'STAGE III', 'STAGE IV')):
            return 1
        return None

    n = 0
    with open(os.path.join(OUT, 'tcga_clin.tsv'), 'w',
              encoding='utf-8', newline='') as out:
        w = csv.writer(out, delimiter='\t', lineterminator='\n')
        w.writerow(['sample', 'stage', 'age', 'male', 'time', 'event'])
        for s, (t, e) in sorted(surv.items()):
            if s not in clin:
                continue
            st = stage_of(s)
            if st is None:
                continue
            try:
                age = float(clin[s].get('age_at_index.demographic', ''))
            except ValueError:
                continue
            g = (clin[s].get('gender.demographic') or '').strip().lower()
            if g not in ('male', 'female'):
                continue
            w.writerow([s, st, age, 1 if g == 'male' else 0, t, e])
            n += 1
    print('TCGA  表达 %d 个基因 | 临床 %d 例' % (n_found, n))
    return n_found, n


# ==================== GSE14520 ====================
def export_gse14520():
    expr_src = os.path.join(RAW, 'GSE14520-GPL3921_series_matrix.txt.gz')
    map_src = os.path.join(RAW, 'annotation_cache', 'gpl3921_probe_symbol.csv')
    if not os.path.exists(expr_src):
        print('GSE14520 表达矩阵未找到，跳过', file=sys.stderr)
        return 0, 0

    # Affy probe → gene symbol
    p2s = {}
    with open(map_src, encoding='utf-8', errors='ignore') as f:
        for row in csv.reader(f):
            if len(row) >= 2 and row[0] and not row[0].startswith('Unnamed'):
                p2s[row[0].strip()] = row[1].strip()

    want = set(GENES)
    rows_by_gene = {}
    with gzip.open(expr_src, 'rt', encoding='utf-8', errors='ignore') as f:
        hdr = None
        seen_begin = False
        for l in f:
            ls = l.strip()
            if ls.startswith('!series_matrix_table_begin'):
                seen_begin = True
                continue
            if ls.startswith('!series_matrix_table_end'):
                break
            if not seen_begin:
                continue
            if not ls:
                continue
            p = [x.strip('"') for x in l.rstrip('\n').split('\t')]
            if hdr is None:
                hdr = p[1:]
                continue
            sym = p2s.get(p[0], '')
            if sym not in want:
                continue
            vals = []
            for v in p[1:]:
                try:
                    vals.append(float(v))
                except ValueError:
                    vals.append(float('nan'))
            rows_by_gene.setdefault(sym, []).append(vals)
    assert hdr is not None and len(hdr) > 100, 'GSE14520 表头未正确解析'

    # 同一基因多个 probe → 取跨样本平均，保证每个基因只有一行
    n_found = 0
    with open(os.path.join(OUT, 'gse14520_expr.tsv'), 'w',
              encoding='utf-8', newline='') as out:
        wr = csv.writer(out, delimiter='\t', lineterminator='\n')
        wr.writerow(['gene_symbol'] + hdr)
        for sym in GENES:
            probes = rows_by_gene.get(sym)
            if not probes:
                continue
            k = len(probes)
            ncol = len(probes[0])
            mean = [sum(probes[j][i] for j in range(k)) / k for i in range(ncol)]
            wr.writerow([sym] + ['%.4f' % v for v in mean])
            n_found += 1

    # 临床 + 生存
    with gzip.open(os.path.join(RAW, 'GSE14520_Extra_Supplement.txt.gz'),
                   'rt', encoding='utf-8', errors='ignore') as f:
        lines = [l.rstrip('\n').split('\t') for l in f if l.strip()]
    hdr2 = lines[0]
    idx = {h: i for i, h in enumerate(hdr2)}
    n = 0
    with open(os.path.join(OUT, 'gse14520_clin.tsv'), 'w',
              encoding='utf-8', newline='') as out:
        w = csv.writer(out, delimiter='\t', lineterminator='\n')
        w.writerow(['sample', 'stage', 'age', 'male', 'time', 'event'])
        for row in lines[1:]:
            gsm = row[idx['Affy_GSM']].strip()
            tt = row[idx['Tissue Type']].strip().lower()
            if 'tumor' not in tt:
                continue
            tnm = row[idx['TNM staging']].strip()
            # TNM 形如 I / II / III / IVA/IVB
            m = re.match(r'\s*(I{1,3}V?|IV[AB]?)', tnm.upper())
            if not m:
                continue
            s = m.group(1)
            stage = 0 if s in ('I', 'II') else 1
            try:
                age = float(row[idx['Age']])
            except ValueError:
                continue
            #🔴 Gender 列在 GSE14520 Extra Supplement 中存的是 **M / F** 单字母，
            #   不是 "male" / "female"。原代码写的是
            #       == 'male'
            #   永远为 False，导致全部 225 例被编码成 male=0，R 里 sex 无方差、
            #   被自动剔除，还让稿件写出「所有肿瘤样本均为女性」这一明显错误结论
            #   （真实分布为 211 男 / 31 女，男占 85%，符合 HBV 相关 HCC 特征）。
            #   审稿人第二轮 Major comment 4 正是质疑这一点。
            #   这里按单字母判定，并显式排除空值。
            g_raw = row[idx['Gender']].strip().upper()
            if g_raw not in ('M', 'F'):
                continue          # 缺失性别不应默认为女性
            male = 1 if g_raw == 'M' else 0
            #生存状态缺失或非数值的排除（与 step3_validate14520 口径一致）
            ss = row[idx['Survival status']].strip()
            if ss in ('', 'NA', '.'):
                continue
            try:
                t = float(row[idx['Survival months']])
            except ValueError:
                continue
            if t <= 0:      # 与 step3 一致：排除零生存时间
                continue
            t *= 30.4375   # 月 → 天
            try:
                e = int(float(ss))
            except ValueError:
                continue
            w.writerow([gsm, stage, age, male, '%.4f' % t, e])
            n += 1
    missing = sorted(want - set(rows_by_gene))
    print('GSE14520  表达 %d 个基因 | 临床 %d 例' % (n_found, n))
    if missing:
        print('  未找到的基因:', missing)
    return n_found, n




def export_raw_expr():
    """额外导出原始（未标准化）表达矩阵，供R 做与稿件同口径的基准验证。
    TCGA 用原始 FPKM，GSE14520 用原始探针强度（多probe 取均值）。"""
    # ---- TCGA ----
    id2sym = load_ensembl_symbol()
    want = set(GENES)
    n = 0
    with gzip.open(os.path.join(RAW, 'TCGA-LIHC.star_fpkm.tsv.gz'),
                   'rt', encoding='utf-8', errors='ignore') as f, \
         open(os.path.join(OUT, 'tcga_expr_raw.tsv'), 'w',
              encoding='utf-8', newline='') as out:
        hdr = f.readline().rstrip('\n').split('\t')
        w = csv.writer(out, delimiter='\t', lineterminator='\n')
        w.writerow(['gene_symbol'] + hdr[1:])
        for l in f:
            q = l.rstrip('\n').split('\t')
            sym = id2sym.get(q[0].split('.')[0])
            if sym in want:
                w.writerow([sym] + q[1:])
                n += 1
    print('tcga_expr_raw.tsv  %d 个基因' % n)

    # ---- GSE14520 ----
    p2s = {}
    with open(os.path.join(RAW, 'annotation_cache', 'gpl3921_probe_symbol.csv'),
              encoding='utf-8', errors='ignore') as f:
        for row in csv.reader(f):
            if len(row) >= 2 and row[0] and not row[0].startswith('Unnamed'):
                p2s[row[0].strip()] = row[1].strip()
    by = {}
    hdr = None
    seen = False
    with gzip.open(os.path.join(RAW, 'GSE14520-GPL3921_series_matrix.txt.gz'),
                   'rt', encoding='utf-8', errors='ignore') as f:
        for l in f:
            ls = l.strip()
            if ls.startswith('!series_matrix_table_begin'):
                seen = True
                continue
            if ls.startswith('!series_matrix_table_end'):
                break
            if not seen or not ls:
                continue
            q = [x.strip('"') for x in l.rstrip('\n').split('\t')]
            if hdr is None:
                hdr = q[1:]
                continue
            sym = p2s.get(q[0], '')
            if sym not in want:
                continue
            vals = []
            for v in q[1:]:
                try:
                    vals.append(float(v))
                except ValueError:
                    vals.append(float('nan'))
            by.setdefault(sym, []).append(vals)
    m = 0
    with open(os.path.join(OUT, 'gse14520_expr_raw.tsv'), 'w',
              encoding='utf-8', newline='') as out:
        w = csv.writer(out, delimiter='\t', lineterminator='\n')
        w.writerow(['gene_symbol'] + hdr)
        for sym in GENES:
            probes = by.get(sym)
            if not probes:
                continue
            k = len(probes)
            mean = [sum(probes[j][i] for j in range(k)) / k
                    for i in range(len(probes[0]))]
            w.writerow([sym] + ['%.6f' % v for v in mean])
            m += 1
    print('gse14520_expr_raw.tsv  %d 个基因' % m)


if __name__ == '__main__':
    print('输出目录:', OUT)
    export_tcga()
    export_gse14520()
    export_raw_expr()
    with open(os.path.join(OUT, 'gene_list.tsv'), 'w',
              encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['gene'])
        for g in GENES:
            w.writerow([g])
    print('完成')
