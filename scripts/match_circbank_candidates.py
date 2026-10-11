# -*- coding: utf-8 -*-
"""
GSE122482 候选 circRNA × circBank 库匹配（按 position 前缀）

输入:
  results/candidate_circRNA.tsv    103 条候选, 坐标格式 hsa_circ:chrN:start-end
  data/circbank_human.jsonl        circBank 抓取结果（每行一个 JSON）
输出:
  results/circbank_match_detail.tsv逐候选匹配明细
  results/circbank_match_summary.json汇总
匹配规则:
  circBank 的 position字段形如 "chr12:56754955-56755968|+"（坐标 + 链方向）
  候选坐标剥掉 "hsa_circ:" 前缀后得到 "chr12:56754955-56755968"
  → 前缀匹配；strict 为「余下部分恰为 |+ 或 |-」，可排除 100-2000 误配 100-200 的情况
"""
import json, os, re, sys
from collections import defaultdict

CAND = 'results/candidate_circRNA.tsv'
CB = 'data/circbank_human.jsonl'
OUT_DETAIL = 'results/circbank_match_detail.tsv'
OUT_SUM = 'results/circbank_match_summary.json'


def norm_pos(p):
    """circBank position -> (chr, start, end, strand)"""
    m = re.match(r'^(chr[0-9XYM]+):(\d+)-(\d+)\|([+-])$', (p or '').strip())
    if not m:
        return None
    return m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)


def bare(cid):
    """hsa_circ:chrN:start-end -> chrN:start-end"""
    return cid.split(':', 1)[1] if cid.startswith('hsa_circ:') else cid


# ---------- 读候选 ----------
cands = []
with open(CAND, encoding='utf-8') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    for line in f:
        line = line.rstrip('\n')
        if not line.strip():
            continue
        p = line.split('\t')
        cands.append(dict(zip(hdr, p)))
print('[CAND] %d 条候选' % len(cands))

# ---------- 建索引：完整 position 串（含链方向）-> [circBank 记录] ----------
idx = defaultdict(list)
total = 0
bad = 0
nopos = 0
with open(CB, encoding='utf-8') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        try:
            x = json.loads(line)
        except Exception:
            bad += 1
            continue
        total += 1
        pp = norm_pos(x.get('position'))
        if pp is None:
            nopos += 1
            continue
        chrom, st, en, strand = pp
        # 用还原后的规范串做键，保证 start-end 边界不被前缀误吃
        idx['%s:%d-%d|%s' % (chrom, st, en, strand)].append(x)

print('[CB] 有效记录 %d 条（坏行 %d，坐标格式不符%d），唯一 position 键 %d 个'
      % (total, bad, nopos, len(idx)))

# ---------- 匹配 ----------
# circBank position = "chr12:56754955-56755968|+"，候选坐标 = "chr12:56754955-56755968"
#前缀匹配后再要求余下部分恰为 "|+" / "|-"，否则 "chr1:100-200" 会误吃 "chr1:100-2000"
rows = []
hit_pref = 0
hit_strict = 0
hit_ambig = 0
for cd in cands:
    coord = bare(cd['circRNA'])
    # 前缀匹配：遍历以 coord 开头的键（键数远小于候选数，成本可忽略）
    pref = []
    for k, v in idx.items():
        if k.startswith(coord):
            pref.extend(v)
    if pref:
        hit_pref += 1
    for x in pref:
        cbpos = x.get('position', '')
        rem = cbpos[len(coord):]
        strict = rem in ('|+', '|-')
        if strict:
            hit_strict += 1
        else:
            hit_ambig += 1
        gene = cd.get('gene', '')
        rows.append(dict(
            circRNA=cd['circRNA'], gene=gene,
            log2FC=cd.get('log2FC', ''), P_adj=cd.get('P_adj', ''),
            regulation=cd.get('regulation', ''),
            circBankID=x.get('circBankID', ''), cb_position=cbpos,
            cb_gene=x.get('geneSymbol', ''), cb_geneID=x.get('geneID', ''),
            circType=x.get('circType', ''), splicedLength=x.get('splicedLength', ''),
            samples=x.get('samples', ''), miRNA=x.get('miRNA', ''),
            m6A=x.get('m6A', ''), ires=x.get('ires', ''), orf=x.get('orf', ''),
            conserved=x.get('conserved', ''), validated=x.get('validated', ''),
            source=x.get('source', ''),
            strict='1' if strict else '0',
        ))
# 同一条候选可能因 + / - 两条链各命中一次，按候选去重计数
uniq_pref = len({r['circRNA'] for r in rows})
uniq_strict = len({r['circRNA'] for r in rows if r['strict'] == '1'})
uniq_ambig = len({r['circRNA'] for r in rows if r['strict'] == '0'})

# ---------- 输出 ----------
os.makedirs('results', exist_ok=True)
cols = ['circRNA', 'gene', 'log2FC', 'P_adj', 'regulation', 'circBankID',
        'cb_position', 'cb_gene', 'cb_geneID', 'circType', 'splicedLength',
        'samples', 'miRNA', 'm6A', 'ires', 'orf', 'conserved', 'validated',
        'source', 'strict']
with open(OUT_DETAIL, 'w', encoding='utf-8', newline='') as f:
    f.write('\t'.join(cols) + '\n')
    for r in rows:
        f.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')

# 未命中的候选也要留痕
hit_ids = {r['circRNA'] for r in rows}
with open(OUT_DETAIL, 'a', encoding='utf-8') as f:
    for cd in cands:
        if cd['circRNA'] not in hit_ids:
            f.write('\t'.join([cd['circRNA'], cd.get('gene', ''),
                               cd.get('log2FC', ''), cd.get('P_adj', ''),
                               cd.get('regulation', '')] + [''] * 14 + ['0']) + '\n')

gene_ok = len({r['circRNA'] for r in rows
               if r['cb_gene'] and r['gene'] and
               r['cb_gene'].upper() == r['gene'].upper()})

summary = dict(
    cb_records_total=total,
    cb_unique_position_keys=len(idx),
    candidates=len(cands),
    matched_prefix=uniq_pref,
    matched_strict=uniq_strict,
    unmatched=len(cands) - uniq_pref,
    match_rate_prefix=round(uniq_pref / len(cands) * 100, 2),
    match_rate_strict=round(uniq_strict / len(cands) * 100, 2),
    matched_and_gene_symbol_same=gene_ok,
    rows_detail=len(rows),
    cb_total_records_in_db=2033961,
    cb_coverage_pct=round(total / 2033961 * 100, 2),
)
json.dump(summary, open(OUT_SUM, 'w', encoding='utf-8'),
          ensure_ascii=False, indent=2)

print('\n=====匹配结果=====')
for k, v in summary.items():
    print('  %-32s %s' % (k, v))
print('\n[SAVED] %s (%d 行)' % (OUT_DETAIL, len(rows) + 1))
print('[SAVED] %s' % OUT_SUM)