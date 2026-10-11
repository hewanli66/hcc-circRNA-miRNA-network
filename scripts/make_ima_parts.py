# -*- coding: utf-8 -*-
"""
把 circBank 人源全量条目转成 TSV 并分卷，供上传到 ima 知识库。

为什么这么做：
  原始 data/circbank_human.jsonl 是 464 MB / 828,808 条 JSON 记录，
  既超出 ima 知识库的单文件上限（TXT/Markdown 10 MB，PDF/Word 200 MB），
  扩展名 .jsonl 也不在支持列表里（preflight 的 EXT_MAP 只认
  pdf/doc/docx/ppt/pptx/xls/xlsx/csv/md/txt/xmind/图片/音频/html/epub），
  zip 与 gz 均不支持，只能传纯文本。

做法：
  1. JSON → TSV，去掉 JSON 包装与空值冗余，体积从 464 MB 降到约 85 MB
  2. 按每卷 ≤ MAX_BYTES 切分，命名 part_01.txt ...
  3. 第一卷带表头，后续卷不重复表头；每卷首行注明它是第几卷/共几卷

保留字段（去掉了 aliasID1-3、geneID、bestTranscriptID、type/typeName 等
在本次分析中未使用的列；如需完整字段可改 COLS）：
  circBankID position seqnames circType splicedLength geneSymbol
  biotype tissue samples miRNA m6A codingProbability conserved
  validated highConfidence ires orf
注意：conserved / validated / highConfidence 三列在 circBank 人源条目中
全部为空，仍予保留以忠实反映源数据（这一点在稿件中已声明）。
"""
import json, os, math

SRC = 'data/circbank_human.jsonl'
OUTDIR = 'ima_parts'
MAX_BYTES = 9 * 1024 * 1024        # 每卷 9 MB，留出余量避免超 10 MB 上限
PREFIX = 'circbank_human_tsv'

COLS = ['circBankID', 'position', 'seqnames', 'circType', 'splicedLength',
        'geneSymbol', 'biotype', 'tissue', 'samples', 'miRNA', 'm6A',
        'codingProbability', 'conserved', 'validated', 'highConfidence',
        'ires', 'orf']


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    total = sum(1 for _ in open(SRC, encoding='utf-8'))
    # 第一卷要放表头+卷说明，先按空卷估个数再修正
    n_parts = math.ceil((total * 120) / MAX_BYTES) + 1   # 120 B/行 粗估
    print('源文件 %s，%d 条记录' % (SRC, total))
    print('预估 %d 卷，每卷 ≤ %.1f MB' % (n_parts, MAX_BYTES / 1048576))

    written = []
    part_no = 0
    rows_in_part = 0
    fh = None
    path = None

    def open_new():
        nonlocal fh, path, part_no, rows_in_part
        part_no += 1
        path = os.path.join(OUTDIR, '%s_part%02d.txt' % (PREFIX, part_no))
        fh = open(path, 'w', encoding='utf-8', newline='\n')
        # 卷说明行 + 表头（仅第一卷）
        fh.write('# circBank human entries (TSV), part %d\n' % part_no)
        if part_no == 1:
            fh.write('\t'.join(COLS) + '\n')
        rows_in_part = 0
        written.append(path)

    def close_part():
        nonlocal fh
        if fh:
            fh.close()
            fh = None

    open_new()
    with open(SRC, encoding='utf-8') as f:
        for line in f:
            d = json.loads(line)
            row = '\t'.join(('' if d.get(c) is None else str(d.get(c)))
                            for c in COLS)
            # +len(row) 估算写入后字节数（含换行）
            if fh.tell() + len(row.encode('utf-8')) + 1 > MAX_BYTES:
                close_part()
                open_new()
            fh.write(row + '\n')
            rows_in_part += 1
    close_part()

    # 回写总卷数到每个文件的首行
    total_parts = len(written)
    for i, p in enumerate(written, 1):
        with open(p, encoding='utf-8') as fh:
            body = fh.read()
        body = body.replace('# circBank human entries (TSV), part %d\n' % i,
                            '# circBank human entries (TSV), part %d of %d\n'
                            % (i, total_parts), 1)
        with open(p, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write(body)

    print('\n生成 %d 卷：' % total_parts)
    tot = 0
    for p in written:
        sz = os.path.getsize(p)
        tot += sz
        print('  %-44s %6.2f MB' % (os.path.basename(p), sz / 1048576))
    print('合计 %.2f MB' % (tot / 1048576))
    return written


if __name__ == '__main__':
    main()