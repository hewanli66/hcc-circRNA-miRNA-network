# -*- coding: utf-8 -*-
"""
解析 GSE122482_RAW.tar：10 对 HBV 相关 HCC 组织 circRNA 表达谱
输出：合并矩阵（circRNA × 20 样本），配对差异分析
"""
import tarfile, gzip, os, math, io, json

RAW = "data/GSE122482_RAW.tar"
OUT = "results"
os.makedirs(OUT, exist_ok=True)

if not os.path.exists(RAW):
    raise SystemExit(f"[ERROR] 未找到 {RAW}")

with tarfile.open(RAW) as tf:
    members = [m for m in tf.getmembers() if m.isfile()]
    print(f"[INFO] tar 内 {len(members)} 个文件")
    for m in members[:30]:
        print(f"    {m.name}  ({m.size} bytes)")

    # ---- 读取每个样本的 RPM 文件 ----
    # 文件名形如 GSM3467136_1232-1_circ_RPM.txt.gz  （-1=肿瘤, -2=癌旁）
    sample_data = {}   # sample_key -> {circRNA: RPM}
    meta = []
    for m in members:
        name = os.path.basename(m.name)
        if not name.endswith(".txt.gz"):
            print(f"[SKIP] 非 RPM 文件: {name}")
            continue
        content = tf.extractfile(m).read()
        text = gzip.decompress(content).decode("utf-8", "ignore")
        lines = text.splitlines()
        # 解析样本元信息
        gsm = None
        pid = None
        region = None
        d = {}
        for l in lines:
            if l.startswith("#") or not l.strip():
                # 元信息行
                if "sample" in l.lower() or "title" in l.lower() or "GSM" in l:
                    pass
                continue
            p = l.split("\t")
            if len(p) < 2:
                continue
            key = p[0].strip()
            try:
                val = float(p[1])
            except ValueError:
                continue
            d[key] = val
        # 从文件名提 GSM 与患者编号
        # 形如 GSM3467136_1232-1_circ_RPM.txt.gz
        base = name[:-len(".txt.gz")]
        if "_" in base:
            gsm, rest = base.split("_", 1)
            # rest 形如 "1232-1_circ_RPM" 或 "1232-2_circ_RPM"
            tok = rest.split("_")[0]          # "1232-1"
            if "-" in tok:
                pid, region = tok.rsplit("-", 1)
            else:
                pid, region = tok, ""
        sample_data[gsm] = d
        meta.append((gsm, pid, region, len(d)))
        print(f"  {gsm}  患者{pid}  区域{region}  circRNA数={len(d)}")

    if not sample_data:
        raise SystemExit("[ERROR] 未解析到任何样本")

    # ---- 合并成矩阵 ----
    all_circ = set()
    for d in sample_data.values():
        all_circ |= set(d.keys())
    all_circ = sorted(all_circ)
    samples = sorted(sample_data.keys())
    print(f"\n[INFO] 合并矩阵：circRNA {len(all_circ)} × 样本 {len(samples)}")

    # 保留在 >=2 个样本中检出的 circRNA
    n_ge2 = 0
    keep = []
    for c in all_circ:
        cnt = sum(1 for s in samples if c in sample_data[s])
        if cnt >= 2:
            keep.append(c)
    print(f"[INFO] 在 >=2 个样本中检出的 circRNA：{len(keep)}")

    # ---- 配对信息 ----
    # 按 GSM 字典：GSM3467136-7145 肿瘤(-1), 7146-7155 癌旁(-2)
    info = {gsm: (pid, region) for gsm, pid, region, _ in meta}
    tumor = [s for s in samples if info[s][1] == "1"]
    normal = [s for s in samples if info[s][1] == "2"]
    print(f"[INFO] 肿瘤样本 {len(tumor)}，癌旁样本 {len(normal)}")

    # 验证配对
    pair_t = {info[s][0]: s for s in tumor}
    pair_n = {info[s][0]: s for s in normal}
    common = sorted(set(pair_t) & set(pair_n))
    print(f"[INFO] 可配对患者：{len(common)} 对 -> {common}")

    # ---- 输出矩阵 ----
    mat = f"{OUT}/GSE122482_circRNA_matrix.tsv"
    with open(mat, "w", encoding="utf-8") as fo:
        fo.write("circRNA\t" + "\t".join(samples) + "\n")
        for c in keep:
            fo.write(c + "\t" + "\t".join(str(sample_data[s].get(c, 0)) for s in samples) + "\n")
    print(f"\n[SAVED] {mat}")

    # ---- 配对差异分析（肿瘤 vs 癌旁）----
    def sf(t, df):
        if df <= 0 or math.isnan(t):
            return float("nan")
        x = df / (df + t * t)
        x = min(max(x, 1e-300), 1 - 1e-15)
        def bc(a, b, x):
            EPS, FP = 3e-16, 1e-300
            qab, qap, qam = a + b, a + 1.0, a - 1.0
            c = 1.0; d = 1.0 - qab * x / qap
            if abs(d) < FP: d = FP
            d = 1.0 / d; h = d
            for m in range(1, 301):
                m2 = 2 * m
                aa = m * (b - m) * x / ((qam + m2) * (a + m2))
                d = 1.0 + aa * d
                if abs(d) < FP: d = FP
                c = 1.0 + aa / c
                if abs(c) < FP: c = FP
                d = 1.0 / d; h *= d * c
                aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
                d = 1.0 + aa * d
                if abs(d) < FP: d = FP
                c = 1.0 + aa / c
                if abs(c) < FP: c = FP
                d = 1.0 / d; de = d * c; h *= de
                if abs(de - 1.0) < EPS: break
            return h
        lb = math.lgamma(df / 2 + .5) - math.lgamma(df / 2) - math.lgamma(.5)
        bt = math.exp(lb + (df / 2) * math.log(x) + .5 * math.log(max(1 - x, 1e-300)))
        a, b = df / 2, .5
        if x < (a + 1) / (a + b + 2):
            return min(1.0, bt * bc(a, b, x) / a)
        return max(0.0, 1 - bt * bc(b, a, 1 - x) / b)

    if common:
        res = []
        for c in keep:
            tv = [sample_data[pair_t[p]].get(c, 0) for p in common]
            nv = [sample_data[pair_n[p]].get(c, 0) for p in common]
            diffs = [tv[k] - nv[k] for k in range(len(common))]
            md = sum(diffs) / len(diffs)
            var = sum((d - md) ** 2 for d in diffs) / (len(diffs) - 1)
            se = math.sqrt(var / len(diffs))
            t = md / se if se > 0 else 0.0
            lfc = math.log2((sum(tv) / len(tv) + 1) / (sum(nv) / len(nv) + 1))
            res.append((c, sum(nv) / len(nv), sum(tv) / len(tv), lfc, t, sf(t, len(diffs) - 1)))

        # BH 校正
        order = sorted(range(len(res)), key=lambda i: res[i][5])
        m_ = len(res); adj = [1.0] * m_; prev = 1.0
        for rank in range(m_, 0, -1):
            v = min(prev, res[order[rank - 1]][5] * m_ / rank)
            adj[order[rank - 1]] = v; prev = v
        padj = {order[i]: adj[i] for i in range(m_)}

        up = dn = 0
        with open(f"{OUT}/GSE122482_circRNA_DEG.txt", "w", encoding="utf-8") as fo:
            fo.write("circRNA\tmean_normal_RPM\tmean_tumor_RPM\tlog2FC\tt\tP\tP_adj\tregulation\n")
            for i, r in enumerate(sorted(res, key=lambda x: -x[3])):
                pa = padj[i]
                reg = "-"
                if abs(r[3]) >= 1.0 and pa < 0.05:
                    reg = "up" if r[3] > 0 else "down"
                    if reg == "up": up += 1
                    else: dn += 1
                fo.write(f"{r[0]}\t{r[1]:.3f}\t{r[2]:.3f}\t{r[3]:.4f}\t{r[4]:.3f}\t{r[5]:.6g}\t{pa:.6g}\t{reg}\n")
        print(f"\n[RESULT] |log2FC|>=1 且 P_adj<0.05 → 上调 {up}，下调 {dn}（总数 {len(res)}）")
        print("[SAVED] results/GSE122482_circRNA_DEG.txt")

        print("\n=== 差异最显著的 circRNA（前 20，按 |log2FC| 排序）===")
        print(f"{'circRNA':<28}{'normal':>11}{'tumor':>11}{'log2FC':>9}{'P':>11}")
        for r in sorted(res, key=lambda x: -abs(x[3]))[:20]:
            print(f"{r[0]:<28}{r[1]:11.2f}{r[2]:11.2f}{r[3]:9.3f}{r[5]:11.3g}")

    json.dump({s: info[s] for s in samples}, open(f"{OUT}/GSE122482_sample_info.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print(f"\n[SAVED] {OUT}/GSE122482_sample_info.json")