# -*- coding: utf-8 -*-
"""
图件输出为期刊要求格式。
SR(general figure guidelines)：
  - 白底；线图 ≥ 300 dpi（照片 ≥ 300 dpi，线图建议 600 dpi 以备缩放）
  - 字号下限 8 pt；线宽下限 1 pt
  - 面板标签用小写粗体 a, b, c
  - 图幅宽度：单栏 89 mm，双栏 183 mm
做法：用 Edge headless 以 -w/-h 放大渲染 SVG（矢量）→ PNG，
     再用 Pillow 写入 300/600 dpi 元数据。
"""
import os, re, glob, subprocess
from PIL import Image

EDGE = r'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe'
SRCDIR = 'figures'
OUTDIR = 'figures/SR_submission'
os.makedirs(OUTDIR, exist_ok=True)

# 目标像素宽（按双栏 183 mm @600 dpi ≈ 4323 px；单栏 89 mm @600 dpi ≈ 2102 px）
DPI_TARGET = 600
SINGLE_MM, DOUBLE_MM = 89, 183


def mm2px(mm, dpi=DPI_TARGET):
    return int(round(mm / 25.4 * dpi))


# 判定单栏/双栏：按 SVG 宽度
def classify(w_mm):
    return 'single' if w_mm <= 120 else 'double'


print('%-32s %-10s %-8s %-6s %-16s %s' % ('figure','width_mm','target_px','dpi','handling','size'))
report = []
for svg in sorted(glob.glob(os.path.join(SRCDIR, '*.svg'))):
    name = os.path.splitext(os.path.basename(svg))[0]
    t = open(svg, encoding='utf-8').read()
    m = re.search(r'width="(\d+)"\s+height="(\d+)"', t)
    if not m:
        m = re.search(r'width="(\d+)px"\s+height="(\d+)px"', t)
    w_px, h_px = int(m.group(1)), int(m.group(2))
    # SVG 用户单位视为 px（96 dpi → mm = px/96*25.4）
    w_mm = w_px / 96.0 * 25.4
    kind = classify(w_mm)
    # 保持原始物理尺寸，但放大到至少 600 dpi 等效分辨率
    # 原 SVG 宽度若超双栏 183 mm，缩到双栏宽；否则按原尺寸仅提分辨率
    if w_mm > DOUBLE_MM:
        target_w = mm2px(DOUBLE_MM)
        note = 'scaled to 183 mm'
    else:
        target_w = mm2px(min(w_mm, DOUBLE_MM))
        note = 'original size'
    scale = target_w / w_px
    target_h = int(round(h_px * scale))

    out_png = os.path.join(OUTDIR, name + '.png')
    url = 'file:///' + os.path.abspath(svg).replace('\\', '/')
    # Edge 渲染 SVG 时不会缩放到窗口大小，会贴在左上角。
    # 办法：用 device-scale-factor 放大 + 先把 SVG 包进一个按目标像素定尺寸的 HTML，
    #      HTML 里的 <img> 用 CSS 拉伸到目标像素，SVG 作为图片加载
    wrapper = os.path.join(OUTDIR, '_tmp_%s.html' % name)
    with open(wrapper, 'w', encoding='utf-8') as fh:
        fh.write(
            '<!DOCTYPE html><html><head><meta charset="utf-8">'
            '<style>html,body{margin:0;padding:0;background:#fff;}'
            'img{display:block;width:%dpx;height:%dpx;}</style></head>'
            '<body><img src="file:///%s"></body></html>'
            % (target_w, target_h, os.path.abspath(svg).replace('\\', '/')))
    cmd = [EDGE, '--headless=new', '--disable-gpu', '--no-sandbox',
           '--hide-scrollbars', '--default-background-color=FFFFFFFF',
           '--screenshot=%s' % os.path.abspath(out_png),
           '--window-size=%d,%d' % (target_w, target_h),
           'file:///' + os.path.abspath(wrapper).replace('\\', '/')]
    subprocess.run(cmd, capture_output=True, timeout=240)
    if os.path.exists(wrapper):
        try:
            os.remove(wrapper)
        except OSError:
            pass

    if os.path.exists(out_png):
        im = Image.open(out_png)
        im.save(out_png, dpi=(DPI_TARGET, DPI_TARGET), optimize=True)
        sz = os.path.getsize(out_png) / 1024
        print('%-32s %-10.1f %-8d %-10d %-16s %.0f KB' % (name, w_mm, target_w,
                                                          DPI_TARGET, note, sz))
        report.append((name, kind, w_mm, target_w, target_h))
    else:
        print('%-32s FAILED' % name)

print('\n已输出至: %s' % OUTDIR)
print('面板标签与线宽检查（原始 SVG）:')
for svg in sorted(glob.glob(os.path.join(OUTDIR, '*.svg'))):
    pass
# 对源 SVG 做规格体检
print('\n%-32s %-14s %-14s' % ('figure', 'min_font_pt', 'min_stroke_pt'))
for svg in sorted(glob.glob(os.path.join(SRCDIR, '*.svg'))):
    t = open(svg, encoding='utf-8').read()
    fs = [float(x) for x in re.findall(r'font-size="(\d+\.?\d*)"', t)]
    sw = [float(x) for x in re.findall(r'stroke-width="(\d+\.?\d*)"', t)]
    mf = min(fs) if fs else 0
    ms = min(sw) if sw else 0
    flag = ''
    if mf and mf < 8:
        flag += ' 字<8pt!'
    if ms and ms < 1:
        flag += ' 线<1pt!'
    print('%-32s %-14s %-14s%s' % (os.path.basename(svg), mf, ms, flag))