"""Compose Fig. 7 (manuscript v0.7) from the 0.6.0 browser screenshots.

Screenshots were taken at deviceScaleFactor 2 by fig7_capture.mjs. Interface clips are
pasted at 1:1 so that interface text prints at about 6 pt at 170 mm width (300 dpi); the
neighbor-evidence legend is redrawn at the same size (review items S13, S14, S41)."""
import os
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw, ImageFont

FIG_DIR = Path(__file__).resolve().parent.parent
S = Path(os.environ.get('FIG7_SCREENS', FIG_DIR / 'fig7_work' / 'fig7_screens'))
OUT = Path(os.environ.get('FIG7_OUTPUT', FIG_DIR / 'output_v07'))
FONT = os.environ.get('FIG7_FONT', '/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf')  # Arial-metric; on macOS set Arial.ttf
BOLD = os.environ.get('FIG7_FONT_BOLD', '/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf')
WIDTH = 2008  # 170 mm at 300 dpi
f_text, f_bold, f_label = ImageFont.truetype(FONT, 26), ImageFont.truetype(BOLD, 26), ImageFont.truetype(BOLD, 56)


def content_box(img, threshold=225, min_count=12, pad=10):
    """Bounding box of the map cells; the faint background dot grid is ignored."""
    a = np.asarray(img.convert('L')) < threshold
    rows, cols = np.where(a.sum(axis=1) >= min_count)[0], np.where(a.sum(axis=0) >= min_count)[0]
    return (max(cols.min() - pad, 0), max(rows.min() - pad, 0), min(cols.max() + pad, img.width), min(rows.max() + pad, img.height))


# a: map of the neighbor-evidence layer with the selected compound (frame line at the top removed).
m = Image.open(S / 'fig7_map_selected.png').convert('RGB')
m = m.crop((6, 6, m.width - 6, m.height - 6))
m = m.crop(content_box(m))
map_w = 1010
m = m.resize((map_w, round(m.height * map_w / m.width)), Image.LANCZOS)

# b: key inspector fields at 1:1; rows from "Structure group" to "Evidence units" omitted.
q = Image.open(S / 'fig7_query_summary.png').convert('RGB')
segments = [q.crop((0, 6, q.width, 100)), q.crop((0, 392, q.width, 592)), q.crop((0, 1000, q.width, q.height))]
b = Image.new('RGB', (q.width, sum(seg.height for seg in segments) + 24 * (len(segments) - 1)), 'white')
d = ImageDraw.Draw(b)
yb = 0
for i, seg in enumerate(segments):
    b.paste(seg, (0, yb)); yb += seg.height
    if i < len(segments) - 1:
        for x in range(8, q.width - 8, 16):
            d.line([(x, yb + 12), (x + 8, yb + 12)], fill='#9aa5b1', width=2)
        yb += 24

# legend column, redrawn from the 0.6.0 interface legend
legend_w = WIDTH - 2 * 30 - map_w - b.width - 2 * 36
leg = Image.new('RGB', (legend_w, max(m.height, b.height)), 'white')
d = ImageDraw.Draw(leg)


def wrap(text, font, width):
    words, lines, line = text.split(), [], ''
    for w in words:
        trial = (line + ' ' + w).strip()
        if d.textlength(trial, font=font) <= width:
            line = trial
        else:
            lines.append(line); line = w
    lines.append(line)
    return lines


def swatch(x, y, frame, style='solid', width=6, ring=True, fill='white'):
    s = 44
    if ring:
        d.rectangle([x, y, x + s, y + s], fill='#1f2933')
        d.rectangle([x + 3, y + 3, x + s - 3, y + s - 3], fill=fill)
        inner = 3
    else:
        d.rectangle([x, y, x + s, y + s], fill=fill); inner = 0
    if style == 'solid':
        for k in range(width):
            d.rectangle([x + inner + k, y + inner + k, x + s - inner - k, y + s - inner - k], outline=frame)
    else:
        dash, gap = (12, 8) if style == 'dashed' else (3, 8)
        for k in range(width):
            o = inner + k
            for t in range(0, s - 2 * o, dash + gap):
                e = min(t + dash, s - 2 * o)
                d.line([(x + o + t, y + o), (x + o + e, y + o)], fill=frame)
                d.line([(x + o + t, y + s - o), (x + o + e, y + s - o)], fill=frame)
                d.line([(x + o, y + o + t), (x + o, y + o + e)], fill=frame)
                d.line([(x + s - o, y + o + t), (x + s - o, y + o + e)], fill=frame)


y = 4
d.text((0, y), 'Neighbor evidence (frame)', font=f_bold, fill='#1f2933'); y += 44
items = [('#2a9d8f', 'solid', True, 'Consistent: active'), ('#e6a23c', 'solid', True, 'Consistent: medium'),
         ('#c44e52', 'solid', True, 'Consistent: inactive'), ('#7a5195', 'dashed', False, 'Mixed evidence'),
         ('#8a99a6', 'dotted', False, 'Insufficient evidence'), ('#b8c2cc', 'solid', False, 'No similar reference molecule')]
for color, style, ring, text in items:
    swatch(0, y, color, style, 6 if text != 'No similar reference molecule' else 2, ring)
    lines = wrap(text, f_text, legend_w - 60)
    for i, line in enumerate(lines):
        d.text((58, y + 8 + 30 * i), line, font=f_text, fill='#1f2933')
    y += max(60, 30 * len(lines) + 24)
y += 10
d.text((0, y), 'Fill', font=f_bold, fill='#1f2933'); y += 40
for line in wrap("The new compound's own measured class (pChEMBL: inactive ≤ 6, medium 6–7, active ≥ 7); white if not measured.", f_text, legend_w):
    d.text((0, y), line, font=f_text, fill='#1f2933'); y += 30
y += 16
d.rectangle([0, y, 44, y + 44], fill='#cfcfcf')
for i, line in enumerate(wrap('Reference compounds: grey (classes in the other layers)', f_text, legend_w - 60)):
    d.text((58, y + 2 + 30 * i), line, font=f_text, fill='#1f2933')
y += 100
d.rectangle([0, y, 44, y + 44], outline='#7a5195', width=3)
for i, line in enumerate(wrap('Five most similar reference compounds (dotted)', f_text, legend_w - 60)):
    d.text((58, y + 2 + 30 * i), line, font=f_text, fill='#1f2933')
y += 100
d.rectangle([0, y, 44, y + 44], outline='#183e92', width=4)
for i, line in enumerate(wrap('Selected compound', f_text, legend_w - 60)):
    d.text((58, y + 8 + 30 * i), line, font=f_text, fill='#1f2933')
y += 70
leg = leg.crop((0, 0, legend_w, max(y, 10)))

# c: new-molecule list (header, counts and the first rows)
lst = Image.open(S / 'fig7_new_list.png').convert('RGB')
c = lst.crop((40, 150, 1880, 1195))

label_h, gap, margin = 76, 36, 30
row1_h = max(m.height, b.height, leg.height)
height = margin + label_h + row1_h + gap + label_h + c.height + margin
fig = Image.new('RGB', (WIDTH, height), 'white')
draw = ImageDraw.Draw(fig)
y0 = margin + label_h
x = margin
fig.paste(m, (x, y0)); draw.text((x, margin), 'a', fill='black', font=f_label)
x += m.width + gap
fig.paste(leg, (x, y0 + 10))
x += leg.width + gap
fig.paste(b, (x, y0)); draw.text((x, margin), 'b', fill='black', font=f_label)
y1 = y0 + row1_h + gap
draw.text((margin, y1), 'c', fill='black', font=f_label)
fig.paste(c, (margin, y1 + label_h))
fig.save(OUT / 'Figure_7.png', dpi=(300, 300))
fig.save(OUT / 'Figure_7.pdf', resolution=300)
print(fig.size, 'mm', round(fig.width / 300 * 25.4), round(fig.height / 300 * 25.4))
