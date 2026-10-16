import os
from pathlib import Path
import pymupdf as fitz
import numpy as np
ROOT = Path(os.environ.get('DJK_PROJECT_ROOT', '.'))
OUT = ROOT / 'Figure5_等高紧凑排版_20260921'
OUT.mkdir(exist_ok=True)
BASE = ROOT / '正文主图单图_大字号_最终版_20260903/Figure_5'
left = fitz.open(ROOT / '蓝藻相对丰度变化曲线_细高版_20260908/Cyanobacteria_relative_abundance_vertical_curve.pdf')
names = [('B', 'Fig5D_WL_leads_Chla'), ('C', 'Fig5B_Cyanobacteria_leads_Chla'), ('D', 'Fig5A_TNTP_mass_ratio_leads_Chla'), ('E', 'Fig5E_WL_leads_Cyanobacteria'), ('F', 'Fig5C_Air_temperature_leads_Cyanobacteria'), ('G', 'Fig5F_Water_temperature_leads_Cyanobacteria')]
docs = [fitz.open(BASE / (n + '.pdf')) for _, n in names]

def ink_box(page, clip=None):
    pix = page.get_pixmap(matrix=fitz.Matrix(2, 2), clip=clip, alpha=False)
    a = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, 3)
    yy, xx = np.where(np.any(a < 230, axis=2))
    return fitz.Rect((pix.x + xx.min()) / 2 - 2, (pix.y + yy.min()) / 2 - 2, (pix.x + xx.max() + 1) / 2 + 2, (pix.y + yy.max() + 1) / 2 + 2)
lc = ink_box(left[0], fitz.Rect(0, 57, 259, 676))
boxes = [ink_box(d[0]) for d in docs]
rc = fitz.Rect(min((b.x0 for b in boxes)), min((b.y0 for b in boxes)), max((b.x1 for b in boxes)), max((b.y1 for b in boxes)))
scale = 0.63
pad = 9
left_scale = 1.0
top = pad + (59.93 - lc.y0) * left_scale
bottom = top + (606.512 - 59.93) * left_scale
axis_h = (254.879 - 17.281) * scale
step = (bottom - top - axis_h) / 2
left_width = lc.width * left_scale
panel_width = rc.width * scale
col_gap = 12
block_width = panel_width * 2 + col_gap
height = lc.height * left_scale + 2 * pad
start_right = pad + left_width + 16
width = start_right + block_width + pad

def place(page, doc, clip, x, y, s):
    page.show_pdf_page(fitz.Rect(x, y, x + clip.width * s, y + clip.height * s), doc, 0, clip=clip)

def draw_right(page, start):
    for i, ((letter, _), doc) in enumerate(zip(names, docs)):
        col, row = divmod(i, 3)
        x = start + col * (panel_width + col_gap)
        y = top + row * step - (17.281 - rc.y0) * scale
        place(page, doc, rc, x, y, scale)
        page.insert_text((x + 1, top + row * step + 5), letter, fontsize=12, fontname='helv')

def export(doc, stem):
    doc.save(OUT / (stem + '.pdf'), garbage=4, deflate=True)
    doc[0].get_pixmap(dpi=600, alpha=False).save(OUT / (stem + '.png'))
    (OUT / (stem + '.svg')).write_text(doc[0].get_svg_image(), encoding='utf-8')
    doc[0].get_pixmap(dpi=130, alpha=False).save(OUT / (stem + '_preview.png'))
full = fitz.open()
p = full.new_page(width=width, height=height)
place(p, left, lc, pad, pad, left_scale)
draw_right(p, start_right)
export(full, 'Figure5_aligned')
right = fitz.open()
p = right.new_page(width=block_width + pad * 2, height=height)
draw_right(p, pad)
export(right, 'Figure5_right_six_aligned')
print('Plot frame top/bottom aligned:', top, bottom)
print('Right panel scale:', scale, 'row step:', step, 'canvas:', width, height)
print(OUT)
