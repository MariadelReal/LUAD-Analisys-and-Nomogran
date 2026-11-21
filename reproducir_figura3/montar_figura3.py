"""Exact raster composition of supplied Cytoscape exports; no network redrawing.
Run: python montar_figura3.py (requires Pillow and reportlab).
"""
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output"
OUT.mkdir(exist_ok=True)
DPI = 600
W, H = round(240 / 25.4 * DPI), round(260 / 25.4 * DPI)
FONT = '/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf'
BOLD = '/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf'
def font(size, bold=False):
    for candidate in ([BOLD if bold else FONT, 'C:/Windows/Fonts/arialbd.ttf' if bold else 'C:/Windows/Fonts/arial.ttf', 'DejaVuSans.ttf']):
        try:
            return ImageFont.truetype(candidate, size)
        except OSError:
            pass
    raise RuntimeError('Install Arial, Liberation Sans or DejaVu Sans.')

im = Image.new('RGB', (W,H), 'white')
d = ImageDraw.Draw(im)
margin, gap = 140, 180
cw = (W-2*margin-gap)//2
for x, title in [(margin,'<60 years'), (margin+cw+gap,'≥60 years')]:
    d.text((x+cw/2,100), title, font=font(92), fill='black', anchor='mt')
panels = [('a','young_module.png',0,300), ('b','old_module.png',1,300),
          ('c','young_hubs.png',0,3190), ('d','old_hubs.png',1,3190)]
for label, filename, col, y in panels:
    x = margin + col*(cw+gap)
    d.text((x,y), label, font=font(110,True), fill='black')
    raw = Image.open(ROOT/'input'/filename).convert('RGBA')
    # Crop only fully transparent margins; retain all visible source pixels.
    raw = raw.crop(raw.getchannel('A').getbbox())
    bg = Image.new('RGBA',raw.size,'white')
    bg.alpha_composite(raw)
    raw = bg.convert('RGB')
    maxw,maxh = cw-70,2630
    scale = min(maxw/raw.width,maxh/raw.height)
    raw = raw.resize((round(raw.width*scale),round(raw.height*scale)),Image.Resampling.LANCZOS)
    im.paste(raw,(x+(cw-raw.width)//2,y+150+(maxh-raw.height)//2))
im.save(OUT/'Fig3_KIF4A_600ppi.tiff',compression='tiff_lzw',dpi=(DPI,DPI))
im.save(OUT/'Fig3_KIF4A_600ppi.png',dpi=(DPI,DPI))
pdf=canvas.Canvas(str(OUT/'Fig3_KIF4A_600ppi.pdf'),pagesize=(W/DPI*72,H/DPI*72))
pdf.drawImage(ImageReader(im),0,0,width=W/DPI*72,height=H/DPI*72)
pdf.showPage();pdf.save()
preview=im.copy();preview.thumbnail((1600,1800));preview.save(OUT/'Fig3_vista_previa.png')
print(f'Exported {W} × {H} pixels; {DPI} ppi; 240 × 260 mm nominal.')
