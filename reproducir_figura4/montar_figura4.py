"""Reassemble Figure 4 from the supplied model-output PNGs; no refitting."""
from pathlib import Path
import io
import numpy as np
from PIL import Image, ImageDraw, ImageFont, ImageChops
import fitz
ROOT=Path(__file__).resolve().parent
FIG=ROOT/'output';FIG.mkdir(exist_ok=True)
QA=FIG
FONT=None;BOLD=None
for regular,bold in [('/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf','/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf'),('C:/Windows/Fonts/arial.ttf','C:/Windows/Fonts/arialbd.ttf')]:
    if Path(regular).exists() and Path(bold).exists(): FONT,BOLD=regular,bold;break
if FONT is None: raise RuntimeError('Install Arial or Liberation Sans.')
def trim(im,pad=20):
    box=ImageChops.difference(im,Image.new('RGB',im.size,'white')).getbbox()
    x0,y0,x1,y1=box
    return im.crop((max(0,x0-pad),max(0,y0-pad),min(im.width,x1+pad),min(im.height,y1+pad)))
# Figure 4c: replace only the crowded text labels with staggered labels.
nom=Image.open(ROOT/'input/Fig4c_Nomogram(1).png').convert('RGB')
assert nom.size==(2048,1303)
before=np.asarray(nom).copy()
dr=ImageDraw.Draw(nom);dr.rectangle((884,900,1735,958),fill='white')
ticks=[(918,'0.95',902),(1137,'0.9',902),(1269,'0.85',902),
       (1364,'0.8',902),(1441,'0.75',902),(1506,'0.7',929),
       (1564,'0.65',902),(1615,'0.6',929),(1663,'0.55',902),(1708,'0.5',929)]
for x,lab,y in ticks:dr.text((x,y),lab,font=ImageFont.truetype(FONT,26),fill='black',anchor='mt')
diff=np.any(np.asarray(nom)!=before,axis=2)
diff[900:959,884:1736]=False
assert not diff.any(), 'Only nomogram tick-label rectangle may change.'

# Figure 4d: keep all plotted pixels and tick values intact; place footer below xlabel.
cal=Image.open(ROOT/'input/Fig4d_Calibration_5y.png').convert('RGB')
assert cal.size==(1300,1100)
cal_new=Image.new('RGB',(1300,1180),'white');cal_new.paste(cal.crop((0,0,1300,1000)),(0,0))
dr=ImageDraw.Draw(cal_new)
dr.text((690,1024),'Predicted 5-year OS',font=ImageFont.truetype(FONT,32),fill='black',anchor='mt')
f=ImageFont.truetype(FONT,22)
dr.text((162,1082),'n=474 d=175 p=9, 80 subjects per group',font=f,fill='black')
dr.text((162,1114),'Gray: ideal',font=f,fill='black')
dr.text((1216,1082),'X - resampling optimism added, B=200',font=f,fill='black',anchor='rt')
dr.text((1216,1114),'Based on observed-predicted',font=f,fill='black',anchor='rt')
assert np.array_equal(np.asarray(cal_new)[:1000],np.asarray(cal)[:1000])

panels=[trim(Image.open(ROOT/'input/Cox_univariate_30genes_forest(20261117-224302).png').convert('RGB')),
        trim(Image.open(ROOT/'input/Fig4b_LASSO_CV(1).png').convert('RGB')),
        trim(nom),trim(cal_new),trim(Image.open(ROOT/'input/Fig4e_Time_ROC(1).png').convert('RGB'))]
W=3600;M=55;G=70;column=(W-2*M-G)//2;label_h=65
def at_width(im,w):return im.resize((w,round(im.height*w/im.width)),Image.Resampling.LANCZOS)
scaled=[at_width(panels[0],column),at_width(panels[1],column),
        at_width(panels[2],W-2*M),at_width(panels[3],column),at_width(panels[4],column)]
row1=max(scaled[0].height,scaled[1].height)+label_h
row2=scaled[2].height+label_h
row3=max(scaled[3].height,scaled[4].height)+label_h
H=M+row1+G+row2+G+row3+M
mont=Image.new('RGB',(W,H),'white');dr=ImageDraw.Draw(mont)
ys=[M,M,M+row1+G,M+row1+G+row2+G,M+row1+G+row2+G]
xs=[M,M+column+G,M,M,M+column+G]
for lab,p,x,y in zip('abcde',scaled,xs,ys):
    dr.text((x,y),lab,font=ImageFont.truetype(BOLD,52),fill='black')
    mont.paste(p,(x,y+label_h))
mont.save(FIG/'Fig4_RSOS.tiff',compression='tiff_lzw',dpi=(300,300))
mont.save(FIG/'Fig4_RSOS.png',dpi=(300,300))
buf=io.BytesIO();mont.save(buf,format='PNG')
pdf=fitz.open();p=pdf.new_page(width=W*72/300,height=H*72/300);p.insert_image(p.rect,stream=buf.getvalue());pdf.save(FIG/'Fig4_RSOS.pdf',deflate=True)
preview=mont.copy();preview.thumbnail((1400,1700));preview.save(QA/'Fig4_RSOS_preview.png')
nom.save(QA/'Nomogram_labels_QA.png');cal_new.save(QA/'Calibration_footer_QA.png')
