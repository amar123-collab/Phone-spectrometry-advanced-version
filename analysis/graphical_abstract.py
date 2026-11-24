"""Graphical abstract: 3D schematic -> phone photographs -> linearity with and without gamma correction."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))) + '/'   # repository root
import json
import numpy as np
from PIL import Image, ImageOps
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

R = json.load(open(ROOT + 'analysis/results.json'))
PH = ROOT + 'photos/raw/ern_photos/'
OUT = ROOT + 'figs/Graphical_abstract'
C_ENC, C_LIN, C_UV = '#8a8a8a', '#c2410c', 'black'
CROP = (600, 3230, 2050, 3700)

sch = Image.open(ROOT + 'figs/schematic_3d.png').convert('RGB'); a = np.asarray(sch)
ys, xs = np.nonzero((a < 245).any(2)); sch = sch.crop((xs.min() - 20, ys.min() - 20, xs.max() + 20, ys.max() + 20))
photo = lambda n: np.asarray(ImageOps.exif_transpose(Image.open(PH + n + '.jpg')).convert('RGB').crop(CROP))

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8})
fig = plt.figure(figsize=(9.0, 4.2))
fig.text(0.015, 0.955, 'Gamma-aware smartphone absorption spectroscopy with a spare-parts prism spectrometer',
         fontsize=11.5, fontweight='bold', va='top')

# 1: instrument
ax = fig.add_axes([0.0, 0.10, 0.56, 0.74]); ax.imshow(sch); ax.set_axis_off()
fig.text(0.015, 0.84, '1  Project the spectrum, photograph it with a fixed exposure', fontsize=8.5, color='0.25')

# 2: photographs
fig.text(0.605, 0.84, '2  Reference and sample photographs', fontsize=8.5, color='0.25')
for k, (n, lab) in enumerate([('20261120_003701', 'Ethanol (I$_0$)'), ('20261120_003349', 'DCM 49.5 µM (I)')]):
    ax = fig.add_axes([0.605, 0.705 - k * 0.115, 0.25, 0.10]); ax.imshow(photo(n)); ax.set_axis_off()
    fig.text(0.865, 0.755 - k * 0.115, lab, fontsize=7.5, va='center')

# 3: linearity
fig.text(0.605, 0.47, '3  Linearize pixel values before taking −log(I/I$_0$)', fontsize=8.5, color='0.25')
ax = fig.add_axes([0.645, 0.11, 0.2, 0.30])
uv = np.array(R['dcm_uvA']); ax.plot([0, 2], [0, 2], color='0.75', lw=0.8, ls=':')
for tag, c, mk, lab in [('lin', C_LIN, 'o', 'sRGB-linearized'), ('enc', C_ENC, 'D', 'encoded\n(no correction)')]:
    d = R['dcm_' + tag]; p = d['fit_all']
    ax.plot(uv, d['mean'], mk, color=c, mec='white', mew=0.5, ms=5)
    ax.plot([0, 2], np.polyval([p['slope'], p['intercept']], [0, 2]), color=c, lw=1)
    fig.text(0.865, 0.36 if tag == 'lin' else 0.21, f"{lab}\nslope {p['slope']:.2f}", color=c, fontsize=7.5, va='center',
             fontweight='bold' if tag == 'lin' else 'normal')
ax.set_xlim(0, 2); ax.set_ylim(0, 2.1); ax.set_xticks([0, 1, 2]); ax.set_yticks([0, 1, 2])
ax.set_xlabel('Absorbance, UV–vis', fontsize=7.5); ax.set_ylabel('Absorbance, phone', fontsize=7.5)
ax.tick_params(labelsize=7); [ax.spines[s].set_visible(False) for s in ('top', 'right')]
fig.text(0.865, 0.10, '1:1 line dotted', color='0.5', fontsize=7)

for (x0, y0), (x1, y1) in [((0.555, 0.55), (0.598, 0.65))]:
    fig.patches.append(FancyArrowPatch((x0, y0), (x1, y1), transform=fig.transFigure, arrowstyle='-|>',
                                       mutation_scale=12, color='0.45', lw=1.2))
fig.savefig(OUT + '.png', dpi=300, facecolor='white')
fig.savefig(OUT + '.pdf', facecolor='white')
print('ok')
