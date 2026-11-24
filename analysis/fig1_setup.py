"""Figure 1: (a) 3D schematic (user's SVG, rendered) and (b) cropped, labelled top-view photograph."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))) + '/'   # repository root
import sys
import numpy as np
from PIL import Image, ImageOps
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

SCH = sys.argv[1] if len(sys.argv) > 1 else ROOT + 'figs/schematic_3d.png'   # rendered schematic PNG
PHOTO = ROOT + 'photos/setup/20260930_172634.jpg'
OUT = ROOT + 'figs/Fig1_setup'

sch = Image.open(SCH).convert('RGB'); a = np.asarray(sch)
ys, xs = np.nonzero((a < 245).any(2)); m = 30
sch = sch.crop((xs.min() - m, ys.min() - m, xs.max() + m, ys.max() + m))

# photo: lamp-to-prism part of the breadboard only (grating mount and loose cables cropped out), beam left -> right
ph = ImageOps.exif_transpose(Image.open(PHOTO)).crop((750, 560, 2000, 2900)).rotate(90, expand=True)
k = 2340 / 1400; X0, Y0 = int(30 * k), int(95 * k)
ph = ph.crop((X0, Y0, int(1400 * k), int(705 * k)))
P = lambda x, y: (x * k - X0, y * k - Y0)            # preview coords -> cropped coords
labels = [('Lamp', 80, 395), ('Condenser\nlens', 310, 390), ('Iris', 430, 395), ('Focusing\nlens', 630, 390),
          ('Sample', 735, 380), ('Lens\nf = 1000 mm', 905, 375), ('Prism', 1140, 330)]

plt.rcParams.update({'font.family': 'DejaVu Sans'})
fig = plt.figure(figsize=(7.2, 6.4))
r1 = sch.height / sch.width; r2 = ph.height / ph.width
gs = fig.add_gridspec(2, 1, height_ratios=[r1, r2], hspace=0.12)
ax = fig.add_subplot(gs[0]); ax.imshow(sch); ax.set_axis_off()
ax.text(0.0, 1.0, 'a', transform=ax.transAxes, fontsize=13, fontweight='bold', va='top')
ax = fig.add_subplot(gs[1]); ax.imshow(ph); ax.set_axis_off()
ax.text(0.0, 1.0, 'b', transform=ax.transAxes, fontsize=13, fontweight='bold', va='bottom')
for i, (t, x, y) in enumerate(labels):
    px, py = P(x, y); ty = P(0, 150 if i % 2 == 0 else 215)[1]
    ax.annotate(t.replace('\n', ' '), xy=(px, py), xytext=(px, ty), ha='center', va='bottom', fontsize=6, color='white',
                bbox=dict(boxstyle='round,pad=0.25', fc='black', ec='none', alpha=0.65),
                arrowprops=dict(arrowstyle='-', color='white', lw=0.6, shrinkA=0, shrinkB=2))
fig.savefig(OUT + '.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(OUT + '.pdf', bbox_inches='tight', facecolor='white')
print('ok')
