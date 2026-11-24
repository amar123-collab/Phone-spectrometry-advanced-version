"""SI figure: the 2025 mixed-dye set as used at Researchers' Night with the basic version (v1.0.0) of the code.
Columns: cuvette | projected spectrum on the screen (automatic exposure) | absorbance from the basic-version script run_spectra.m
(v1_as_used, run unchanged except for plotting; pixel axis, its own window baseline) | laboratory UV-vis spectrum (normalized).
Photos and UV-vis files: example_data of the basic record (10.5281/zenodo.23400260), to be copied to basic_example_data/;
cuvette photo photos/setup/20260930_195347.jpg. basic2025/v1_A.csv: absorbance saved by run_spectra.m (v1_as_used) run on those photos."""
import numpy as np, cv2, csv
from PIL import Image, ImageOps
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))) + '/'
BASIC = ROOT + 'basic_example_data/'          # example_data folder of the basic record (10.5281/zenodo.23400260), copied here
CUV = ROOT + 'photos/setup/20260930_195347.jpg'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.5, 'axes.linewidth': 0.7, 'axes.spines.top': False, 'axes.spines.right': False})
COL = [(0.17, 0.31, 0.84), (0.83, 0.66, 0), (0.04, 0.16, 0.42), (0.48, 0.18, 0.82), (0.91, 0.25, 0.54)]   # colors of the basic version
NAMES = ['blue', 'yellow', 'dark blue', 'violet', 'pink']
FILES = ['01_reference_white', '02_sample1_a', '04_sample2', '05_sample3', '06_sample4', '07_sample5']
SHIFT = [(0, 0), (-78, -261), (-22, -51), (-31, 199), (0, -72), (20, 44)]     # (dy, dx) reported by run_spectra.m
BOX = (820, 2010, 1860, 2400)                                                  # crop on the reference photo (x0, y0, x1, y1)
POLY = [[(328, 1000), (792, 1000), (790, 2190), (770, 2265), (365, 2265), (345, 2190)],
        [(1105, 985), (1512, 985), (1508, 2190), (1490, 2262), (1130, 2262), (1112, 2190)],
        [(1945, 990), (2305, 990), (2300, 2185), (2285, 2252), (1965, 2252), (1950, 2185)],
        [(2705, 995), (3110, 995), (3105, 2215), (3085, 2292), (2740, 2292), (2718, 2215)],
        [(3400, 995), (3822, 995), (3818, 2225), (3795, 2302), (3430, 2302), (3410, 2225)]]

def cuvettes():
    img = np.asarray(ImageOps.exif_transpose(Image.open(CUV)).convert('RGB')).astype(float); out = []
    for poly in POLY:
        p = np.array(poly); x0, y0 = p.min(0) - 60; x1, y1 = p.max(0) + 60
        crop = img[y0:y1, x0:x1]; h, w = crop.shape[:2]
        m = np.zeros((h, w), np.uint8); cv2.fillPoly(m, [(p - [x0, y0]).astype(np.int32)], 255)
        bg = crop[cv2.dilate(m, np.ones((25, 25))) == 0]; c = np.clip(crop * (238.0 / np.median(bg, 0)), 0, 255)
        a = cv2.GaussianBlur(m.astype(float) / 255, (0, 0), 3)[..., None]; comp = c * a + 255 * (1 - a)
        ys, xs = np.nonzero(m); out.append(comp[ys.min():ys.max() + 1, xs.min() - 2:xs.max() + 3].astype(np.uint8))
    return out

def screen(k):
    im = ImageOps.exif_transpose(Image.open(BASIC + 'photos/' + FILES[k] + '.jpg')).convert('RGB')
    dy, dx = SHIFT[k]; x0, y0, x1, y1 = BOX
    return np.asarray(im.crop((x0 + dx, y0 + dy, x1 + dx, y1 + dy)))     # follow the spectrum as the script did

def lab(k):
    rows = list(csv.reader(open(BASIC + f'uvvis/sample{k}.csv', encoding='utf-8')))[1:]
    w, a = np.array([[float(r[0]), float(r[1])] for r in rows]).T
    a = a - a[(w >= 750) & (w <= 800)].mean(); m = (w >= 400) & (w <= 700); return w[m], a[m] / a[m].max()

V = np.loadtxt(ROOT + 'analysis/basic2025/v1_A.csv', delimiter=','); x = V[:, 0]
T = cuvettes()
fig = plt.figure(figsize=(7.2, 7.4))
gs = fig.add_gridspec(6, 4, width_ratios=[0.3, 1.45, 1.0, 1.0], hspace=0.55, wspace=0.22, left=0.0, right=0.985, top=0.93, bottom=0.07)
for r in range(6):
    a = fig.add_subplot(gs[r, 1]); a.imshow(screen(r)); a.set_xticks([]); a.set_yticks([])
    for sp in a.spines.values(): sp.set_visible(True); sp.set_linewidth(0.5)
    a.set_title('No sample (reference)' if r == 0 else f'Sample {r} ({NAMES[r-1]})', loc='left', fontsize=7.5, pad=2)
    if r == 0: continue
    c = fig.add_subplot(gs[r, 0]); c.imshow(T[r - 1]); c.set_axis_off()
    p = fig.add_subplot(gs[r, 2]); A = V[:, r]; p.plot(x, A / np.nanmax(A), color=COL[r - 1], lw=1.3)
    p.axhline(0, color='0.75', lw=0.6, ls=':'); p.set_xlim(0, x.max()); p.set_ylim(-0.4, 1.15); p.set_yticks([0, 1])
    q = fig.add_subplot(gs[r, 3]); w, u = lab(r); q.plot(w, u, color=COL[r - 1], lw=1.3); q.set_xlim(400, 700); q.set_ylim(-0.1, 1.15); q.set_yticks([0, 1])
    if r < 5: p.set_xticklabels([]); q.set_xticklabels([])
    if r == 3: p.set_ylabel('Normalized absorbance'); q.set_ylabel('Normalized absorbance')
p.set_xlabel('Position along the spectrum (pixel,\nblue → red)'); q.set_xlabel('Wavelength (nm)')
fig.text(0.02, 0.985, '(a) Samples and the projected spectrum', fontsize=8.5, va='top')
x2 = gs[0, 2].get_position(fig).x0; x3 = gs[0, 3].get_position(fig).x0
fig.text(x2 - 0.04, 0.985, '(b) Basic version\n(phone, pixel axis)', fontsize=8.5, va='top')
fig.text(x3 - 0.03, 0.985, '(c) Laboratory\nUV–vis', fontsize=8.5, va='top')
fig.savefig(ROOT + 'figs/S_basic2025.png', dpi=300); fig.savefig(ROOT + 'figs/S_basic2025.pdf')
peaks = {r: float(x[np.nanargmax(V[:, r])]) for r in range(1, 6)}; print('basic-version peak positions (px):', peaks)
