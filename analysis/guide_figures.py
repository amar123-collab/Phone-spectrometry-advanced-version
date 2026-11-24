"""Pictures for docs/PYTHON_GUIDE.md: one figure per step of the analysis, made from the real 2026 photographs.
Run after analyse_paperB.py (needs results.json and arrays.npz):   python guide_figures.py
Writes docs/guide/step01_photos.png ... step10_corrections.png"""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))) + '/'   # repository root
import json, sys
import numpy as np
from PIL import Image, ImageOps
import exifread
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, ROOT + 'uvvis'); sys.path.insert(0, ROOT + 'analysis')
from analyse import load as load_uv
import baselines as BL

PH = ROOT + 'photos/raw/ern_photos/20261120_'
OUT = ROOT + 'docs/guide/'; _os.makedirs(OUT, exist_ok=True)
R = json.load(open(ROOT + 'analysis/results.json')); Z = np.load(ROOT + 'analysis/arrays.npz')
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
px = Z['px']; ok = Z['ok_lin']; NM = R['cal']['offset'] + R['cal']['nm_per_px'] * px
ORANGE, VIOLET, GREY = '#a63603', '#6d3fb0', '0.55'
CROP = (500, 3000, 2150, 3900)                      # (x0, y0, x1, y1): region around the spectrum in the 3000 x 4000 photos
FILES = {'dark': '002855', 'reference (ethanol)': '003701', 'DCM 49.5 µM': '003349', 'sample 1 (violet)': '003411'}
srgb_lin = lambda a: np.where(a <= 0.04045, a / 12.92, ((a + 0.055) / 1.055) ** 2.4)

def photo(f): return ImageOps.exif_transpose(Image.open(PH + f + '.jpg')).convert('RGB')
def settings(f):
    t = exifread.process_file(open(PH + f + '.jpg', 'rb'), details=False); e = t['EXIF ExposureTime'].values[0]
    return f"1/{round(e.den / e.num)} s, ISO {t['EXIF ISOSpeedRatings']}"
def save(fig, name): fig.savefig(OUT + name, dpi=150, bbox_inches='tight', facecolor='white'); plt.close(fig); print('wrote', name)

# 1 -- the photographs
fig, ax = plt.subplots(1, 5, figsize=(13, 3.2), gridspec_kw=dict(width_ratios=[0.75, 1, 1, 1, 1]))
ax[0].imshow(photo(FILES['reference (ethanol)']).resize((300, 400))); ax[0].set_title('Whole photo (reference)', fontsize=9)
ax[0].add_patch(plt.Rectangle((CROP[0] / 10, CROP[1] / 10), (CROP[2] - CROP[0]) / 10, (CROP[3] - CROP[1]) / 10, fill=False, ec='w', lw=1))
for a, (lab, f) in zip(ax[1:], FILES.items()):
    a.imshow(photo(f).crop(CROP)); a.set_title(f'{lab}\n{f}.jpg, {settings(f)}', fontsize=8.5)
for a in ax: a.set_xticks([]); a.set_yticks([])
save(fig, 'step01_photos.png')

# 2 -- linearizing pixel values
v = np.arange(256) / 255; RESP = np.load(ROOT + 'analysis/response_L.npy')
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
ax[0].plot(v * 255, v, color=GREY, label='encoded value used directly (wrong)'); ax[0].plot(v * 255, srgb_lin(v), color=ORANGE, lw=2, label='sRGB inverse (used)')
ax[0].plot(v[5:251] * 255, RESP[5:251], color='#1d4ed8', ls='--', label='measured camera response (check)')
ax[0].plot([128, 128, 0], [0, srgb_lin(128 / 255), srgb_lin(128 / 255)], ':', color='k', lw=0.8)
ax[0].set_xlabel('Pixel value in the JPG (0–255)'); ax[0].set_ylabel('Relative light intensity'); ax[0].set_xlim(0, 255); ax[0].set_ylim(0, 1)
ax[0].set_title('(a) Pixel value 128 is only 22% of full light, not 50%', fontsize=9, loc='left'); ax[0].legend(frameon=False, fontsize=8, loc='upper left')
crop = np.asarray(photo(FILES['reference (ethanol)']).crop(CROP)).astype(float) / 255
row = crop[:, :, :].sum(2)[(3450 - CROP[1]) - 40:(3450 - CROP[1]) + 40].mean(0) / 3
rowl = srgb_lin(crop)[:, :, :].sum(2)[(3450 - CROP[1]) - 40:(3450 - CROP[1]) + 40].mean(0) / 3
xx = np.arange(CROP[0], CROP[2])
ax[1].plot(xx, row / row.max(), color=GREY, label='encoded'); ax[1].plot(xx, rowl / rowl.max(), color=ORANGE, lw=2, label='linearized')
ax[1].set_xlabel('Pixel column'); ax[1].set_ylabel('Brightness (relative to maximum)'); ax[1].set_ylim(0, 1.05)
ax[1].set_title('(b) Same spectrum: linearizing deepens the contrast', fontsize=9, loc='left'); ax[1].legend(frameon=False, fontsize=8)
fig.tight_layout(); save(fig, 'step02_linearize.png')

# 3 -- finding the strip (same procedure as analyse_paperB.py, on the four end-of-series ethanol references)
refimg = sum(srgb_lin(np.asarray(photo(f)).astype(np.float32) / 255).sum(2) for f in ['003654', '003657', '003659', '003701']) / 4
xs = np.arange(700, 1950); yc = []
for x in xs:
    col = refimg[3300:3650, x]; w = np.clip(col - 0.2 * col.max(), 0, None)
    yc.append(3300 + (np.arange(col.size) * w).sum() / w.sum() if w.sum() > 0 else np.nan)
yc = np.array(yc); good = ~np.isnan(yc) & (refimg[3300:3650, 700:1950].max(0) > 0.2 * refimg.max())
tilt = np.polyfit(xs[good], yc[good], 1); cen = np.polyval(tilt, xs)
fig, ax = plt.subplots(figsize=(10, 4.2))
ax.imshow(np.asarray(photo(FILES['reference (ethanol)']).crop(CROP)), extent=(CROP[0], CROP[2], CROP[3], CROP[1]))
ax.plot(xs[good][::15], yc[good][::15], '.', color='w', ms=3, label='brightness-weighted centre of each column')
ax.plot(xs, cen, color='#facc15', lw=1.2, label=f'straight-line fit (tilt {R["tilt_deg"]:.1f}°)')
ax.plot(xs, cen - 40, color='#facc15', lw=0.8, ls='--', label='strip: ±40 rows around the line')
ax.plot(xs, cen + 40, color='#facc15', lw=0.8, ls='--')
ax.set_xlabel('Pixel column'); ax.set_ylabel('Pixel row'); ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), frameon=False, fontsize=8)
ax.set_title('Step 3: the strip follows the slightly tilted spectrum', fontsize=9, loc='left')
save(fig, 'step03_strip.png')

# 4 -- intensity along the strip
def prof(f):
    a = srgb_lin(np.asarray(photo(f)).astype(np.float32) / 255).sum(2); c = np.round(cen).astype(int)
    return np.array([a[ci - 40:ci + 40, x].mean() for x, ci in zip(xs, c)])
D = prof(FILES['dark']); I0 = prof(FILES['reference (ethanol)']); I = prof(FILES['DCM 49.5 µM'])
fig, ax = plt.subplots(figsize=(8, 3.6))
ax.plot(xs, I0, 'k', lw=1.8, label='reference (ethanol)'); ax.plot(xs, I, color=ORANGE, lw=1.5, label='DCM 49.5 µM')
ax.plot(xs, D, color=GREY, lw=1.2, label='dark (lamp blocked)')
ax.set_xlabel('Pixel column (blue → red)'); ax.set_ylabel('Linear intensity (sum of R+G+B)')
ax.set_title('Step 4: light along the strip, one value per column', fontsize=9, loc='left'); ax.legend(frameon=False, fontsize=8)
save(fig, 'step04_profiles.png')

# 5 -- dark subtraction, averaged reference, usable range
ref = Z['ref_lin']
fig, ax = plt.subplots(figsize=(8, 3.6))
ax.plot(px, ref, 'k', lw=1.8, label='reference: mean of 15 ethanol photos, dark subtracted')
ax.axhline(0.15 * ref.max(), color=ORANGE, ls='--', lw=1, label='15% of the maximum')
ax.axvspan(px[ok].min(), px[ok].max(), color='#fde68a', alpha=0.4, lw=0, label=f'usable range: columns {R["usable_px"][0]}–{R["usable_px"][1]}')
ax.set_xlabel('Pixel column (blue → red)'); ax.set_ylabel('Linear intensity')
ax.set_title('Step 5: only columns with enough light are used', fontsize=9, loc='left'); ax.legend(frameon=False, fontsize=8, loc='upper right')
save(fig, 'step05_usable_range.png')

# 6 -- absorbance from three photographs
A = Z['dcm_lin'][-1]; m = np.where(ok, A.mean(0), np.nan); sd = np.where(ok, A.std(0, ddof=1), np.nan)
fig, ax = plt.subplots(figsize=(8, 3.6))
for k, a in enumerate(A): ax.plot(px, np.where(ok, a, np.nan), color=ORANGE, lw=0.6, alpha=0.6, label='each of the three photos' if k == 0 else None)
ax.plot(px, m, color=ORANGE, lw=2, label='mean'); ax.fill_between(px, m - sd, m + sd, color=ORANGE, alpha=0.2, lw=0, label='± standard deviation')
w0, w1 = R['dcm_window_nm']; c0, c1 = [(x - R['cal']['offset']) / R['cal']['nm_per_px'] for x in (w0, w1)]
ax.axvspan(c0, c1, color='0.85', lw=0, label='band-maximum window (31 columns)')
ax.set_xlabel('Pixel column (blue → red)'); ax.set_ylabel('Absorbance  A = −log₁₀(I / I₀)')
ax.set_title('Step 6: absorbance of 49.5 µM DCM', fontsize=9, loc='left'); ax.legend(frameon=False, fontsize=8)
save(fig, 'step06_absorbance.png')

# 7 -- blanks and detection limit
B = R['blank_lin']; vals = B['values']
ETH_A = ['002937', '002950', '003005', '003012']
ETH_B = ['003045', '003051', '003059', '003104', '003109', '003116', '003122', '003126', '003131', '003134', '003138']
fig, ax = plt.subplots(figsize=(8, 3.4))
ax.plot(range(4), [vals[f] for f in ETH_A], 'o', color=GREY, label=f'not moved: SD {B["sd_notmoved"]:.3f}')
ax.plot(range(5, 16), [vals[f] for f in ETH_B], 'o', color=VIOLET, label=f'taken out and put back: SD {B["sd_replaced"]:.3f}')
ax.axhline(0, color='0.8', lw=0.8)
ax.set_xticks([]); ax.set_ylabel('Absorbance of ethanol\nat the DCM window')
ax.set_title(f'Step 7: blank scatter → detection limit 3·SD/slope = A {R["lod_A"]:.2f} ({R["lod_uM"]:.1f} µM DCM)', fontsize=9, loc='left')
ax.legend(frameon=False, fontsize=8, loc='upper left', bbox_to_anchor=(1.01, 1))
save(fig, 'step07_blanks_lod.png')

# 8 -- wavelength calibration
b = R['bands']; names = ['DCM', 'ERN4', 'ERN3', 'ERN1']; lab = {'DCM': 'DCM', 'ERN1': 'sample 1', 'ERN3': 'sample 3', 'ERN4': 'sample 4'}
fig, ax = plt.subplots(figsize=(6.5, 3.8))
xx = np.array(R['usable_px']); ax.plot(xx, R['cal']['nm_per_px'] * xx + R['cal']['offset'], color=GREY)
for n in names: ax.plot(*b[n], 'o', ms=7, label=f'{lab[n]}: column {b[n][0]:.0f} → {b[n][1]:.0f} nm')
ax.set_xlabel('Pixel column of the band maximum (phone)'); ax.set_ylabel('Band maximum, laboratory (nm)')
ax.set_title(f'Step 8: nm = {R["cal"]["nm_per_px"]:.4f} × column + {R["cal"]["offset"]:.1f}  (leave-one-out error {R["cal"]["loo_rms"]:.1f} nm)', fontsize=9, loc='left')
ax.legend(frameon=False, fontsize=8, loc='upper left', bbox_to_anchor=(1.01, 1))
save(fig, 'step08_calibration.png')

# 9 -- comparison with the laboratory: linearity
uvA = np.array(R['dcm_uvA']); d = R['dcm_lin']; p = d['fit_all']
fig, ax = plt.subplots(figsize=(5.2, 4.2))
ax.plot([0, 2], [0, 2], ':', color='0.7', label='1:1')
ax.errorbar(uvA, d['mean'], yerr=d['sd'], fmt='o', color=ORANGE, capsize=0, label=f'phone (slope {p["slope"]:.2f}, R² {p["R2"]:.4f})')
ax.plot([0, 2], np.polyval([p['slope'], p['intercept']], [0, 2]), color=ORANGE, lw=0.8)
ax.set_xlim(0, 2.05); ax.set_ylim(0, 2.05); ax.set_aspect('equal')
ax.set_xlabel('Absorbance, laboratory (465–475 nm)'); ax.set_ylabel('Absorbance, phone (band window)')
ax.set_title('Step 9: six DCM concentrations', fontsize=9, loc='left'); ax.legend(frameon=False, fontsize=8, loc='upper left')
save(fig, 'step09_linearity.png')

# 10 -- corrections for comparing whole spectra (sample 1)
A1 = Z['ern_lin'][0].mean(0); bl, kind = BL.baseline(A1, 'real1')
w, u = load_uv('real1')[:2]
def blur(w, a, fwhm):
    T = 10 ** -a; g = np.arange(380, 800.5, 0.5); Tg = np.interp(g, w, T); s = fwhm / 2.355; k = np.exp(-0.5 * (np.arange(-4 * s, 4 * s + 0.5, 0.5) / s) ** 2)
    return g, -np.log10(np.convolve(Tg, k / k.sum(), 'same'))
g, ub = blur(w, u, 25); mm = (g >= 430) & (g <= 675)
fig, ax = plt.subplots(1, 2, figsize=(11, 3.6), sharex=True)
ax[0].plot(NM, np.where(ok, A1, np.nan), color=VIOLET, lw=1.8, label='phone, as measured')
ax[0].plot(NM, np.where(ok, bl, np.nan), color='k', ls='--', lw=1, label=f'baseline ({kind} through the transparent parts)')
ax[0].set_title('(a) The cuvette adds a sloping background', fontsize=9, loc='left')
ax[1].plot(g[mm], ub[mm] / ub[mm].max(), color='k', lw=1, label='laboratory, broadened to 25 nm')
c = np.where(ok, A1 - bl, np.nan); ax[1].plot(NM, c / np.nanmax(c), color=VIOLET, lw=1.8, label='phone, baseline removed')
ax[1].set_title('(b) After the correction, the shapes agree', fontsize=9, loc='left')
for a in ax: a.set_xlabel('Wavelength (nm)'); a.legend(frameon=False, fontsize=8, loc='upper left'); a.set_xlim(425, 680)
ax[0].set_ylabel('Absorbance'); ax[1].set_ylabel('Normalized absorbance')
fig.tight_layout(); save(fig, 'step10_corrections.png')
