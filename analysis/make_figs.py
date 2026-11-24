"""Figures for Paper B (main text + supplement) from results.json / arrays.npz and the photos.
Style: all text outside the plot areas, colour progression for concentration series."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))) + '/'   # repository root
import json, sys
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image, ImageOps
sys.path.insert(0, ROOT + 'uvvis'); from analyse import load as load_uv
import exifread

OUT = ROOT + 'figs/'
import os; os.makedirs(OUT, exist_ok=True)
R = json.load(open(ROOT + 'analysis/results.json'))
Z = np.load(ROOT + 'analysis/arrays.npz')
NM, OK, PX, CEN = Z['nm'], Z['ok_lin'], Z['px'], Z['centre']
PH = ROOT + 'photos/raw/ern_photos/'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8.5, 'axes.linewidth': 0.8,
                     'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False})
ORANGES = plt.cm.Oranges(np.linspace(0.35, 0.95, 6))
C_ENC, C_LIN, C_MEAS, C_UV = '#9a9a9a', '#c2410c', '#1d4ed8', 'black'
ERN_COL = ['#6d3fb0', '#1f78b4', '#c2185b', '#e0457b', '#9aa800']

def photo(name, crop=None):
    im = ImageOps.exif_transpose(Image.open(PH + name + '.jpg')).convert('RGB')
    return np.asarray(im.crop(crop) if crop else im)
def save(fig, name):
    fig.savefig(OUT + name + '.png', dpi=300, bbox_inches='tight'); fig.savefig(OUT + name + '.pdf', bbox_inches='tight'); plt.close(fig)
def masked(a):
    a = a.copy(); a[~OK] = np.nan; return a
CROP = (600, 3230, 2050, 3700)
uvA = np.array(R['dcm_uvA']); conc = np.array(R['dcm_conc_uM'])

# ---------- Fig 2: samples ----------
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.3), gridspec_kw={'width_ratios': [1, 1]})
for a, n, t in zip(ax, ['20261119_231654', '20261119_233828'], ['(a) Ethanol and the DCM series', '(b) Mixed-dye event samples 1–5']):
    a.imshow(photo(n)); a.set_axis_off(); a.set_title(t, loc='left', fontsize=8.5)
fig.tight_layout(); save(fig, 'Fig2_samples')

# ---------- Fig 3: from photo to absorbance ----------
fig = plt.figure(figsize=(7.2, 5.0))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.25])
for k, (n, t) in enumerate([('20261120_003701', '(a) Reference (ethanol)'), ('20261120_003349', '(b) DCM, 49.5 µM')]):
    a = fig.add_subplot(gs[0, k]); a.imshow(photo(n, CROP)); a.set_axis_off(); a.set_title(t, loc='left', fontsize=8.5)
    xx = PX - CROP[0]; a.plot(xx, CEN - 40 - CROP[1], color='white', lw=0.6, ls='--'); a.plot(xx, CEN + 40 - CROP[1], color='white', lw=0.6, ls='--')
a = fig.add_subplot(gs[1, 0])
ref = Z['ref_lin']; a.plot(NM, ref / ref.max(), color='0.2', lw=1.4, label='reference')
smp = ref * 10 ** (-Z['dcm_lin'][-1].mean(0)); a.plot(NM, smp / ref.max(), color=ORANGES[-1], lw=1.4, label='DCM 49.5 µM')
a.axvspan(*R['usable_nm'], color='0.93', zorder=0)
a.set_xlabel('Wavelength (nm)'); a.set_ylabel('Linear intensity (norm.)'); a.set_title('(c) Strip profiles', loc='left', fontsize=8.5)
a.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), fontsize=7.5); a.set_xlim(420, 690)
a = fig.add_subplot(gs[1, 1])
w, u, _ = load_uv('a5'); m = (w >= 420) & (w <= 690)
a.plot(w[m], u[m], color=C_UV, lw=1.2, label='UV-Vis'); a.plot(NM, masked(Z['dcm_lin'][-1].mean(0)), color=C_LIN, lw=1.4, label='phone, sRGB-linearized')
a.plot(NM, masked(Z['dcm_enc'][-1].mean(0)), color=C_ENC, lw=1.4, label='phone, encoded')
a.set_xlabel('Wavelength (nm)'); a.set_ylabel('Absorbance'); a.set_title('(d) DCM 49.5 µM', loc='left', fontsize=8.5); a.set_xlim(420, 690)
a.legend(loc='upper left', bbox_to_anchor=(0.0, -0.22), ncol=3, fontsize=7.5)
fig.tight_layout(); save(fig, 'Fig3_photo_to_absorbance')

# ---------- Fig 4: wavelength calibration from sample bands ----------
cal = R['cal']; B = R['bands']
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.7), gridspec_kw={'width_ratios': [1.3, 1]})
names = ['DCM', 'ERN4', 'ERN3', 'ERN1']; lab = {'DCM': 'DCM', 'ERN1': 'Sample 1', 'ERN3': 'Sample 3', 'ERN4': 'Sample 4'}
cols = {'DCM': ORANGES[-1], 'ERN1': ERN_COL[0], 'ERN3': ERN_COL[2], 'ERN4': ERN_COL[3]}
xx = np.linspace(R['usable_px'][0], R['usable_px'][1], 50)
ax[0].plot(xx, cal['nm_per_px'] * xx + cal['offset'], color='0.5', lw=1)
for n in names: ax[0].plot(B[n][0], B[n][1], 'o', ms=7, color=cols[n], mec='white', label=lab[n])
for n, c in (('ERN2_edge', ERN_COL[1]), ('ERN5_edge', ERN_COL[4])):
    ax[0].plot(B[n][0], B[n][1], 'o', ms=7, mfc='white', mec=c, label=('Sample 2' if '2' in n else 'Sample 5') + ' (at strip edge, not used)')
ax[0].set_xlabel('Position along the strip (pixel)'); ax[0].set_ylabel('UV-Vis band maximum (nm)'); ax[0].set_title('(a) Reference-based calibration', loc='left', fontsize=8.5)
ax[0].legend(loc='upper left', bbox_to_anchor=(0.0, -0.22), ncol=2, fontsize=7)
loo = cal['loo']; ax[1].bar(range(4), [loo[n] for n in names], color=[cols[n] for n in names], width=0.6)
ax[1].axhline(0, color='0.3', lw=0.8); ax[1].set_xticks(range(4)); ax[1].set_xticklabels([lab[n] for n in names], fontsize=7.5)
ax[1].set_ylabel('Leave-one-out error (nm)'); ax[1].set_title('(b) Cross-validation', loc='left', fontsize=8.5); ax[1].set_ylim(-15, 15)
fig.tight_layout(); save(fig, 'Fig4_calibration')

# ---------- Fig 5: linearity ----------
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={'width_ratios': [1.25, 1]})
for i in range(6):
    w, u, _ = load_uv(['b1', 'b2', 'a2', 'a3', 'a4', 'a5'][i]); m = (w >= 425) & (w <= 680)
    ax[0].plot(w[m], u[m], color=ORANGES[i], lw=0.9, ls='--')
    ax[0].plot(NM, masked(Z['dcm_lin'][i].mean(0)), color=ORANGES[i], lw=1.5, label=f'{conc[i]:.1f} µM')
ax[0].set_xlim(425, 680); ax[0].set_xlabel('Wavelength (nm)'); ax[0].set_ylabel('Absorbance')
ax[0].set_title('(a) Phone (solid) and UV-Vis (dashed)', loc='left', fontsize=8.5)
ax[0].legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), fontsize=7.5)
lim = [0, 2.05]; ax[1].plot(lim, lim, color='0.75', lw=0.8, ls=':')
for tag, c, mk, lab in (('lin', C_LIN, 'o', 'sRGB-linearized'), ('meas', C_MEAS, 's', 'measured response'), ('enc', C_ENC, 'D', 'encoded')):
    d = R['dcm_' + tag]; ax[1].errorbar(uvA, d['mean'], yerr=d['sd'], xerr=np.array(R['dcm_conc_err_uM']) * R['uv_beer']['eps'] / 1e6, fmt=mk, ms=5, color=c, mec='white', mew=0.6, elinewidth=1, capsize=0,
                                       label=f"{lab} (slope {d['fit_all']['slope']:.2f})")
    p = d['fit_all']; ax[1].plot(lim, np.polyval([p['slope'], p['intercept']], lim), color=c, lw=0.8)
ax[1].set_xlim(lim); ax[1].set_ylim(lim); ax[1].set_aspect('equal')
ax[1].set_xlabel('Absorbance, UV-Vis (465–475 nm)'); ax[1].set_ylabel('Absorbance, phone'); ax[1].set_title('(b) Linearity', loc='left', fontsize=8.5)
ax[1].legend(loc='upper left', bbox_to_anchor=(0.0, -0.2), fontsize=7.2)
fig.tight_layout(); save(fig, 'Fig5_linearity')

# ---------- Fig 6: camera response ----------
sh = R['shutter']; e = np.array([s['exposure_s'] for s in sh]); sat = np.array([s['sat_pct'] for s in sh])
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
for tag, c, mk, lab in (('lin', C_LIN, 'o', 'sRGB-linearized'), ('meas', C_MEAS, 's', 'measured response'), ('enc', C_ENC, 'D', 'encoded')):
    y = np.array([s['signal_' + tag] for s in sh]); y = y / y[0] * e[0] / e[0]
    ok_ = sat < 0.5
    ax[0].loglog(e[ok_], y[ok_], mk + '-', color=c, mec='white', ms=6, lw=1, label=f"{lab} (slope {R['shutter_loglog_slope_' + tag]:.2f})")
    ax[0].loglog(e[~ok_], y[~ok_], mk, color=c, mfc='none', ms=6)
ee = np.array([e[0], e[2]]); ax[0].loglog(ee, ee / e[0], color='0.75', lw=0.8, ls=':')
ax[0].set_xlabel('Shutter time (s)'); ax[0].set_ylabel('Signal relative to 1/500 s'); ax[0].set_title('(a) Shutter series (open: clipped)', loc='left', fontsize=8.5)
ax[0].legend(loc='upper left', bbox_to_anchor=(0.0, -0.22), fontsize=7.2)
L = np.load(ROOT + 'analysis/response_L.npy'); v = np.arange(256) / 255
srgb = np.where(v <= 0.04045, v / 12.92, ((v + 0.055) / 1.055) ** 2.4)
ax[1].plot(v * 255, v, color=C_ENC, lw=1.3, label='encoded (no correction)'); ax[1].plot(v * 255, srgb, color=C_LIN, lw=1.5, label='sRGB inverse')
ax[1].plot(v[5:251] * 255, L[5:251], color=C_MEAS, lw=1.5, label='measured, Galaxy S25 FE')
ax[1].set_xlabel('Pixel value (0–255)'); ax[1].set_ylabel('Relative linear intensity'); ax[1].set_title('(b) Inverse tone curves', loc='left', fontsize=8.5)
ax[1].legend(loc='upper left', bbox_to_anchor=(0.0, -0.22), fontsize=7.2); ax[1].set_xlim(0, 255); ax[1].set_ylim(0, 1)
fig.tight_layout(); save(fig, 'Fig6_camera_response')

# ---------- Fig 7: event samples ----------
fig, axs = plt.subplots(1, 5, figsize=(7.2, 2.1), sharey=True)
uv = ['real1', 'real2', 'real3', 'real4', 'real5']
for i, a in enumerate(axs):
    w, u, _ = load_uv(uv[i]); m = (w >= 425) & (w <= 690); um = u[m].max()
    ph = masked(Z['ern_lin'][i].mean(0)); sd = masked(Z['ern_lin'][i].std(0, ddof=1))
    a.plot(w[m], u[m] / um, color=C_UV, lw=1); a.plot(NM, ph / um, color=ERN_COL[i], lw=1.4)
    a.fill_between(NM, (ph - sd) / um, (ph + sd) / um, color=ERN_COL[i], alpha=0.25, lw=0)
    a.axvspan(425, R['usable_nm'][0], color='0.92', zorder=0); a.axvspan(R['usable_nm'][1], 690, color='0.92', zorder=0)
    a.set_xlim(425, 690); a.set_title(f'Sample {i+1}', fontsize=8); a.set_xticks([450, 550, 650])
    if i == 0: a.set_ylabel('Absorbance / UV-Vis max.')
fig.supxlabel('Wavelength (nm)', fontsize=8.5, y=0.02)
fig.tight_layout(); save(fig, 'Fig7_event_samples')

# ================= Supplement =================
# S: clipping test photos and profiles
fig = plt.figure(figsize=(7.2, 4.6)); gs = fig.add_gridspec(2, 3)
EX = [('20261120_001159', '1/50 s'), ('20261120_002033', '1/100 s'), ('20261120_002315', '1/125 s')]
for k, (n, t) in enumerate(EX):
    a = fig.add_subplot(gs[0, k]); im = photo(n); I = im.sum(2); ys, xs_ = np.nonzero(I > 0.2 * I.max())
    c = (max(xs_.min() - 150, 0), max(ys.min() - 120, 0), xs_.max() + 150, ys.max() + 120)
    a.imshow(im[c[1]:c[3], c[0]:c[2]]); a.set_axis_off()
    a.set_title(f"({'abc'[k]}) {t}: {R['clipping'][n[9:]]['clipped_pct_of_spectrum']:.1f} % clipped", loc='left', fontsize=8)
a = fig.add_subplot(gs[1, :])
for (n, t), c in zip(EX, ['#7f1d1d', '#c2410c', '#fb923c']):
    im = photo(n).astype(float); I = im.sum(2); ys, xs_ = np.nonzero(I > 0.2 * I.max()); yc = int(np.median(ys))
    for ch, ls in zip(range(3), ['-', '--', ':']):
        p = im[yc - 20:yc + 20, :, ch].mean(0); sel = slice(xs_.min() - 50, xs_.max() + 50)
        a.plot(np.arange(sel.stop - sel.start), p[sel], color=c, ls=ls, lw=1, label=f'{t}, {"RGB"[ch]}')
a.axhline(255, color='0.6', lw=0.6); a.set_xlabel('Position along the spectrum (pixel, from blue end)'); a.set_ylabel('JPG value')
a.set_title('(d) Channel profiles; flat tops at 255 are clipped', loc='left', fontsize=8)
a.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), fontsize=6.5)
fig.tight_layout(); save(fig, 'S_clipping_test')

# S: auto mode
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.3))
for a, (n, t) in zip(ax, [('20261120_011859', 'Auto mode (ISO 1600, 1/50 s)'), ('20261120_011830', 'Pro mode (ISO 100, 1/90 s)')]):
    im = photo(n); a.imshow(im[3350:3850, 700:1950]); a.set_axis_off()
    a.set_title(f"{t}: {R['clipping'].get(n[9:], {}).get('clipped_pct_of_spectrum', float('nan')):.0f} % clipped" if n[9:] in R['clipping'] else t, loc='left', fontsize=8)
fig.tight_layout(); save(fig, 'S_auto_mode')

# S: shutter series photos
fig, axs = plt.subplots(2, 3, figsize=(7.2, 3.0))
for a, s in zip(axs.flat, R['shutter']):
    a.imshow(photo('20261120_' + s['file'])[3350:3850, 700:1950]); a.set_axis_off()
    a.set_title(f"1/{1/s['exposure_s']:.0f} s, {s['sat_pct']:.0f} % clipped" if s['exposure_s'] < 1 else f"{s['exposure_s']} s", loc='left', fontsize=7.5)
fig.tight_layout(); save(fig, 'S_shutter_series')

# S: blanks, lamp drift, no cuvette
from datetime import datetime
t0 = datetime.strptime('00:29:37', '%H:%M:%S')
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.6))
names = list(R['blank_lin']['values'].keys()); T = R['ref_times']; tot = R['ref_total_lin']
tt = [(datetime.strptime(T[n], '%H:%M:%S') - t0).total_seconds() / 60 for n in names]
grp = ['not moved'] * 4 + ['taken out and put back'] * 11 + ['end of series'] * 4 + ['no cuvette?']
gc = {'not moved': '#2563eb', 'taken out and put back': '#c2410c', 'end of series': '#15803d', 'no cuvette?': '0.5'}
for g in gc:
    idx = [i for i, x in enumerate(grp) if x == g]
    ax[0].plot([tt[i] for i in idx], [tot[names[i]] / tot[names[0]] for i in idx], 'o', color=gc[g], ms=5, mec='white', label=g)
    ax[1].plot([tt[i] for i in idx], [R['blank_lin']['values'][names[i]] for i in idx], 'o', color=gc[g], ms=5, mec='white')
ax[0].set_xlabel('Time (min)'); ax[0].set_ylabel('Reference signal (rel.)'); ax[0].set_title('(a) Lamp and placement', loc='left', fontsize=8)
ax[1].axhline(0, color='0.6', lw=0.6); ax[1].set_xlabel('Time (min)'); ax[1].set_ylabel('Blank absorbance (466–473 nm)'); ax[1].set_title('(b) Blanks', loc='left', fontsize=8)
ax[0].legend(loc='upper left', bbox_to_anchor=(0.0, -0.25), ncol=2, fontsize=6.8)
for i in range(4, 15): ax[2].plot(NM, masked(Z['blank_lin'][i]), color=gc['taken out and put back'], lw=0.6, alpha=0.7)
for i in range(4): ax[2].plot(NM, masked(Z['blank_lin'][i]), color=gc['not moved'], lw=0.6)
ax[2].axhline(0, color='0.6', lw=0.6); ax[2].set_ylim(-0.3, 0.3); ax[2].set_xlabel('Wavelength (nm)'); ax[2].set_ylabel('Blank absorbance')
ax[2].set_title('(c) Blank spectra', loc='left', fontsize=8)
fig.tight_layout(); save(fig, 'S_blanks_drift')

# S: UV-Vis processing (raw vs processed) for DCM 49.5 and sample 2
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.6))
for a, f, t in zip(ax, ['a5', 'real2'], ['DCM 49.5 µM', 'Sample 2']):
    d = np.genfromtxt(ROOT + 'uvvis/' + f + '.csv', delimiter=',', skip_header=1)
    w, u, b = load_uv(f); m = (w >= 380) & (w <= 800)
    a.plot(d[:, 0][m], d[:, 1][m], color='0.6', lw=1, label='as recorded'); a.plot(w[m], u[m], color='black', lw=1.1, label='baseline-corrected, artefacts patched')
    for lo, hi in [(484, 488), (630, 641), (653, 659)]: a.axvspan(lo, hi, color='#fde68a', zorder=0)
    a.set_xlabel('Wavelength (nm)'); a.set_ylabel('Absorbance'); a.set_title(t, loc='left', fontsize=8)
ax[0].legend(loc='upper left', bbox_to_anchor=(0.0, -0.22), ncol=2, fontsize=7)
fig.tight_layout(); save(fig, 'S_uvvis_processing')

# S: all UV-Vis spectra + Beer-Lambert
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.6), gridspec_kw={'width_ratios': [1.2, 1, 1.2]})
for i, f in enumerate(['b1', 'b2', 'a2', 'a3', 'a4', 'a5']):
    w, u, _ = load_uv(f); m = (w >= 380) & (w <= 750); ax[0].plot(w[m], u[m], color=ORANGES[i], lw=1.2)
w, u, _ = load_uv('a1'); m = (w >= 380) & (w <= 750); ax[0].plot(w[m], u[m], color='0.6', lw=1)
ax[0].set_xlabel('Wavelength (nm)'); ax[0].set_ylabel('Absorbance'); ax[0].set_title('(a) DCM and ethanol', loc='left', fontsize=8)
ub = R['uv_beer']; ax[1].errorbar(conc, uvA, xerr=R['dcm_conc_err_uM'], fmt='o', color=ORANGES[-1], mec='white', ms=5, elinewidth=1, capsize=0)
xx = np.array([0, 52]); ax[1].plot(xx, ub['eps'] / 1e6 * xx + ub['intercept'], color='0.5', lw=0.8)
ax[1].set_xlabel('Concentration (µM)'); ax[1].set_ylabel('Absorbance (465–475 nm)'); ax[1].set_title('(b) Beer–Lambert', loc='left', fontsize=8)
for i, f in enumerate(uv):
    w, u, _ = load_uv(f); m = (w >= 380) & (w <= 750); ax[2].plot(w[m], u[m], color=ERN_COL[i], lw=1.2, label=f'Sample {i+1}')
ax[2].set_xlabel('Wavelength (nm)'); ax[2].set_ylabel('Absorbance'); ax[2].set_title('(c) Event samples', loc='left', fontsize=8)
ax[2].legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), fontsize=7)
fig.tight_layout(); save(fig, 'S_uvvis_all')

# S: laser photos
fig, axs = plt.subplots(1, 5, figsize=(7.2, 1.9))
for a, (n, t) in zip(axs, [('20261120_003701', 'Spectrum (00:37)'), ('20261120_005842', '405 nm'), ('20261120_010040', '532 nm'), ('20261120_011358', '650 nm'), ('20261120_011839', 'Spectrum (01:18)')]):
    a.imshow(photo(n)[3250:3900, 600:2050]); a.set_axis_off(); a.set_title(t, loc='left', fontsize=7.5)
fig.tight_layout(); save(fig, 'S_lasers')

# S: replicate phone spectra for every DCM concentration (3 photos each)
fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.0), sharex=True)
for i, a in enumerate(axs.flat):
    for r in range(3): a.plot(NM, masked(Z['dcm_lin'][i][r]), color=ORANGES[i], lw=0.9, alpha=0.9)
    w, u, _ = load_uv(['b1', 'b2', 'a2', 'a3', 'a4', 'a5'][i]); m = (w >= 425) & (w <= 680); a.plot(w[m], u[m], color='black', lw=0.9, ls='--')
    a.set_title(f'{conc[i]:.1f} µM', loc='left', fontsize=8); a.set_xlim(425, 680)
fig.supxlabel('Wavelength (nm)', fontsize=8.5); fig.supylabel('Absorbance', fontsize=8.5)
fig.tight_layout(); save(fig, 'S_dcm_replicates')

# S: power supply / setup photo
fig, ax = plt.subplots(figsize=(3.0, 2.3)); ax.imshow(photo('20261119_181007')); ax.set_axis_off(); fig.tight_layout(); save(fig, 'S_power_supply')
print('figures written:', sorted(os.listdir(OUT)))
