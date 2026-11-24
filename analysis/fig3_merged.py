"""Main Figure 3: gamma encoding (a,b), camera response (c,d), wavelength calibration (e,f) in one figure.
Replaces the former Figures 4, 6 and 7."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))) + '/'   # repository root
import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
A_DIR = ROOT + 'analysis/'; OUT = ROOT + 'figs/'
R = json.load(open(A_DIR + 'results.json')); G = json.load(open(A_DIR + 'gamma_2026.json')); GP = np.load(A_DIR + 'gamma_2026.npz')
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8, 'axes.linewidth': 0.8, 'axes.spines.top': False,
                     'axes.spines.right': False, 'legend.frameon': False})
C_ENC, C_LIN, C_MEAS = '#9a9a9a', '#c2410c', '#1d4ed8'
ORANGES = plt.cm.Oranges(np.linspace(0.35, 0.95, 6)); ERN_COL = ['#6d3fb0', '#1f78b4', '#c2185b', '#e0457b', '#9aa800']
enc = lambda x: np.where(x <= 0.0031308, 12.92 * x, 1.055 * x ** (1 / 2.4) - 0.055)
LEG = dict(loc='upper left', bbox_to_anchor=(0.0, -0.24), fontsize=6.8, handlelength=1.8)

fig, ax = plt.subplots(3, 2, figsize=(7.2, 9.4))
# (a) calculated effect of encoding
At = np.linspace(0, 2.5, 200); greys = plt.cm.Greys(np.linspace(0.95, 0.45, 4))
for f, c, ls in zip([1, 0.3, 0.1, 0.02], greys, ['-', '--', '-.', ':']):
    ax[0, 0].plot(At, -np.log10(enc(f * 10 ** -At) / enc(f)), color=c, ls=ls, lw=1.3, label=f'reference at {f:g} of full scale')
ax[0, 0].plot([0, 2.5], [0, 2.5], color='0.8', lw=0.8)
ax[0, 0].set_xlim(0, 2.5); ax[0, 0].set_ylim(0, 2.5); ax[0, 0].set_xlabel('True absorbance'); ax[0, 0].set_ylabel('Absorbance from encoded values')
ax[0, 0].set_title('(a) Effect of encoding (calculated)', loc='left', fontsize=8.5); ax[0, 0].legend(ncol=2, **LEG)
# (b) 2026 set, Pro mode
ax[0, 1].plot(GP['xD'], GP['yD'], 'o', ms=2, color=ORANGES[-1], alpha=0.7, mew=0, label=f'DCM 49.5 µM ({G["dcm_49uM"]:.2f})')
for k in range(1, 6):
    ax[0, 1].plot(GP[f'x{k}'], GP[f'y{k}'], 'o', ms=2, color=ERN_COL[k - 1], alpha=0.7, mew=0, label=f'Sample {k} ({G["samples"][str(k)]:.2f})')
xx = np.array([0, 2.0]); ax[0, 1].plot(xx, xx, color='0.8', lw=0.8); ax[0, 1].plot(xx, G['slope_all'] * xx, 'k--', lw=0.9, label=f'samples 1–5 ({G["slope_all"]:.2f})')
ax[0, 1].set_xlim(-0.1, 2.0); ax[0, 1].set_ylim(-0.1, 2.0)
ax[0, 1].set_xlabel('Absorbance from linearized values'); ax[0, 1].set_ylabel('Absorbance from encoded values')
ax[0, 1].set_title('(b) Measured: 2026 set and DCM (Pro mode)', loc='left', fontsize=8.5); ax[0, 1].legend(ncol=2, markerscale=2.5, **LEG)
# (c) shutter series
sh = R['shutter']; e = np.array([s['exposure_s'] for s in sh]); sat = np.array([s['sat_pct'] for s in sh]); okc = sat < 0.5
for tag, c, mk, lab in (('lin', C_LIN, 'o', 'sRGB-linearized'), ('meas', C_MEAS, 's', 'measured response'), ('enc', C_ENC, 'D', 'encoded')):
    y = np.array([s['signal_' + tag] for s in sh]); y = y / y[0]
    ax[1, 0].loglog(e[okc], y[okc], mk + '-', color=c, mec='white', ms=5, lw=1, label=f"{lab} (slope {R['shutter_loglog_slope_' + tag]:.2f})")
    ax[1, 0].loglog(e[~okc], y[~okc], mk, color=c, mfc='none', ms=5)
ee = np.array([e[0], e[2]]); ax[1, 0].loglog(ee, ee / e[0], color='0.75', lw=0.8, ls=':')
ax[1, 0].set_xlabel('Shutter time (s)'); ax[1, 0].set_ylabel('Signal relative to 1/500 s')
ax[1, 0].set_title('(c) Shutter series (open symbols: clipped)', loc='left', fontsize=8.5); ax[1, 0].legend(**LEG)
# (d) inverse tone curves
L = np.load(A_DIR + 'response_L.npy'); v = np.arange(256) / 255
srgb = np.where(v <= 0.04045, v / 12.92, ((v + 0.055) / 1.055) ** 2.4)
ax[1, 1].plot(v * 255, v, color=C_ENC, lw=1.3, label='encoded (no correction)'); ax[1, 1].plot(v * 255, srgb, color=C_LIN, lw=1.5, label='sRGB inverse')
ax[1, 1].plot(v[5:251] * 255, L[5:251], color=C_MEAS, lw=1.5, label='measured, Galaxy S25 FE')
ax[1, 1].set_xlabel('Pixel value (0–255)'); ax[1, 1].set_ylabel('Relative linear intensity'); ax[1, 1].set_xlim(0, 255); ax[1, 1].set_ylim(0, 1)
ax[1, 1].set_title('(d) Inverse tone curves', loc='left', fontsize=8.5); ax[1, 1].legend(**LEG)
# (e) wavelength calibration, (f) cross-validation
cal = R['cal']; B = R['bands']; names = ['DCM', 'ERN4', 'ERN3', 'ERN1']
lab = {'DCM': 'DCM', 'ERN1': 'Sample 1', 'ERN3': 'Sample 3', 'ERN4': 'Sample 4'}
cols = {'DCM': ORANGES[-1], 'ERN1': ERN_COL[0], 'ERN3': ERN_COL[2], 'ERN4': ERN_COL[3]}
xx = np.linspace(R['usable_px'][0], R['usable_px'][1], 50); ax[2, 0].plot(xx, cal['nm_per_px'] * xx + cal['offset'], color='0.5', lw=1)
for n in names: ax[2, 0].plot(B[n][0], B[n][1], 'o', ms=6, color=cols[n], mec='white', label=lab[n])
for n, c in (('ERN2_edge', ERN_COL[1]), ('ERN5_edge', ERN_COL[4])):
    ax[2, 0].plot(B[n][0], B[n][1], 'o', ms=6, mfc='white', mec=c, label=('Sample 2' if '2' in n else 'Sample 5') + ' (strip edge, not used)')
ax[2, 0].set_xlabel('Position along the strip (pixel)'); ax[2, 0].set_ylabel('Laboratory band maximum (nm)')
ax[2, 0].set_title('(e) Wavelength calibration', loc='left', fontsize=8.5); ax[2, 0].legend(ncol=2, **LEG)
loo = cal['loo']; ax[2, 1].bar(range(4), [loo[n] for n in names], color=[cols[n] for n in names], width=0.6)
ax[2, 1].axhline(0, color='0.3', lw=0.8); ax[2, 1].set_xticks(range(4)); ax[2, 1].set_xticklabels([lab[n] for n in names], fontsize=7.5)
ax[2, 1].set_ylabel('Leave-one-out error (nm)'); ax[2, 1].set_ylim(-15, 15)
ax[2, 1].set_title(f'(f) Cross-validation (RMS {cal["loo_rms"]:.1f} nm)', loc='left', fontsize=8.5)
fig.subplots_adjust(left=0.09, right=0.98, top=0.97, bottom=0.1, hspace=0.8, wspace=0.3)
fig.savefig(OUT + 'Fig3_gamma_response_calibration.png', dpi=300); fig.savefig(OUT + 'Fig3_gamma_response_calibration.pdf')
