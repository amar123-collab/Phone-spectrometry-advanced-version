"""Figure 5 in the style of Figure 8: one panel per DCM concentration, commercial UV-vis (black) vs home-built prism + phone (colored),
plus the linearity plot below."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))) + '/'   # repository root
import json, sys, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, ROOT + 'uvvis')
from analyse import load as load_uv
sys.path.insert(0, ROOT + 'analysis'); import baselines as BL
A_DIR = ROOT + 'analysis/'; FIG = ROOT + 'figs/'
Z = np.load(A_DIR + 'arrays.npz'); R = json.load(open(A_DIR + 'results.json'))
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8.5, 'axes.linewidth': 0.8,
                     'axes.spines.top': False, 'axes.spines.right': False})
C_ENC, C_LIN, C_MEAS = '#9a9a9a', '#c2410c', '#1d4ed8'
ORANGES = plt.cm.Oranges(np.linspace(0.35, 0.95, 6))
px = Z['px']; ok = Z['ok_lin']; NM = R['cal']['offset'] + R['cal']['nm_per_px'] * px
conc = R['dcm_conc_uM']; uvA = np.array(R['dcm_uvA'])
mask = lambda y: np.where(ok, y, np.nan)
UVF = ['b1', 'b2', 'a2', 'a3', 'a4', 'a5']

from cuvette_panel import tiles
TC = tiles('20261119_231654')[1:]                       # the six DCM cuvettes (ethanol reference omitted)
fig = plt.figure(figsize=(7.2, 6.6))
gs = fig.add_gridspec(3, 6, height_ratios=[0.62, 0.85, 1.45], hspace=0.4, wspace=0.12, left=0.09, right=0.98, top=0.95, bottom=0.08)
r_vals = []; a0 = None
for i in range(6):
    c = fig.add_subplot(gs[0, i]); c.imshow(TC[i]); c.set_axis_off()
    a = fig.add_subplot(gs[1, i], sharey=a0); a0 = a0 or a
    w, u, _ = load_uv(UVF[i]); m = (w >= 425) & (w <= 690)
    A = Z['dcm_lin'][i]; A = np.array([x - BL.baseline(x, UVF[i], dcm=True)[0] for x in A]); ph = mask(A.mean(0)); sd = mask(A.std(0, ddof=1))
    r = np.corrcoef(np.interp(NM[ok], w, u), A.mean(0)[ok])[0, 1]; r_vals.append(float(r))
    a.plot(w[m], u[m], color='black', lw=1)
    a.plot(NM, ph, color=ORANGES[i], lw=1.4); a.fill_between(NM, ph - sd, ph + sd, color=ORANGES[i], alpha=0.3, lw=0)
    a.axvspan(425, R['usable_nm'][0], color='0.92', zorder=0); a.axvspan(R['usable_nm'][1], 690, color='0.92', zorder=0)
    a.set_xlim(425, 690); a.set_xticks([450, 550, 650]); a.tick_params(labelsize=7)
    a.set_ylim(-0.15, 2.1)
    a.set_title(f'{conc[i]:.1f} µM\nr = {r:.3f}', fontsize=7.5)
    if i == 0: a.set_ylabel('Absorbance')
    else: plt.setp(a.get_yticklabels(), visible=False)
ya = gs[1, 0].get_position(fig).y0
fig.text(0.535, ya - 0.055, 'Wavelength (nm)', ha='center', fontsize=8.5)
fig.text(0.01, 0.99, '(a) DCM in ethanol: commercial UV–vis (black) and home-built spectrometer (colored)', fontsize=8.5, va='top')

ax = fig.add_subplot(gs[2, 0:6])
lim = [0, 2.05]; ax.plot(lim, lim, color='0.75', lw=0.8, ls=':')
xerr = np.array(R['dcm_conc_err_uM']) * R['uv_beer']['eps'] / 1e6
for tag, c, mk, lab in (('lin', C_LIN, 'o', 'sRGB-linearized'), ('meas', C_MEAS, 's', 'measured response'), ('enc', C_ENC, 'D', 'encoded (no correction)')):
    d = R['dcm_' + tag]; p = d['fit_all']
    ax.errorbar(uvA, d['mean'], yerr=d['sd'], xerr=xerr, fmt=mk, ms=5, color=c, mec='white', mew=0.6, elinewidth=1, capsize=0,
                label=f"{lab}: slope {p['slope']:.2f}, R² = {p['R2']:.4f}")
    ax.plot(lim, np.polyval([p['slope'], p['intercept']], lim), color=c, lw=0.8)
RW = json.load(open(A_DIR + 'raw_results.json'))['dcm_raw']
ax.errorbar(uvA, RW['mean'], yerr=RW['sd'], xerr=xerr, fmt='^', ms=5.5, color='#0f766e', mec='white', mew=0.6, elinewidth=1, capsize=0,
            label=f"RAW (DNG), linear sensor data: slope {RW['slope']:.2f}, R² = {RW['R2']:.4f}")
ax.plot(lim, np.polyval([RW['slope'], RW['intercept']], lim), color='#0f766e', lw=0.8)
ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect('equal')
ax.set_xlabel('Absorbance, commercial UV–vis (465–475 nm)'); ax.set_ylabel('Absorbance, home-built (466–473 nm)')
ax.set_title('(b) Linearity at the band maximum', loc='left', fontsize=8.5)
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 0.95), fontsize=7.5, frameon=False)
ax.set_anchor('W')
fig.savefig(FIG + 'Fig5_linearity.png', dpi=300); fig.savefig(FIG + 'Fig5_linearity.pdf')
json.dump(dict(r_dcm=r_vals), open(A_DIR + 'fig5_r.json', 'w'))
print('r per concentration', np.round(r_vals, 3))
