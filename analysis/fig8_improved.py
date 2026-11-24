"""Figure 8 (two rows: as measured / after cuvette-baseline correction and bandpass matching) and SI figure on calibration options."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))) + '/'   # repository root
import json, sys, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, ROOT + 'uvvis'); sys.path.insert(0, ROOT + 'analysis')
from analyse import load as load_uv
import calib_compare as cc                      # runs the comparison and exposes helpers/results

R = cc.R; z = cc.z; px = cc.px; ok = cc.ok; res = cc.res
FIG = ROOT + 'figs/'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8.5, 'axes.linewidth': 0.8,
                     'axes.spines.top': False, 'axes.spines.right': False})
ERN_COL = ['#6d3fb0', '#1f78b4', '#c2185b', '#e0457b', '#9aa800']
NM = R['cal']['offset'] + R['cal']['nm_per_px'] * px
def save(fig, n): fig.savefig(FIG + n + '.png', dpi=300); fig.savefig(FIG + n + '.pdf')
mask = lambda y: np.where(ok, y, np.nan)

# ---------- main Figure 8 ----------
import baselines as BL
# sample 5 below the JPG range: RAW files with a stray-light correction (stray level from the 49.5 uM DCM solution)
W = np.load(ROOT + 'analysis/raw_arrays.npz'); nmR = W['nm']; I0R = W['I0']
winR = (nmR >= 466) & (nmR <= 473); wD, uD = load_uv('a5')[:2]
Tm = 10 ** -W['dcm'][-1].mean(0)[winR].mean(); Tt = 10 ** -uD[(wD >= 465) & (wD <= 475)].mean()
STRAY = (Tm - Tt) / (1 - Tm) * I0R[winR].mean()
I5 = 10 ** (-W['ern'][4]) * I0R
A5R = -np.log10(np.clip(I5 - STRAY, 1e-9, None) / np.clip(I0R - STRAY, 1e-9, None))
good = I0R > 10 * STRAY; low = ~good & (nmR < 450) & (nmR > 412)
# red end of the axis: linear extrapolation beyond the last calibration band (sample 1, pixel 1478) is not constrained by any point.
# A quadratic bend beyond that pixel is fitted to the shape of sample 2 (the only band there); bands at or below 587 nm are unchanged.
PX1 = R['bands']['ERN1'][0]
def nm_bent(c): return NM + c * np.clip(px - PX1, 0, None) ** 2
def r_corr(i, nmx):
    w, u = cc.UV[f'S{i+1}']; g, ub = cc.blur(w, u, 25); BL.NM = nmx
    A = z['ern_lin'][i].mean(0); b, _ = BL.baseline(A, f'real{i+1}'); BL.NM = NM
    sel = ok & (nmx >= 450) if i == 4 else ok
    return np.corrcoef(np.interp(nmx[sel], g, ub), (A - b)[sel])[0, 1]
grid = -np.arange(0, 1.3e-4, 1e-5)
mean_r = [np.mean([r_corr(i, nm_bent(c)) for i in range(5)]) for c in grid]
CB = float(grid[int(np.argmax(mean_r))]); NMB = nm_bent(CB); END_B = float(NMB[ok].max())
print('bend', CB, 'red end', END_B, 'mean r', np.round(mean_r, 4))

def panel(a, row, i, norm, kinds):
    k = f'S{i+1}'; A = z['ern_lin'][i]; w, u = cc.UV[k]
    if row == 0:
        g, ub = w, u; reps = A; nmx = NM
    else:
        g, ub = cc.blur(w, u, 25); nmx = NMB; BL.NM = NMB
        out = [BL.baseline(x, f'real{i+1}') for x in A]; reps = np.array([x - b for x, (b, _) in zip(A, out)]); kinds.append(out[0][1])
        BL.NM = NM
    m = (g >= 400) & (g <= 690)
    okp = ok & (nmx >= 450) if (row == 1 and i == 4) else ok
    ph = np.where(okp, reps.mean(0), np.nan); sd = np.where(okp, reps.std(0, ddof=1), np.nan)
    r = np.corrcoef(np.interp(nmx[okp], g, ub), reps.mean(0)[okp])[0, 1]
    extra = []
    if row == 1 and i == 4:
        mr = good & (nmR < 460); ml = low
        sm9 = lambda y: np.convolve(y, np.ones(25) / 25, 'same')       # 25-column (about 5 nm) running mean for display
        extra = [(nmR[mr][12:-12], sm9(A5R.mean(0))[mr][12:-12], '--'), (nmR[ml][12:-12], sm9(W['ern'][4].mean(0))[ml][12:-12], ':')]
    if norm == 'self':                      # lab: own maximum; phone: maximum within its range matched to the lab maximum over that range
        sl = ub[m].max(); fin = np.isfinite(ph)
        lab_in = np.interp(nmx[fin], g, ub).max(); sp = np.nanmax(ph) * sl / lab_in
        spr = (extra[0][1].max() * sl / np.interp(extra[0][0], g, ub).max()) if extra else 1.0   # same rule for the RAW stream (stray-corrected part)
    else:
        sl = sp = spr = 1.0
    a.plot(g[m], ub[m] / sl, color='black', lw=1); a.plot(nmx, ph / sp, color=ERN_COL[i], lw=1.4)
    a.fill_between(nmx, (ph - sd) / sp, (ph + sd) / sp, color=ERN_COL[i], alpha=0.25, lw=0)
    for xx, yy, ls in extra: a.plot(xx, yy / spr, color='#0f766e', lw=1.4, ls=ls)   # RAW: own scale factor, same rule (dotted part uses the dashed part's factor)
    lo, hi = nmx[ok].min(), nmx[ok].max()
    a.axvspan(400, lo, color='0.92', zorder=0); a.axvspan(hi, 690, color='0.92', zorder=0)
    a.set_xlim(400, 690); a.set_xticks([450, 550, 650])
    if norm == 'self': a.set_ylim(-0.35, 1.15)
    return r

def fig8(norm, name):
    """SI version: both rows (as measured / corrected)."""
    fig, axs = plt.subplots(2, 5, figsize=(7.2, 4.6), sharex=True, sharey=(norm == 'self'))
    r_rows = {}; kinds = []
    for row in (0, 1):
        for i in range(5):
            a = axs[row, i]; r = panel(a, row, i, norm, kinds); r_rows.setdefault(row, []).append(r)
            lab = f'r = {r:.2f}' + ('' if row == 0 else f'\n({kinds[-1]})')
            a.set_title((f'Sample {i+1}\n' if row == 0 else '') + lab, fontsize=7)
            if i == 0: a.set_ylabel('Normalized absorbance' if norm == 'self' else 'Absorbance')
    fig.supxlabel('Wavelength (nm)', fontsize=8.5, y=0.01)
    fig.subplots_adjust(left=0.1, right=0.99, top=0.86, bottom=0.11, hspace=0.72, wspace=0.14 if norm == 'self' else 0.32)
    fig.text(0.01, 0.975, '(a) As measured', fontsize=8.5, va='top')
    fig.text(0.01, 0.50, '(b) Corrected: baseline removed, red end of axis adjusted, laboratory at 25 nm bandpass, sample 5 < 450 nm from RAW files (teal)',
             fontsize=7.5, va='top')
    save(fig, name); return r_rows, kinds

def fig_main(name):
    """Main-text figure: the five cuvettes above their corrected, normalized spectra."""
    from cuvette_panel import tiles
    T = tiles('20261119_233828')
    fig = plt.figure(figsize=(7.2, 4.0))
    gs = fig.add_gridspec(2, 5, height_ratios=[0.8, 1], hspace=0.42, wspace=0.14, left=0.1, right=0.99, top=0.87, bottom=0.12)
    names = ['violet', 'cyan', 'magenta', 'pink', 'lime']; kinds = []; sh = None
    for i in range(5):
        c = fig.add_subplot(gs[0, i]); c.imshow(T[i]); c.set_axis_off()
        c.set_title(f'Sample {i+1} ({names[i]})', fontsize=7.5)
        a = fig.add_subplot(gs[1, i], sharey=sh); sh = sh or a
        r = panel(a, 1, i, 'self', kinds)
        a.set_title(f'r = {r:.2f}', fontsize=7)
        if i == 0: a.set_ylabel('Normalized absorbance')
        else: plt.setp(a.get_yticklabels(), visible=False)
    fig.text(0.55, 0.01, 'Wavelength (nm)', ha='center', fontsize=8.5)
    fig.text(0.01, 0.985, '(a) The 2026 set of mixed-dye solutions', fontsize=8.5, va='top')
    yb = gs[1, 0].get_position(fig).y1
    fig.text(0.01, yb + 0.075, '(b) Phone (colored; sample 5 below 450 nm from RAW files, teal) and laboratory spectra (black)', fontsize=8.5, va='bottom')
    save(fig, name)

r_rows, kinds = fig8('self', 'S_event_samples_normalized')
fig8('none', 'S_event_samples_absolute')
fig_main('Fig5_event_samples')
json.dump({'bend_c': CB, 'bent_red_end_nm': END_B, 'stray_fraction_of_I0max': float(STRAY / I0R.max()), 'raw_good_from_nm': float(nmR[good].min()), 'raw_low_from_nm': float(nmR[low].min()) if low.any() else None,
           'S5_raw_at': {str(l): float(A5R.mean(0)[np.argmin(abs(nmR - l))]) for l in (420, 430, 440, 450)}}, open(ROOT + 'analysis/s5_blue.json', 'w'), indent=1)

# ---------- SI figure: calibration options ----------
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.3), gridspec_kw=dict(width_ratios=[1, 1.5]))
P = np.array([p[2] for p in res['points']]); Lm = np.array([p[3] for p in res['points']])
kinds = {'max': 'o', 'blue half-max': '<', 'red half-max': '>'}
cols = {'DCM': '#c2410c', 'S1': ERN_COL[0], 'S2': ERN_COL[1], 'S3': ERN_COL[2], 'S4': ERN_COL[3], 'S5': ERN_COL[4]}
for (k, kind, p, l) in res['points']:
    ax[0].plot(p, l, kinds[kind], color=cols[k], mec='white', mew=0.5, ms=6)
xx = np.linspace(760, 1910, 10)
names = list(cc.CAL)
for name, ls in zip(names, ['-', '--', ':']):
    ax[0].plot(xx, np.polyval(cc.CAL[name], xx), color='0.3', lw=1, ls=ls, label=f'{name} (LOO {res["loo_rms_nm"][name]:.1f} nm)')
ax[0].set_xlabel('Pixel'); ax[0].set_ylabel('Wavelength, laboratory (nm)')
ax[0].set_title('(a) Calibration points', loc='left', fontsize=8.5)
ax[0].legend(frameon=False, fontsize=6.5, loc='upper center', bbox_to_anchor=(0.5, -0.2))
cases = [(f'{names[0]} | bandpass none | baseline no', 'As in Figure 8a'),
         (f'{names[1]} | bandpass none | baseline no', f'{names[1]}'),
         (f'{names[0]} | bandpass 25 nm | baseline yes', 'Baseline + 25 nm bandpass (Figure 8b)'),
         (f'{names[1]} | bandpass 25 nm | baseline yes', f'{names[1]} + baseline + bandpass')]
keys = ['DCM', 'S1', 'S2', 'S3', 'S4', 'S5']; wdt = 0.2; shades = ['0.75', '0.55', '#c2410c', '#7c2d12']
for j, (key, lab) in enumerate(cases):
    v = [res['cases'][key][k][0] for k in keys]
    ax[1].bar(np.arange(6) + (j - 1.5) * wdt, v, wdt, color=shades[j], label=lab)
ax[1].set_xticks(range(6)); ax[1].set_xticklabels(['DCM'] + [f'Sample {i}' for i in range(1, 6)], fontsize=7)
ax[1].set_ylim(0.5, 1.0); ax[1].set_ylabel('Correlation, r')
ax[1].set_title('(b) Agreement of spectral shape', loc='left', fontsize=8.5)
ax[1].legend(frameon=False, fontsize=6.5, loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)
fig.tight_layout(); save(fig, 'S_calibration_options')
json.dump({'fig8_r_measured': r_rows[0], 'fig8_r_corrected': r_rows[1], 'methods': kinds}, open(ROOT + 'analysis/fig8_r.json', 'w'))
print('Fig 8 r measured', np.round(r_rows[0], 2), 'corrected', np.round(r_rows[1], 2))
