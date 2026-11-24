"""Baseline correction from the repeated ethanol blanks.
The 19 ethanol photographs (same cuvette, removed and replaced) vary in three smooth shapes (principal components).
Each sample spectrum is fitted with these shapes (plus an offset) ONLY where the sample is transparent; the fitted
baseline is then subtracted everywhere. The laboratory spectrum is used only to decide where the sample is transparent."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))) + '/'   # repository root
import json, sys, numpy as np
sys.path.insert(0, ROOT + 'uvvis')
from analyse import load as load_uv
A_DIR = ROOT + 'analysis/'
Z = np.load(A_DIR + 'arrays.npz'); R = json.load(open(A_DIR + 'results.json'))
px = Z['px']; ok = Z['ok_lin']; NM = R['cal']['offset'] + R['cal']['nm_per_px'] * px
win = (NM >= R['dcm_window_nm'][0]) & (NM <= R['dcm_window_nm'][1])
B = Z['blank_lin'][:19]                                   # ETH_A + ETH_B + ETH_C, each against the mean reference
Bd = B[:, ok] - B[:, ok].mean(0)
U, s, Vt = np.linalg.svd(Bd, full_matrices=False)
frac = (s ** 2 / (s ** 2).sum())
PCS = np.zeros((3, len(px))); PCS[:, ok] = Vt[:3]
DCM_UV = ['b1', 'b2', 'a2', 'a3', 'a4', 'a5']

def transparent(uvfile, A_shape=None):
    w, u = load_uv(uvfile)[:2]; uu = np.interp(NM, w, u)
    return ok & (uu < 0.05 * uu[ok].max()), uu

def basis(kind):
    cols = [np.ones_like(NM)]
    if kind == 'line': cols.append(NM - 550)
    if kind.startswith('pc'): cols += list(PCS[:int(kind[2])])
    return np.array(cols).T

def correct(A, tr, kind):
    if kind == 'none': return A.copy(), 0
    X = basis(kind); c, *_ = np.linalg.lstsq(X[tr], A[tr], rcond=None); b = X @ c
    return A - b, b

KINDS = ['none', 'offset', 'line', 'pc1', 'pc2', 'pc3']
out = {'pc_variance': frac[:5].tolist(), 'dcm': {}, 'ern': {}}
uvA = np.array(R['dcm_uvA'])
for kind in KINDS:
    rs, win_vals, sds = [], [], []
    for i, f in enumerate(DCM_UV):
        tr, uu = transparent(f)
        reps = np.array([correct(a, tr, 'offset' if kind == 'none' else kind)[0] if kind != 'none' else a for a in Z['dcm_lin'][i]])
        m = reps.mean(0); rs.append(float(np.corrcoef(uu[ok], m[ok])[0, 1]))
        v = reps[:, win].mean(1); win_vals.append(float(v.mean())); sds.append(float(v.std(ddof=1)))
    p = np.polyfit(uvA, win_vals, 1); rr = np.array(win_vals) - np.polyval(p, uvA)
    out['dcm'][kind] = dict(r=rs, window=win_vals, sd=sds, slope=float(p[0]), intercept=float(p[1]),
                            R2=float(1 - (rr ** 2).sum() / ((np.array(win_vals) - np.mean(win_vals)) ** 2).sum()))
    rs = []
    for i in range(5):
        tr, uu = transparent(f'real{i+1}')
        m = Z['ern_lin'][i].mean(0)
        if tr.sum() < 20 or kind == 'none': c = m
        else: c = correct(m, tr, kind)[0]
        rs.append(float(np.corrcoef(uu[ok], c[ok])[0, 1]))
    out['ern'][kind] = dict(r=rs, n_transparent=[int(transparent(f'real{i+1}')[0].sum()) for i in range(5)])
# blank-based check: does the PC baseline also remove the scatter of the blanks themselves? (leave-one-out PCs)
loo_sd = {}
for kind in ('none', 'pc1', 'pc2', 'pc3'):
    vals = []
    for j in range(4, 15):                                # ETH_B photos (removed and replaced)
        others = np.delete(np.arange(19), j); Bo = B[others][:, ok] - B[others][:, ok].mean(0)
        _, _, V = np.linalg.svd(Bo, full_matrices=False)
        P = np.zeros((3, len(px))); P[:, ok] = V[:3]
        a = B[j]; tr = ok & (NM > 560)                     # pretend only the DCM-transparent region is known
        if kind == 'none': c = a
        else:
            X = np.array([np.ones_like(NM)] + list(P[:int(kind[2])])).T; cc, *_ = np.linalg.lstsq(X[tr], a[tr], rcond=None); c = a - X @ cc
        vals.append(c[win].mean())
    loo_sd[kind] = float(np.std(vals, ddof=1))
out['blank_window_sd_after_correction'] = loo_sd
np.save(A_DIR + 'blank_pcs.npy', PCS)
json.dump(out, open(A_DIR + 'blank_baseline.json', 'w'), indent=1)
print('PC variance', np.round(frac[:4], 3))
for kind in KINDS:
    d = out['dcm'][kind]; e = out['ern'][kind]
    print(f'{kind:7s} DCM r', np.round(d['r'], 3), f"slope {d['slope']:.3f} int {d['intercept']:+.3f} R2 {d['R2']:.4f}", '| ERN r', np.round(e['r'], 2))
print('blank SD in DCM window after correction (fit >560 nm only):', {k: round(v, 3) for k, v in loo_sd.items()})

# ---------------- figure ----------------
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8.5, 'axes.linewidth': 0.8, 'axes.spines.top': False, 'axes.spines.right': False})
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.7), gridspec_kw=dict(width_ratios=[1.2, 1.2, 1]))
cols = ['#1d4ed8'] * 4 + ['#c2410c'] * 11 + ['#6b7280'] * 4
for j in range(19):
    ax[0].plot(NM[ok], B[j, ok], color=cols[j], lw=0.7, alpha=0.8)
ax[0].set_xlabel('Wavelength (nm)'); ax[0].set_ylabel('Absorbance'); ax[0].set_ylim(-0.4, 0.4)
ax[0].set_title('(a) 19 ethanol blanks', loc='left', fontsize=8.5)
for j, c in enumerate(['#111827', '#c2410c', '#1d4ed8']):
    ax[1].plot(NM[ok], PCS[j, ok] * np.sqrt(ok.sum()), color=c, lw=1.2, label=f'Shape {j+1} ({100*frac[j]:.0f}%)')
ax[1].axhline(0, color='0.8', lw=0.6); ax[1].set_xlabel('Wavelength (nm)'); ax[1].set_ylabel('Relative amplitude')
ax[1].set_title('(b) Shapes in the blanks', loc='left', fontsize=8.5)
ax[1].legend(frameon=False, fontsize=7, loc='upper center', bbox_to_anchor=(0.5, -0.25), ncol=1)
labs = ['none', 'pc1', 'pc2', 'pc3']; names = ['None', '1 shape', '2 shapes', '3 shapes']
ax[2].bar(range(4), [loo_sd[k] for k in labs], color=['0.7', '#fdba74', '#f97316', '#c2410c'], width=0.6)
ax[2].set_xticks(range(4)); ax[2].set_xticklabels(names, fontsize=7.5, rotation=20)
ax[2].set_ylabel('SD at 466–473 nm'); ax[2].set_title('(c) Placement scatter', loc='left', fontsize=8.5)
fig.tight_layout(); fig.savefig(ROOT + 'figs/S_blank_baseline.png', dpi=300); fig.savefig(ROOT + 'figs/S_blank_baseline.pdf')
