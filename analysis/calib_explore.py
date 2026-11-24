"""Explore better wavelength calibrations for the phone spectra (more points / whole-spectrum fit), with honest cross-validation."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))) + '/'   # repository root
import sys, json, numpy as np
from scipy.optimize import minimize
sys.path.insert(0, ROOT + 'uvvis')
from analyse import load as load_uv

z = np.load('arrays.npz'); px = z['px'].astype(float); ref = z['ref_lin']
R = json.load(open('results.json'))
TH = float(sys.argv[1]) if len(sys.argv) > 1 else 0.15
ok = ref > TH * ref.max()
SPEC = {'DCM': (z['dcm_lin'][-1].mean(0), 'a5')}
for i in range(5): SPEC[f'S{i+1}'] = (z['ern_lin'][i].mean(0), f'real{i+1}')
UV = {k: load_uv(f)[:2] for k, (_, f) in SPEC.items()}

def blur(w, a, fwhm):
    g = np.arange(330, 851, 0.5); T = 10 ** (-np.interp(g, w, a, left=0, right=0)); s = max(fwhm, 0.6) / 2.3548
    k = np.exp(-0.5 * (np.arange(-4 * s, 4 * s + 0.5, 0.5) / s) ** 2); k /= k.sum()
    return g, -np.log10(np.clip(np.convolve(T, k, 'same'), 1e-6, None))

xn = (px - 1335) / 565.0                          # normalized pixel (-1..1 over usable range)
def lam(p, model):
    if model == 'lin': return p[0] + p[1] * xn
    if model == 'quad': return p[0] + p[1] * xn + p[2] * xn ** 2
    if model == 'cubic': return p[0] + p[1] * xn + p[2] * xn ** 2 + p[3] * xn ** 3

def sample_cost(nm, fwhm, key, extra=None):
    A, _ = SPEC[key]; w, a = UV[key]; g, ab = blur(w, a, fwhm)
    m = ok & (nm > 380) & (nm < 780)
    u = np.interp(nm[m], g, ab); y = A[m]
    X = np.c_[u, np.ones_like(u)]                 # phone = k*uv + c (scale + offset only; r is invariant to these)
    coef, *_ = np.linalg.lstsq(X, y, rcond=None); res = y - X @ coef
    r = np.corrcoef(u, y)[0, 1]
    return float(np.mean(res ** 2) / np.var(y)), float(r)

def fit(keys, model, p0):
    def f(q):
        p, fw = q[:-1], q[-1]; nm = lam(p, model)
        if np.any(np.diff(nm[ok]) <= 0) or fw < 1 or fw > 60: return 10.0
        return np.mean([sample_cost(nm, fw, k)[0] for k in keys])
    best = None
    for fw0 in (10, 25, 40):
        r = minimize(f, list(p0) + [fw0], method='Nelder-Mead', options=dict(maxiter=6000, xatol=1e-3, fatol=1e-7))
        if best is None or r.fun < best.fun: best = r
    return best.x

keys = list(SPEC)
c = R['cal']; base = [c['offset'] + c['nm_per_px'] * 1335, c['nm_per_px'] * 565]
nm0 = base[0] + base[1] * xn
out = {'threshold': TH, 'usable_px': [float(px[ok][0]), float(px[ok][-1])]}
out['peaks_linear'] = {k: dict(zip(['nrmse', 'r'], sample_cost(nm0, 0.6, k))) for k in keys}
for model, p0 in [('lin', base), ('quad', base + [0]), ('cubic', base + [0, 0])]:
    q = fit(keys, model, p0); nm = lam(q[:-1], model)
    allfit = {k: sample_cost(nm, q[-1], k)[1] for k in keys}
    cv = {}
    for k in keys:                                 # leave-one-spectrum-out: calibrate on the other five, test r on k
        qq = fit([j for j in keys if j != k], model, q[:-1]); cv[k] = sample_cost(lam(qq[:-1], model), qq[-1], k)[1]
    out[model] = dict(params=q[:-1].tolist(), fwhm=float(q[-1]), range_nm=[float(nm[ok][0]), float(nm[ok][-1])],
                      r_fit=allfit, r_cv=cv, nm_per_px_blue=float(np.gradient(nm, px)[ok][0]), nm_per_px_red=float(np.gradient(nm, px)[ok][-1]))
    np.save(f'nm_{model}_{TH}.npy', nm)
json.dump(out, open(f'calib_explore_{TH}.json', 'w'), indent=1)
print(json.dumps({k: (v if not isinstance(v, dict) or 'r_fit' not in v else {kk: v[kk] for kk in ('fwhm', 'range_nm', 'r_fit', 'r_cv', 'nm_per_px_blue', 'nm_per_px_red')}) for k, v in out.items()}, indent=1, default=lambda x: round(x, 3)))
