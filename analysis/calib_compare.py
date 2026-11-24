"""Compare ways to improve the agreement of phone and laboratory spectra of the mixed-dye solutions (2026 set)."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))) + '/'   # repository root
import numpy as np, json, sys
sys.path.insert(0, ROOT + 'uvvis')
from analyse import load as load_uv
z = np.load('arrays.npz'); px = z['px'].astype(float); ref = z['ref_lin']; R = json.load(open('results.json'))
ok = ref > 0.15 * ref.max()
SPEC = {'DCM': (z['dcm_lin'][-1], 'a5')}
for i in range(5): SPEC[f'S{i+1}'] = (z['ern_lin'][i], f'real{i+1}')
UV = {k: load_uv(f)[:2] for k, (_, f) in SPEC.items()}
def blur(w, a, fwhm):
    g = np.arange(330, 851, 0.5); T = 10 ** (-np.interp(g, w, a, left=0, right=0)); s = max(fwhm, 0.6) / 2.3548
    k = np.exp(-0.5 * (np.arange(-4 * s, 4 * s + 0.5, 0.5) / s) ** 2); k /= k.sum()
    return g, -np.log10(np.clip(np.convolve(T, k, 'same'), 1e-6, None))
sm = lambda y: np.convolve(y, np.ones(21) / 21, 'same')

# ---- calibration points: band maxima + half-maximum edges (phone vs UV blurred to 25 nm) ----
def half_edges(x, y, i, h, okm):
    out = []
    for step in (-1, 1):
        j = i
        while 0 < j < len(y) - 1 and okm[j] and y[j] > h: j += step
        out.append(x[j] if okm[j] and y[j] <= h else None)
    return out
pts = []
for k, (A, _) in SPEC.items():
    a = sm(A.mean(0)); a = np.where(ok, a, -9); i = int(np.argmax(a)); h = a[i] / 2
    pb, pr = half_edges(px, a, i, h, ok)
    w, u = UV[k]; g, ub = blur(w, u, 25); mm = (g > 400) & (g < 700); g, ub = g[mm], ub[mm]
    j = int(np.argmax(ub)); ubb, urr = half_edges(g, ub, j, ub[j] / 2, np.ones_like(g, bool))
    edge = k in ('S2', 'S5')                           # maxima outside the usable range
    if not edge: pts.append((k, 'max', px[i], g[j]))
    if pb is not None and ubb is not None: pts.append((k, 'blue half-max', pb, ubb))
    if pr is not None and urr is not None: pts.append((k, 'red half-max', pr, urr))
P = np.array([p[2] for p in pts]); Lm = np.array([p[3] for p in pts])
CAL = {'4 maxima (paper)': np.array([R['cal']['nm_per_px'], R['cal']['offset']])}
CAL[f'{len(pts)} points, linear'] = np.polyfit(P, Lm, 1)
CAL[f'{len(pts)} points, quadratic'] = np.polyfit(P, Lm, 2)
loo = {}
for name, deg in ((f'{len(pts)} points, linear', 1), (f'{len(pts)} points, quadratic', 2)):
    e = []
    for k in SPEC:
        mk = np.array([p[0] != k for p in pts]); c = np.polyfit(P[mk], Lm[mk], deg)
        e += list(np.polyval(c, P[~mk]) - Lm[~mk])
    loo[name] = float(np.sqrt(np.mean(np.square(e))))
loo['4 maxima (paper)'] = R['cal']['loo_rms']

def baseline(k, A, nm):
    """Straight line through the points where the laboratory spectrum is transparent (<5% of max) on BOTH sides of the band;
    offset only if transparent points exist on one side; no correction if none."""
    w, a = UV[k]; u = np.interp(nm, w, a); tr = ok & (u < 0.05 * u[ok].max())
    if tr.sum() < 20: return np.zeros_like(A), 'none'
    pk = nm[ok][np.argmax(u[ok])]; both = (nm[tr] < pk - 10).any() and (nm[tr] > pk + 10).any()
    c = np.polyfit(nm[tr], A[tr], 1 if both else 0); return np.polyval(c, nm), 'line' if both else 'offset'

def rr(k, nm, fw, base):
    A = SPEC[k][0].mean(0); w, a = UV[k]; g, ab = blur(w, a, fw); u = np.interp(nm, g, ab)
    b, kind = baseline(k, A, nm) if base else (0, '-')
    if kind != 'line': b = 0
    return float(np.corrcoef(u[ok], (A - b)[ok])[0, 1]), kind

res = {'points': [(p[0], p[1], float(p[2]), float(p[3])) for p in pts], 'loo_rms_nm': loo, 'cases': {}}
for cname, c in CAL.items():
    nm = np.polyval(c, px)
    for fw in (0.6, 25):
        for base in (False, True):
            key = f'{cname} | bandpass {"none" if fw < 1 else "25 nm"} | baseline {"yes" if base else "no"}'
            res['cases'][key] = {k: rr(k, nm, fw, base) for k in SPEC}
    res.setdefault('range', {})[cname] = [float(nm[ok][0]), float(nm[ok][-1])]
    res.setdefault('maxima_nm', {})[cname] = {k: float(np.polyval(c, R['bands'][kk][0])) for k, kk in
                                             [('DCM', 'DCM'), ('S1', 'ERN1'), ('S3', 'ERN3'), ('S4', 'ERN4')]}
json.dump(res, open('calib_compare.json', 'w'), indent=1)
print('points:'); [print('  ', p) for p in res['points']]
print('LOO rms', {k: round(v, 1) for k, v in loo.items()}); print('range', {k: np.round(v, 0).tolist() for k, v in res['range'].items()})
print('maxima', {k: {kk: round(vv) for kk, vv in v.items()} for k, v in res['maxima_nm'].items()})
for key, v in res['cases'].items(): print(f'{key:70s}', ' '.join(f'{k}:{x[0]:.2f}' for k, x in v.items()))
