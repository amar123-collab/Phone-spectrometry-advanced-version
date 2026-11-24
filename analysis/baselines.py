"""Baseline correction rules used for Figures 5a and 8b.
Transparent region: laboratory absorbance below 5% of its maximum (and below 0.02 for DCM).
- Transparent over less than 50 nm: no correction (the baseline cannot be estimated).
- Transparent on both sides of the band: straight line through the transparent points.
- Transparent on one side only: non-negative combination of the 19 measured ethanol blanks plus an offset,
  fitted to the transparent points (the blanks supply realistic baseline shapes to extrapolate with)."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))) + '/'   # repository root
import json, sys, numpy as np
from scipy.optimize import lsq_linear
sys.path.insert(0, ROOT + 'uvvis')
from analyse import load as load_uv
A_DIR = ROOT + 'analysis/'
Z = np.load(A_DIR + 'arrays.npz'); R = json.load(open(A_DIR + 'results.json'))
px = Z['px']; ok = Z['ok_lin']; NM = R['cal']['offset'] + R['cal']['nm_per_px'] * px
BK = Z['blank_lin'][:19]
STEP = float(np.median(np.diff(NM)))

def baseline(A, uvfile, dcm=False):
    w, u = load_uv(uvfile)[:2]; uu = np.interp(NM, w, u)
    thr = min(0.02, 0.05 * uu[ok].max()) if dcm else 0.05 * uu[ok].max()
    tr = ok & (uu < thr)
    if tr.sum() * STEP < 50: return np.zeros_like(A), 'none'
    pk = NM[ok][np.argmax(uu[ok])]
    if (NM[tr] < pk - 10).any() and (NM[tr] > pk + 10).any():
        c = np.polyfit(NM[tr], A[tr], 1); return np.polyval(c, NM), 'line'
    X = np.hstack([BK.T, np.ones((len(A), 1))])
    res = lsq_linear(X[tr], A[tr], bounds=(np.r_[np.zeros(19), -np.inf], np.inf))
    return X @ res.x, 'ethanol blanks'
