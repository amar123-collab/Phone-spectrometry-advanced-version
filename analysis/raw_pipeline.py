"""Same analysis from the RAW (DNG) files: linear 12-bit sensor data, no gamma, no tone curve.
RAW columns are mapped onto the JPG pixel axis with the band maxima of DCM and samples 1, 3, 4, so the same wavelength calibration applies."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))) + '/'   # repository root
import json, sys, numpy as np, rawpy
sys.path.insert(0, ROOT + 'uvvis')
from analyse import load as load_uv
PH = ROOT + 'photos/raw/ern_photos/20261120_'
A_DIR = ROOT + 'analysis/'
z = np.load(A_DIR + 'arrays.npz'); R = json.load(open(A_DIR + 'results.json'))
cal = R['cal']; nm_of = lambda p: cal['offset'] + cal['nm_per_px'] * p

DARK = ['002855', '002857']
ETH_A = ['002937', '002950', '003005', '003012']
ETH_B = ['003045', '003051', '003059', '003104', '003109', '003116', '003122', '003126', '003131', '003134', '003138']
ETH_C = ['003654', '003657', '003659', '003701']
DCM = [['003202', '003208', '003211'], ['003227', '003230', '003232'], ['003246', '003248', '003251'],
       ['003310', '003312', '003314'], ['003326', '003329', '003332'], ['003349', '003352', '003355']]
ERN = [['003411', '003414', '003417'], ['003429', '003432', '003435'], ['003450', '003453', '003455'],
       ['003509', '003512', '003514'], ['003526', '003529', '003530']]

def dng(f):
    with rawpy.imread(PH + f + '.dng') as r:
        im = r.postprocess(gamma=(1, 1), no_auto_bright=True, output_bps=16, user_wb=[1, 1, 1, 1],
                           output_color=rawpy.ColorSpace.raw, use_camera_wb=False)
    return im.astype(float).sum(2)

ref_img = np.mean([dng(f) for f in ETH_C], 0)
xs = np.arange(600, 2000); yc = []
for x in xs:
    col = ref_img[3350:3750, x]; w = np.clip(col - 0.2 * col.max(), 0, None)
    yc.append(3350 + (np.arange(col.size) * w).sum() / w.sum() if w.sum() > 0 else np.nan)
yc = np.array(yc); good = ~np.isnan(yc) & (ref_img[3350:3750, 600:2000].max(0) > 0.2 * ref_img.max())
cen = np.polyval(np.polyfit(xs[good], yc[good], 1), xs).round().astype(int)
P = {}
def prof(f):
    if f not in P: P[f] = np.array([dng(f)[c - 41:c + 41, x].mean() for x, c in zip(xs, cen)]) if False else None
    return P[f]
# faster: decode each file once and take the strip
def strip(f):
    a = dng(f); return np.array([a[c - 41:c + 41, x].mean() for x, c in zip(xs, cen)])
ALL = DARK + ETH_A + ETH_B + ETH_C + sum(DCM, []) + sum(ERN, [])
for f in ALL: P[f] = strip(f)
D = np.mean([P[f] for f in DARK], 0)
I0 = np.mean([P[f] for f in ETH_B + ETH_C], 0) - D
Afun = lambda f: -np.log10(np.clip(P[f] - D, 1e-3, None) / np.clip(I0, 1e-3, None))
ok = I0 > 0.15 * I0.max()
sm = lambda y: np.convolve(y, np.ones(21) / 21, 'same')
def peak(A, lo, hi):
    s = sm(A); m = ok & (xs > lo) & (xs < hi); return float(xs[m][np.argmax(s[m])])
dcm = np.array([[Afun(f) for f in g] for g in DCM]); ern = np.array([[Afun(f) for f in g] for g in ERN])
# map RAW columns -> JPG px with the four band maxima used for the JPG calibration
raw_pk = [peak(dcm[-1].mean(0), 850, 1150), peak(ern[0].mean(0), 1350, 1650), peak(ern[2].mean(0), 1300, 1600), peak(ern[3].mean(0), 1200, 1500)]
jpg_pk = [R['bands'][k][0] for k in ('DCM', 'ERN1', 'ERN3', 'ERN4')]
mp = np.polyfit(raw_pk, jpg_pk, 1)
uv_pk = [R['bands'][k][1] for k in ('DCM', 'ERN1', 'ERN3', 'ERN4')]
cr = np.polyfit(raw_pk, uv_pk, 1); NM = np.polyval(cr, xs)        # RAW-specific calibration on the same four bands
loo = [float(np.polyval(np.polyfit(np.delete(raw_pk, i), np.delete(uv_pk, i), 1), raw_pk[i]) - uv_pk[i]) for i in range(4)]
win = (NM >= R['dcm_window_nm'][0]) & (NM <= R['dcm_window_nm'][1])
uvA = np.array(R['dcm_uvA']); m = dcm[:, :, win].mean(2)
mean, sd = m.mean(1), m.std(1, ddof=1); p = np.polyfit(uvA, mean, 1); rr = mean - np.polyval(p, uvA)
out = dict(raw_peaks=raw_pk, cal_raw=cr.tolist(), cal_raw_loo=loo, cal_raw_loo_rms=float(np.sqrt(np.mean(np.square(loo)))), map_raw_to_jpg=mp.tolist(), map_resid_px=(np.polyval(mp, raw_pk) - np.array(jpg_pk)).tolist(),
           usable_nm={str(t): [float(NM[I0 > t * I0.max()].min()), float(NM[I0 > t * I0.max()].max())] for t in (0.15, 0.05, 0.02)},
           dcm_raw=dict(mean=mean.tolist(), sd=sd.tolist(), slope=float(p[0]), intercept=float(p[1]),
                        R2=float(1 - (rr ** 2).sum() / ((mean - mean.mean()) ** 2).sum()), ratio=(mean / uvA).tolist()))
B = np.array([Afun(f)[win].mean() for f in ETH_B]); Aa = np.array([Afun(f)[win].mean() for f in ETH_A])
out['blank_raw'] = dict(sd_replaced=float(B.std(ddof=1)), sd_notmoved=float(Aa.std(ddof=1)))
out['lod_raw_A'] = 3 * out['blank_raw']['sd_replaced'] / p[0]
r_ = {}
for i in range(5):
    w, u = load_uv(f'real{i+1}')[:2]; uu = np.interp(NM, w, u); A = ern[i].mean(0)
    r_[f'S{i+1}'] = float(np.corrcoef(uu[ok], A[ok])[0, 1])
w, u = load_uv('a5')[:2]; r_['DCM'] = float(np.corrcoef(np.interp(NM, w, u)[ok], dcm[-1].mean(0)[ok])[0, 1])
out['r_raw'] = r_
np.savez(A_DIR + 'raw_arrays.npz', nm=NM, I0=I0, ok=ok, dcm=dcm, ern=ern, xs=xs)
json.dump(out, open(A_DIR + 'raw_results.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
