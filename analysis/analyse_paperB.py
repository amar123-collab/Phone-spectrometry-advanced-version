"""Paper B analysis on the real data of 20 Nov 2026 (Galaxy S25 FE, prism spectrometer).
Mirrors the MATLAB pipeline: sRGB linearisation, dark subtraction, strip profile following the
spectrum's tilt, absorbance against ethanol references. Writes results.json and arrays.npz."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))) + '/'   # repository root
import json, sys
import numpy as np
from datetime import datetime
from PIL import Image, ImageOps
import exifread

PH = ROOT + 'photos/raw/ern_photos/20261120_'
UV = ROOT + 'uvvis/'
sys.path.insert(0, UV)
from analyse import load as load_uv          # baseline 700-800 nm, lamp-artefact patching

COLS = (700, 1950)
HALF = 40                                     # strip half-height (rows) around the spectrum centre line

_cache = {}
RESP = np.load(ROOT + 'analysis/response_L.npy')   # measured inverse tone curve (response.py)
TAG = lambda lin: 'meas' if lin == 'meas' else ('lin' if lin else 'enc')
def img(f, lin=True):
    key = (f, lin)
    if key not in _cache:
        if len(_cache) >= 2: _cache.pop(next(iter(_cache)))
        a = np.asarray(ImageOps.exif_transpose(Image.open(PH + f + '.jpg')).convert('RGB')).astype(np.float32) / 255
        if lin == 'meas':
            a = RESP[np.round(a * 255).astype(np.int16)].astype(np.float32)
        elif lin:
            a = np.where(a <= 0.04045, a / 12.92, ((a + 0.055) / 1.055) ** 2.4)
        _cache[key] = a
    return _cache[key]

def exif(f):
    t = exifread.process_file(open(PH + f + '.jpg', 'rb'), details=False)
    e = t.get('EXIF ExposureTime').values[0]
    return dict(time=str(t.get('EXIF DateTimeOriginal')), exposure=float(e.num) / float(e.den),
                iso=int(str(t.get('EXIF ISOSpeedRatings'))), wb=str(t.get('EXIF WhiteBalance')))

# ---------------- sets ----------------
DARK = ['002855', '002857']
ETH_A = ['002937', '002950', '003005', '003012']                       # not moved (start)
ETH_B = ['003045', '003051', '003059', '003104', '003109', '003116', '003122', '003126',
         '003131', '003134', '003138']                                 # taken out and put back
ETH_C = ['003654', '003657', '003659', '003701']                       # end of series
NOCUV = ['003707']
DCM = [['003202', '003208', '003211'], ['003227', '003230', '003232'], ['003246', '003248', '003251'],
       ['003310', '003312', '003314'], ['003326', '003329', '003332'], ['003349', '003352', '003355']]
DCM_UL = np.array([3, 4, 5, 10, 15, 20.]); DCM_C = 5000 * DCM_UL / (2000 + DCM_UL)
DCM_UV = ['b1', 'b2', 'a2', 'a3', 'a4', 'a5']
ERN = [['003411', '003414', '003417'], ['003429', '003432', '003435'], ['003450', '003453', '003455'],
       ['003509', '003512', '003514'], ['003526', '003529', '003530']]
ERN_UV = ['real1', 'real2', 'real3', 'real4', 'real5']
SHUTTER = ['011823', '011828', '011830', '011832', '011834', '011839']
EXPTEST = ['001124', '001159', '002033', '002315']
AUTO = ['011859']

# ---------------- strip following the tilt (from ethanol references) ----------------
refimg = sum(img(f).sum(2) for f in ETH_C) / len(ETH_C)
xs = np.arange(*COLS); yc = []
for x in xs:
    col = refimg[3300:3650, x]; w = np.clip(col - 0.2 * col.max(), 0, None)
    yc.append(3300 + (np.arange(col.size) * w).sum() / w.sum() if w.sum() > 0 else np.nan)
yc = np.array(yc); good = ~np.isnan(yc) & (refimg[3300:3650, COLS[0]:COLS[1]].max(0) > 0.2 * refimg.max())
tilt = np.polyfit(xs[good], yc[good], 1)
centre = np.polyval(tilt, xs).round().astype(int)
TILT_DEG = float(np.degrees(np.arctan(tilt[0])))

_pc = {}
def profile(f, lin=True):
    if (f, lin) in _pc: return _pc[(f, lin)].copy()
    a = img(f, lin).sum(2)
    _pc[(f, lin)] = np.array([a[c - HALF:c + HALF, x].mean() for x, c in zip(xs, centre)])
    return _pc[(f, lin)].copy()

def dark(lin):   return np.mean([profile(f, lin) for f in DARK], 0)

R = {}; ARR = {}
R['tilt_deg'] = TILT_DEG
# ---- alignment along the dispersion direction (cuvette placement shifts the image by a few px)
BLUE_W = (xs >= 775) & (xs <= 835); RED_W = (xs >= 1835) & (xs <= 1895)
def est_shift(p, ref, wins):
    best = (np.inf, 0.0)
    for sft in np.arange(-15, 15.01, 0.25):
        q = np.interp(xs + sft, xs, p); cost = 0
        for w in wins:
            d = np.log(np.clip(q[w], 1e-7, None)) - np.log(np.clip(ref[w], 1e-7, None)); cost += ((d - d.mean()) ** 2).sum()
        best = min(best, (cost, sft))
    return best[1]
D0 = dark(True)
ALL_ETH = ETH_A + ETH_B + ETH_C + NOCUV
ref0 = np.mean([profile(f) - D0 for f in ETH_B + ETH_C], 0)
SHIFT = {}
for it in range(2):                                   # align references to their mean, twice
    for f in ETH_B + ETH_C: SHIFT[f] = est_shift(profile(f) - D0, ref0, [BLUE_W, RED_W])
    ref0 = np.mean([np.interp(xs + SHIFT[f], xs, profile(f) - D0) for f in ETH_B + ETH_C], 0)
for f in ETH_A + NOCUV: SHIFT[f] = est_shift(profile(f) - D0, ref0, [BLUE_W, RED_W])
WINS = {}
for i, s_ in enumerate(DCM):  WINS.update({f: [RED_W] if i >= 3 else [BLUE_W, RED_W] for f in s_})
for i, s_ in enumerate(ERN):  WINS.update({f: {0: [BLUE_W, RED_W], 1: [BLUE_W], 2: [BLUE_W, RED_W], 3: [BLUE_W, RED_W], 4: [RED_W]}[i] for f in s_})
for f, w in WINS.items(): SHIFT[f] = est_shift(profile(f) - D0, ref0, w)
R['shift_px'] = SHIFT
ALIGN = False   # edge-based alignment was unreliable for coloured samples (shifts hit the ±15 px search limit); not used
def aprof(f, lin, D):
    p = profile(f, lin) - D
    return np.interp(xs + SHIFT.get(f, 0.0), xs, p) if ALIGN else p
for lin in (True, False, 'meas'):
    tag = TAG(lin)
    D = dark(lin)
    P = {f: aprof(f, lin, D) for f in ALL_ETH}
    ref = np.mean([P[f] for f in ETH_B + ETH_C], 0)
    ok = ref > 0.15 * ref.max()
    ARR[f'ref_{tag}'] = ref; ARR[f'ok_{tag}'] = ok
    absorb = lambda f: -np.log10(np.clip(aprof(f, lin, D), 1e-7, None) / ref)
    S = np.array([[absorb(f) for f in s] for s in DCM]); ARR[f'dcm_{tag}'] = S
    E = np.array([[absorb(f) for f in s] for s in ERN]); ARR[f'ern_{tag}'] = E
    bl = {}
    for f in ALL_ETH:
        others = [P[g] for g in ETH_B + ETH_C if g != f]
        bl[f] = -np.log10(np.clip(P[f], 1e-7, None) / np.mean(others, 0))
    ARR[f'blank_{tag}'] = np.array([bl[f] for f in ALL_ETH])
    R[f'ref_total_{tag}'] = {f: float(P[f][ok].sum()) for f in ALL_ETH}
ok = ARR['ok_lin']; px = xs
R['usable_px'] = [int(px[ok].min()), int(px[ok].max())]

# ---------------- DCM analytical window (same for both modes) ----------------
top = ARR['dcm_lin'][-1].mean(0).copy(); top[~ok] = np.nan
sm = np.convolve(np.nan_to_num(top), np.ones(21) / 21, 'same'); sm[~ok] = np.nan
pk = int(np.nanargmax(sm)); WIN = slice(pk - 15, pk + 16)
R['dcm_peak_px'] = int(px[pk])
uvA = []
for f in DCM_UV:
    w, a, _ = load_uv(f); uvA.append(a[(w >= 465) & (w <= 475)].mean())
uvA = np.array(uvA); R['dcm_uvA'] = uvA.tolist(); R['dcm_conc_uM'] = DCM_C.tolist()
R['dcm_conc_err_uM'] = (5000 * 2000 / (2000 + DCM_UL) ** 2 * 0.3).tolist()
pu = np.polyfit(DCM_C, uvA, 1); r = uvA - np.polyval(pu, DCM_C)
R['uv_beer'] = dict(eps=float(pu[0] * 1e6), intercept=float(pu[1]), R2=float(1 - (r ** 2).sum() / ((uvA - uvA.mean()) ** 2).sum()))
for tag in ('lin', 'enc', 'meas'):
    M = np.nanmean(ARR[f'dcm_{tag}'][:, :, WIN], 2)                # conc x rep
    m, sd = M.mean(1), M.std(1, ddof=1)
    out = dict(mean=m.tolist(), sd=sd.tolist(), reps=M.tolist())
    for lab, sel in (('all', slice(None)), ('le1', uvA <= 1.0)):
        p = np.polyfit(uvA[sel], m[sel], 1); rr = m[sel] - np.polyval(p, uvA[sel])
        out[f'fit_{lab}'] = dict(slope=float(p[0]), intercept=float(p[1]),
                                 R2=float(1 - (rr ** 2).sum() / ((m[sel] - m[sel].mean()) ** 2).sum()))
    out['ratio_to_uv'] = (m / uvA).tolist()
    R[f'dcm_{tag}'] = out
R['meas_slope'] = R['dcm_meas']['fit_all']['slope']
R['enc_over_lin_slope'] = R['dcm_enc']['fit_all']['slope'] / R['dcm_lin']['fit_all']['slope']

# ---------------- blanks / repeatability / LOD ----------------
names = ALL_ETH
for tag in ('lin', 'enc', 'meas'):
    B = np.nanmean(ARR[f'blank_{tag}'][:, WIN], 1)
    d = dict(zip(names, B.tolist()))
    nb = np.array([d[f] for f in ETH_B]); na = np.array([d[f] for f in ETH_A])
    R[f'blank_{tag}'] = dict(values=d, sd_replaced=float(nb.std(ddof=1)), mean_replaced=float(nb.mean()),
                             sd_notmoved=float(na.std(ddof=1)), n_replaced=len(nb), n_notmoved=len(na))
slope = R['dcm_lin']['fit_all']['slope']
R['lod_A'] = 3 * R['blank_lin']['sd_replaced'] / slope
R['lod_uM'] = R['lod_A'] / (R['uv_beer']['eps'] / 1e6)
R['nocuvette_A'] = R['blank_lin']['values'][NOCUV[0]]
# lamp drift: total reference signal vs time
R['ref_times'] = {f: exif(f)['time'][11:] for f in ETH_A + ETH_B + ETH_C + NOCUV}

# ---------------- wavelength calibration from sample bands (reference-based) ----------------
def peak_px(A):
    A = A.copy(); A[~ok] = np.nan
    s = np.convolve(np.nan_to_num(A), np.ones(21) / 21, 'same'); s[~ok] = np.nan
    i = int(np.nanargmax(s))
    if 2 <= i < len(s) - 2 and np.all(np.isfinite(s[i - 2:i + 3])):     # parabolic refinement on smoothed data
        y0, y1, y2 = s[i - 1], s[i], s[i + 1]; dd = y0 - 2 * y1 + y2
        if dd < 0: return px[i] + 0.5 * (y0 - y2) / dd
    return float(px[i])
uv_peak = {}
for f in DCM_UV + ERN_UV:
    w, a, _ = load_uv(f); m = (w >= 380) & (w <= 750); uv_peak[f] = float(w[m][np.argmax(a[m])])
bands = {'DCM': (peak_px(ARR['dcm_lin'][-1].mean(0)), 470.0)}
for i in (0, 2, 3):
    bands[f'ERN{i+1}'] = (peak_px(ARR['ern_lin'][i].mean(0)), uv_peak[ERN_UV[i]])
for i in (1, 4):
    bands[f'ERN{i+1}_edge'] = (peak_px(ARR['ern_lin'][i].mean(0)), uv_peak[ERN_UV[i]])
R['bands'] = bands
names_c = ['DCM', 'ERN1', 'ERN3', 'ERN4']
P = np.array([bands[n][0] for n in names_c]); L = np.array([bands[n][1] for n in names_c])
cal = np.polyfit(P, L, 1)
R['cal'] = dict(nm_per_px=float(cal[0]), offset=float(cal[1]), resid=dict(zip(names_c, (np.polyval(cal, P) - L).round(1).tolist())))
loo = {}
for i, n in enumerate(names_c):
    m = np.arange(len(P)) != i; c = np.polyfit(P[m], L[m], 1); loo[n] = float(np.polyval(c, P[i]) - L[i])
R['cal']['loo'] = loo; R['cal']['loo_rms'] = float(np.sqrt(np.mean(np.array(list(loo.values())) ** 2)))
NM = np.polyval(cal, px); ARR['nm'] = NM
R['usable_nm'] = [float(NM[ok].min()), float(NM[ok].max())]
R['dcm_window_nm'] = [float(NM[WIN][0]), float(NM[WIN][-1])]

# ---------------- effective resolution: blur UV-Vis to match phone spectra ----------------
def blur(w, a, fwhm):
    g = np.arange(380, 801, 0.5); T = 10 ** (-np.interp(g, w, a)); s = fwhm / 2.3548
    k = np.exp(-0.5 * (np.arange(-4 * s, 4 * s + 0.5, 0.5) / s) ** 2); k /= k.sum()
    Tb = np.convolve(T, k, 'same'); return g, -np.log10(np.clip(Tb, 1e-6, None))
def match(phoneA, uvf, fwhm):
    w, a, _ = load_uv(uvf); g, ab = blur(w, a, fwhm)
    u = np.interp(NM, g, ab); m = ok & np.isfinite(phoneA)
    x, y = u[m], phoneA[m]; k = (x @ y) / (x @ x)
    return float(np.sqrt(np.mean((y - k * x) ** 2))), float(np.corrcoef(x, y)[0, 1]), k
fw = np.arange(2, 41, 1.0); res = {}
for lab, A, f in [('DCM49', ARR['dcm_lin'][-1].mean(0), 'a5'), ('DCM25', ARR['dcm_lin'][3].mean(0), 'a3')] + \
        [(f'ERN{i+1}', ARR['ern_lin'][i].mean(0), ERN_UV[i]) for i in range(5)]:
    rms = [match(A, f, x)[0] for x in fw]; best = float(fw[int(np.argmin(rms))])
    r_sharp = match(A, f, 0.5)[1]; r_best = match(A, f, best)[1]
    res[lab] = dict(best_fwhm_nm=best, r_unblurred=r_sharp, r_blurred=r_best)
R['spectral_match'] = res
fw_ern = [res[k]['best_fwhm_nm'] for k in ('ERN1', 'ERN3', 'ERN4')]
R['resolution_estimate_nm'] = dict(median=float(np.median(fw_ern)), values=fw_ern)

# ---------------- shutter series (separate geometry): camera response ----------------
sh = []
base = img(SHUTTER[2]).sum(2)[3400:3800, 900:1750]          # 1/90 s photo defines the spectrum pixels
MASK = base > 0.2 * base.max()
for f in SHUTTER:
    e = exif(f)['exposure']; out = {'file': f, 'exposure_s': e}
    enc = img(f, False); sat = (enc >= 254 / 255).any(2)[3400:3800, 900:1750]
    out['sat_pct'] = float((sat & MASK).sum() / MASK.sum() * 100)
    for lin in (True, False, 'meas'):
        I = img(f, lin).sum(2)[3400:3800, 900:1750]
        bg = np.median(I[~MASK])
        out['signal_' + TAG(lin)] = float((I[MASK] - bg).mean())
    sh.append(out)
e = np.array([s_['exposure_s'] for s_ in sh])
for tag in ('lin', 'enc', 'meas'):
    y = np.array([s_['signal_' + tag] for s_ in sh])
    sel = np.array([s_['sat_pct'] < 0.5 for s_ in sh])
    p = np.polyfit(np.log10(e[sel]), np.log10(y[sel]), 1)
    R[f'shutter_loglog_slope_{tag}'] = float(p[0]); R[f'shutter_n_{tag}'] = int(sel.sum())
R['shutter'] = sh
# exposure tests and auto-mode: clipped fraction inside the spectrum
clip = {}
for f in EXPTEST + AUTO + ['002315', '003104']:
    enc = img(f, False); I = img(f).sum(2); m = I > 0.2 * I.max()
    clip[f] = dict(exif(f), clipped_pct_of_spectrum=float(((enc >= 254 / 255).any(2) & m).sum() / m.sum() * 100))
R['clipping'] = clip

json.dump(R, open(ROOT + 'analysis/results.json', 'w'), indent=1, default=float)
np.savez(ROOT + 'analysis/arrays.npz', px=px, centre=centre, **ARR)
print('shifts', {k: round(v,2) for k,v in SHIFT.items()})
print(json.dumps({k: R[k] for k in ['tilt_deg', 'usable_px', 'usable_nm', 'dcm_peak_px', 'dcm_window_nm', 'uv_beer',
                  'enc_over_lin_slope', 'lod_A', 'lod_uM', 'nocuvette_A', 'cal', 'resolution_estimate_nm',
                  'shutter_loglog_slope_lin', 'shutter_loglog_slope_enc', 'shutter_loglog_slope_meas']}, indent=1, default=float))
for t in ('lin', 'enc', 'meas'):
    d = R[f'dcm_{t}']; print(t, np.round(d['mean'], 3), np.round(d['sd'], 3), d['fit_all'], d['fit_le1'])
    print('  blank', {k: round(v, 4) for k, v in R[f'blank_{t}'].items() if k != 'values'})
print('spectral match', json.dumps(R['spectral_match'], indent=0))
print('shutter', [(round(s['exposure_s'], 4), round(s['sat_pct'], 2)) for s in sh])
print('clipping', {k: (round(v['exposure'], 4), v['iso'], round(v['clipped_pct_of_spectrum'], 1)) for k, v in clip.items()})
