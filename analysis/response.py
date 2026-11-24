"""Camera response (inverse tone curve) of the Galaxy S25 FE in Pro mode, from the shutter series
(Debevec & Malik 1997 method: g(Z) = ln E + ln dt, smoothness-regularised least squares)."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))) + '/'   # repository root
import numpy as np, exifread
from PIL import Image, ImageOps
PH = ROOT + 'photos/raw/ern_photos/20261120_'
FS = ['011823', '011828', '011830', '011832']
def enc(f): return np.asarray(ImageOps.exif_transpose(Image.open(PH + f + '.jpg')).convert('RGB'))[3400:3800, 900:1750]
def ex(f):
    v = exifread.process_file(open(PH + f + '.jpg', 'rb'), details=False)['EXIF ExposureTime'].values[0]; return v.num / v.den
def solve(Z, B, lam=40.0):
    n = 256; N, P = Z.shape
    wt = lambda z: 0.0 if (z >= 250 or z <= 4) else 1.0
    eq = []; rhs = []
    for i in range(N):
        for j in range(P):
            if wt(Z[i, j]) == 0: continue
            r = np.zeros(n + N); r[Z[i, j]] = 1; r[n + i] = -1; eq.append(r); rhs.append(B[j])
    r = np.zeros(n + N); r[128] = 1; eq.append(r); rhs.append(0.0)
    for z in range(1, n - 1):
        r = np.zeros(n + N); r[z - 1] = lam; r[z] = -2 * lam; r[z + 1] = lam; eq.append(r); rhs.append(0.0)
    x = np.linalg.lstsq(np.array(eq), np.array(rhs), rcond=None)[0]; return x[:n]
def fit():
    imgs = [enc(f) for f in FS]; B = np.log([ex(f) for f in FS])
    rng = np.random.default_rng(0)
    m = imgs[2].sum(2) > 60; ys, xs = np.nonzero(m); sel = rng.choice(len(ys), 300, replace=False)
    Z = np.array([[im[ys[s], xs[s], c] for im in imgs] for s in sel for c in range(3)]).astype(int)   # each channel value
    g = solve(Z, B)
    L = np.exp(g - g[250]); L = L / L[250] * ((250/255 + 0.055) / 1.055) ** 2.4; L[:5] = np.linspace(0, L[5], 5)   # scale matched to sRGB at 250
    return L                                   # linear value for each encoded level 0..255, L[255] = 1
if __name__ == '__main__':
    L = fit(); np.save(ROOT + 'analysis/response_L.npy', L)
    s = np.arange(256) / 255; srgb = np.where(s <= 0.04045, s / 12.92, ((s + 0.055) / 1.055) ** 2.4)
    for z in (10, 30, 60, 100, 150, 200, 240): print(z, round(L[z], 4), round(srgb[z], 4))
