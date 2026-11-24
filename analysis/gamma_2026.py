"""Encoded vs linearized absorbance for the 2026 set (Pro mode, fixed exposure; mean of three photographs per sample).
Slopes through the origin for linearized absorbance above 0.2, over the usable range (reference above 15% of its maximum)."""
import json, numpy as np
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))) + '/'
Z = np.load(ROOT + 'analysis/arrays.npz'); ok = Z['ok_lin']
out = {'samples': {}, 'dcm_49uM': None}; pts = {}; X, Y = [], []
for i in range(5):
    l = Z['ern_lin'][i].mean(0); e = Z['ern_enc'][i].mean(0); m = ok & np.isfinite(l) & np.isfinite(e); f = m & (l > 0.2)
    out['samples'][str(i + 1)] = float((e[f] * l[f]).sum() / (l[f] ** 2).sum()); pts[f'x{i+1}'] = l[m]; pts[f'y{i+1}'] = e[m]
    X.append(l[f]); Y.append(e[f])
l = Z['dcm_lin'][-1].mean(0); e = Z['dcm_enc'][-1].mean(0); f = ok & (l > 0.2)
out['dcm_49uM'] = float((e[f] * l[f]).sum() / (l[f] ** 2).sum()); pts['xD'] = l[ok]; pts['yD'] = e[ok]
X, Y = np.concatenate(X), np.concatenate(Y); out['slope_all'] = float((X * Y).sum() / (X ** 2).sum())
np.savez(ROOT + 'analysis/gamma_2026.npz', **pts); json.dump(out, open(ROOT + 'analysis/gamma_2026.json', 'w'), indent=1)
print(json.dumps(out))
