"""Clean, individually cut-out cuvette photos assembled into Figure 2 (and a standalone version).
Polygons are hand-picked cuvette outlines in the coordinates of the crops made by cuvettes.py
(box from cuvettes.BOX, padded by 40 px, full resolution)."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))) + '/'   # repository root
import numpy as np, cv2, os
from PIL import Image
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from cuvettes import BOX, S, load

POLY = {
 '20261119_231654': [
  [(72, 88), (250, 88), (368, 112), (368, 765), (262, 835), (98, 835), (68, 770)],
  [(85, 85), (345, 85), (347, 770), (285, 822), (110, 822), (80, 760)],
  [(98, 95), (300, 95), (302, 760), (268, 822), (110, 822), (94, 760)],
  [(90, 85), (290, 85), (292, 770), (258, 827), (95, 827), (85, 765)],
  [(80, 70), (262, 70), (264, 800), (247, 842), (103, 842), (78, 800)],
  [(80, 95), (350, 95), (352, 760), (302, 832), (110, 832), (63, 760)],
  [(85, 95), (370, 95), (377, 745), (332, 812), (125, 812), (63, 750)]],
 '20261119_233828': [
  [(72, 130), (130, 62), (418, 60), (448, 830), (400, 925), (160, 928), (140, 880), (95, 400)],
  [(145, 140), (165, 78), (330, 78), (345, 135), (348, 890), (330, 945), (180, 945), (160, 900)],
  [(100, 72), (318, 72), (322, 860), (312, 920), (118, 920), (105, 860)],
  [(85, 80), (310, 75), (360, 120), (362, 840), (332, 937), (140, 937), (65, 860), (70, 150)],
  [(65, 158), (480, 158), (472, 880), (402, 962), (200, 962), (58, 880)]],
}
LABEL = {
 '20261119_231654': ['Ethanol', '7.5 µM', '10.0 µM', '12.5 µM', '24.9 µM', '37.2 µM', '49.5 µM'],
 '20261119_233828': ['Sample 1\nviolet', 'Sample 2\ncyan', 'Sample 3\nmagenta', 'Sample 4\npink', 'Sample 5\nlime'],
}
PAD, FEATHER, TILE_H = 40, 3, 900          # crop padding; edge softening (px); common cuvette height (px)

def tiles(n):
    img = load(n); out = []
    for b, poly in zip(BOX[n], POLY[n]):
        x0, y0, x1, y1 = [int(round(v * S)) for v in b]
        crop = img[y0 - PAD:y1 + PAD, x0 - PAD:x1 + PAD].astype(float)
        h, w = crop.shape[:2]
        m = np.zeros((h, w), np.uint8); cv2.fillPoly(m, [np.array(poly, np.int32)], 255)
        # white balance and exposure from the wall/table around the cuvette (neutral grey -> light grey)
        bg = crop[(cv2.dilate(m, np.ones((25, 25))) == 0)]
        gain = 238.0 / np.median(bg, 0)
        c = np.clip(crop * gain, 0, 255)
        a = cv2.GaussianBlur(m.astype(float) / 255, (0, 0), FEATHER)[..., None]
        comp = c * a + 255 * (1 - a)
        ys, xs = np.nonzero(m); comp = comp[ys.min():ys.max() + 1, xs.min() - 2:xs.max() + 3]
        sc = TILE_H / comp.shape[0]
        out.append(cv2.resize(comp, (int(comp.shape[1] * sc), TILE_H), interpolation=cv2.INTER_AREA).astype(np.uint8))
    return out

def strip(ts, gap=330, margin=80):
    W = sum(t.shape[1] for t in ts) + gap * (len(ts) - 1) + 2 * margin
    canvas = np.full((TILE_H + 2 * margin, W, 3), 255, np.uint8); centres = []; x = margin
    for t in ts:
        canvas[margin:margin + TILE_H, x:x + t.shape[1]] = t; centres.append(x + t.shape[1] / 2); x += t.shape[1] + gap
    return canvas, centres

if __name__ == '__main__':
    FIG = ROOT + 'figs/'
    panels = {n: strip(tiles(n)) for n in POLY}
    for n, (cv, _) in panels.items(): Image.fromarray(cv).save(FIG + f'cuvettes_{n}.png')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8})
    (A, ca), (Bimg, cb) = panels['20261119_231654'], panels['20261119_233828']
    fig = plt.figure(figsize=(7.2, 4.6))
    gs = fig.add_gridspec(2, 1, hspace=0.32)
    for k, (im, cs, labs, title) in enumerate([(A, ca, LABEL['20261119_231654'], '(a) Ethanol and DCM in ethanol'),
                                                (Bimg, cb, LABEL['20261119_233828'], '(b) Mixed-dye event samples')]):
        ax = fig.add_subplot(gs[k]); ax.imshow(im); ax.set_axis_off(); ax.set_anchor('W')
        ax.set_title(title, loc='left', fontsize=8.5)
        for c, l in zip(cs, labs):
            ax.text(c, im.shape[0] + (75 if k == 0 else 20), l, ha='center', va=('baseline' if k == 0 else 'top'), fontsize=8, linespacing=1.15)
    fig.savefig(FIG + 'Fig2_samples.png', dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(FIG + 'Fig2_samples.pdf', bbox_inches='tight', facecolor='white')
    print('done', A.shape, Bimg.shape)
