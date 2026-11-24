"""Cut out each cuvette from the two sample photos (GrabCut + convex hull), place on a clean white
background, and assemble labelled panels for Figure 2."""
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))) + '/'   # repository root
import cv2, numpy as np, os, sys
from PIL import Image, ImageOps

SRC = ROOT + 'photos/raw/ern_photos/'
DBG = sys.argv[1] if len(sys.argv) > 1 else '/tmp'
S = 4000 / 1400                      # preview -> full-resolution scale

# generous boxes in preview coordinates (x0, y0, x1, y1)
BOX = {
 '20261119_231654': [(55, 520, 190, 812), (248, 525, 368, 812), (445, 530, 550, 815), (645, 532, 750, 818),
                     (858, 532, 958, 818), (1023, 538, 1145, 820), (1220, 540, 1345, 822)],
 '20261119_233828': [(65, 505, 240, 820), (380, 495, 497, 822), (675, 488, 785, 810), (928, 483, 1057, 813),
                     (1175, 480, 1343, 815)],
}

def load(n):
    return np.array(ImageOps.exif_transpose(Image.open(SRC + n + '.jpg')).convert('RGB'))

def cutout(img, box, pad=40):
    x0, y0, x1, y1 = [int(round(v * S)) for v in box]
    x0 -= pad; y0 -= pad; x1 += pad; y1 += pad
    crop = img[y0:y1, x0:x1].copy()
    h, w = crop.shape[:2]
    mask = np.full((h, w), cv2.GC_PR_BGD, np.uint8)
    mask[:pad, :] = mask[-pad // 2:, :] = cv2.GC_BGD; mask[:, :pad] = mask[:, -pad:] = cv2.GC_BGD
    iw, ih = w - 2 * pad, h - 2 * pad
    mask[pad + int(0.06 * ih):pad + int(0.80 * ih), pad + int(0.12 * iw):pad + int(0.88 * iw)] = cv2.GC_PR_FGD
    mask[pad + int(0.10 * ih):pad + int(0.75 * ih), pad + int(0.35 * iw):pad + int(0.65 * iw)] = cv2.GC_FGD
    rect = None; bgd = np.zeros((1, 65), np.float64); fgd = np.zeros((1, 65), np.float64)
    bgr = cv2.cvtColor(crop, cv2.COLOR_RGB2BGR)
    cv2.grabCut(bgr, mask, rect, bgd, fgd, 8, cv2.GC_INIT_WITH_MASK)
    fg = np.where((mask == 1) | (mask == 3), 255, 0).astype(np.uint8)
    fg = cv2.morphologyEx(fg, cv2.MORPH_OPEN, np.ones((9, 9), np.uint8))
    n, lab, stats, _ = cv2.connectedComponentsWithStats(fg)
    if n > 1:
        k = 1 + np.argmax(stats[1:, cv2.CC_STAT_AREA]); fg = np.where(lab == k, 255, 0).astype(np.uint8)
    cnts, _ = cv2.findContours(fg, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    hull = cv2.convexHull(np.vstack(cnts))
    m = np.zeros((h, w), np.uint8); cv2.fillConvexPoly(m, hull, 255)
    return crop, m

if __name__ == '__main__':
    out = []
    for n, boxes in BOX.items():
        img = load(n)
        for i, b in enumerate(boxes):
            crop, m = cutout(img, b)
            dbg = crop.copy(); cnt, _ = cv2.findContours(m, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            cv2.drawContours(dbg, cnt, -1, (255, 0, 0), 6)
            np.save(f'{DBG}/{n}_{i}_mask.npy', m); Image.fromarray(crop).save(f'{DBG}/{n}_{i}_crop.png')
            out.append(Image.fromarray(dbg).resize((dbg.shape[1] // 4, dbg.shape[0] // 4)))
    W = sum(o.width for o in out); H = max(o.height for o in out)
    sheet = Image.new('RGB', (W, H), 'white'); x = 0
    for o in out: sheet.paste(o, (x, 0)); x += o.width
    sheet.save(f'{DBG}/debug_masks.png'); print('ok', sheet.size)
