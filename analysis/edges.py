import numpy as np
def edge_pos(p, px, side, frac=0.5):
    """sub-pixel x where profile crosses frac*max on the red (right) or blue (left) flank"""
    i0 = int(np.argmax(np.convolve(p, np.ones(15)/15, 'same'))); h = frac * p[i0]
    if side == 'red':
        j = i0
        while j < len(p) - 1 and p[j] > h: j += 1
        return px[j-1] + (p[j-1] - h) / (p[j-1] - p[j]) * (px[j] - px[j-1])
    j = i0
    while j > 0 and p[j] > h: j -= 1
    return px[j] + (h - p[j]) / (p[j+1] - p[j]) * (px[j+1] - px[j])
