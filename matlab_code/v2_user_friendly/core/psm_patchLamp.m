function A = psm_patchLamp(wl, A, win)
%PSM_PATCHLAMP  Smooth over lamp/filter-change artefacts in a laboratory UV-Vis spectrum.
%   A = psm_patchLamp(wl, A)        default windows of the instrument used for the paper
%   A = psm_patchLamp(wl, A, win)   win: one row [lo hi] (nm) per artefact; [] = no patching
%   Inside each window the absorbance is replaced by a quadratic fitted to the 6 nm on either
%   side (same as uvvis/analyse.py). Change the windows for your own instrument; a single
%   column of centres c is also accepted and treated as [c-2 c+2].
if nargin < 3, win = [484 488; 630 641; 653 659]; end
if isempty(win), return; end
if size(win, 2) == 1, win = [win(:) - 2, win(:) + 2]; end
wl = wl(:);  sz = size(A);  A = A(:);
for i = 1:size(win, 1)
    lo = win(i, 1);  hi = win(i, 2);
    bad = wl >= lo & wl <= hi;
    nb = (wl >= lo - 6 & wl < lo) | (wl > hi & wl <= hi + 6);
    if nnz(bad) > 0 && nnz(nb) >= 4
        c = (lo + hi) / 2;  pp = polyfit(wl(nb) - c, A(nb), 2);  A(bad) = polyval(pp, wl(bad) - c);
    end
end
A = reshape(A, sz);
end
