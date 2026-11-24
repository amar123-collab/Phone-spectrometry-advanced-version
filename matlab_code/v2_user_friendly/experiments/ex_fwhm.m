function [w, lo, hi] = ex_fwhm(pr, pk, halfWin)
%EX_FWHM  Full width at half maximum of the peak at index pk (pixels).
%   The local baseline is the minimum within +/- halfWin pixels (default 40).
%   lo, hi are the interpolated half-maximum crossing positions.
if nargin < 3, halfWin = 40; end
p = double(pr(:)');  n = numel(p);
i1 = max(1, pk - halfWin);  i2 = min(n, pk + halfWin);
base = min(p(i1:i2));
h = base + (p(pk) - base) / 2;
j = pk;  while j > i1 && p(j) > h, j = j - 1; end
lo = j + (h - p(j)) / (p(j + 1) - p(j));
j = pk;  while j < i2 && p(j) > h, j = j + 1; end
hi = j - (h - p(j)) / (p(j - 1) - p(j));
w = hi - lo;
if p(i1) > h || p(i2) > h, w = NaN; end          % peak not resolved inside window
end
