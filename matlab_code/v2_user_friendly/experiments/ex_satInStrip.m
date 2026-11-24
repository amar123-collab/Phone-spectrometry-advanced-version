function s = ex_satInStrip(P, strip, idx, halfWin)
%EX_SATINSTRIP  True for each strip position idx that has saturated pixels
%   within +/- halfWin positions (default 3), across the strip width.
if nargin < 4, halfWin = 3; end
r = round(strip.roi(1:2));  c = round(strip.roi(3:4));
S = P.sat(r(1):r(2), c(1):c(2));
if strip.vertical, prof = any(S, 2)'; else, prof = any(S, 1); end
if strip.reverse, prof = fliplr(prof); end
s = false(size(idx));
for j = 1:numel(idx)
    i = max(1, idx(j) - halfWin):min(numel(prof), idx(j) + halfWin);
    s(j) = any(prof(i));
end
end
