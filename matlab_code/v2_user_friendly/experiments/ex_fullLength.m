function strip = ex_fullLength(strip, sz)
%EX_FULLLENGTH  Extend a strip along the rainbow to (almost) the full photo.
%   A margin of 3 % at each end leaves room for alignment shifts.
if strip.vertical, d = 1; else, d = 2; end
m = max(20, round(0.03 * sz(d)));
k = 2 * d - 1;                                % roi index of the start along the rainbow
strip.roi(k:k + 1) = [1 + m, sz(d) - m];
end
