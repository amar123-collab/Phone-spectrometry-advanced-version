function [A, P] = ex_phoneA(file, linearise, strip, anchor, align, smoothPx, Idark, Iref, minRefFrac)
%EX_PHONEA  Absorbance spectrum of one photo along a fixed strip.
P = psm_loadPhoto(file, 'full', linearise);
[sh, ok] = psm_shift(P, anchor, align);
if ~ok, warning('%s: no blue end found, not aligned.', P.name); end
I = psm_profile(P, strip, sh, smoothPx) - Idark;
A = psm_absorbance(I, Iref, minRefFrac);
end
