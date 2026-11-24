function S = ex_blurSpectrum(S, fwhmNm)
%EX_BLURSPECTRUM  Convolve a spectrum [wl, A] with a Gaussian of given FWHM (nm).
%   Used to bring a UV-Vis spectrum to the resolution of the phone spectrometer.
%   Works on transmittance, as the instrument blurs light, not absorbance.
[wl, o] = sort(S(:, 1));  A = S(o, 2);
g  = (min(wl):0.25:max(wl))';
T  = 10 .^ (-interp1(wl, A, g));
sg = fwhmNm / 2.3548;
k  = exp(-0.5 * ((-ceil(4 * sg / 0.25):ceil(4 * sg / 0.25))' * 0.25 / sg) .^ 2);  k = k / sum(k);
Tb = conv(T, k, 'same');
e  = numel(k);  Tb(1:e) = T(1:e);  Tb(end - e + 1:end) = T(end - e + 1:end);   % keep the ends
S  = [g, -log10(max(Tb, 1e-6))];
end
