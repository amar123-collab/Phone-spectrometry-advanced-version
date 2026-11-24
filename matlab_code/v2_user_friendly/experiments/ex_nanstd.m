function s = ex_nanstd(X, dim)
%EX_NANSTD  Sample standard deviation ignoring NaN (MATLAB and Octave).
if nargin < 2, dim = find(size(X) > 1, 1); if isempty(dim), dim = 1; end, end
ok = ~isnan(X);  n = sum(ok, dim);  X(~ok) = 0;
mu = sum(X, dim) ./ n;
d = (X - mu) .* ok;
s = sqrt(sum(d.^2, dim) ./ max(n - 1, 1));
s(n < 2) = NaN;
end
