function F = ex_linfit(x, y)
%EX_LINFIT  Least-squares line y = m x + c with standard errors and R^2.
x = x(:);  y = y(:);  ok = ~isnan(x) & ~isnan(y);  x = x(ok);  y = y(ok);
n = numel(x);  X = [x, ones(n, 1)];
b = X \ y;  r = y - X * b;
F.m = b(1);  F.c = b(2);  F.n = n;
F.R2 = 1 - sum(r.^2) / sum((y - mean(y)).^2);
if n > 2
    s2 = sum(r.^2) / (n - 2);  C = s2 * inv(X' * X);
    F.sm = sqrt(C(1, 1));  F.sc = sqrt(C(2, 2));
else
    F.sm = NaN;  F.sc = NaN;
end
end
