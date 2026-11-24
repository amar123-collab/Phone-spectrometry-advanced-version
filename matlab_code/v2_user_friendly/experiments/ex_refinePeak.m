function pk = ex_refinePeak(pr, pk)
%EX_REFINEPEAK  Sub-pixel peak positions by a parabola through 3 points.
p = double(pr(:)');
for j = 1:numel(pk)
    i = round(pk(j));
    if i > 1 && i < numel(p)
        d = p(i - 1) - 2 * p(i) + p(i + 1);
        if d < 0, pk(j) = i + 0.5 * (p(i - 1) - p(i + 1)) / d; end
    end
end
end
