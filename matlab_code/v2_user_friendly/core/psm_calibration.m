function cal = psm_calibration(type, pix, nm, order)
%PSM_CALIBRATION  Build a pixel -> wavelength calibration.
%   cal = psm_calibration('none')
%   cal = psm_calibration('laser', pix, nm, order)
%         known lines (laser spots, lamp lines) at positions pix along the
%         strip; order 1 = grating (straight line), 2 = prism (curve)
%   cal = psm_calibration('prism', pix, nm)
%         prism dispersion model  pix = a + b / nm^2  fitted to independent
%         lines (laser spots, lamp lines); 3+ lines give a residual check
%   cal = psm_calibration('model', pix, nm)
%         straight line nm = a + b * pix fitted to the absorption maxima of
%         samples with known lambda_max (as in the paper). A straight line is
%         used because such peaks usually span only part of the range, and a
%         curved prism model extrapolates poorly beyond them. NOT independent:
%         use only when no laser/lamp calibration exists, and say so when reporting.
%
%   Use x = psm_wavelength(cal, 1:n) to convert positions to nm.
cal = struct('type', type, 'p', [], 'resid', [], 'info', '');
switch lower(type)
    case 'none'
        cal.info = 'No calibration: x-axis in pixels (blue -> red).';
    case 'laser'
        if nargin < 4, order = 1; end
        order = min(order, numel(pix) - 1);
        cal.p = polyfit(pix(:), nm(:), order);
        cal.resid = nm(:)' - polyval(cal.p, pix(:)');
        if numel(pix) > order + 1
            cal.info = sprintf('Laser calibration (order %d), residuals %s nm.', order, mat2str(round(10 * cal.resid) / 10));
        else
            cal.info = sprintf('Laser calibration (order %d), exact fit - add a line to check accuracy.', order);
        end
    case 'prism'
        X = [ones(numel(pix), 1), 1 ./ nm(:).^2];
        cal.p = X \ pix(:);                       % [a; b]
        cal.resid = psm_wavelength(cal, pix(:)') - nm(:)';
        cal.info = sprintf('Prism-model calibration (pix = a + b/nm^2) on %d lines, residuals %s nm.', ...
                           numel(pix), mat2str(round(10 * cal.resid) / 10));
    case 'model'
        cal.p = polyfit(pix(:), nm(:), 1);
        cal.resid = polyval(cal.p, pix(:)') - nm(:)';
        cal.range = [min(nm), max(nm)];
        cal.info = sprintf(['Estimated axis (straight line through %d known peaks, %.0f-%.0f nm), residuals %s nm. ', ...
                            'Not an independent calibration; outside %.0f-%.0f nm the axis is extrapolated.'], ...
                           numel(pix), min(nm), max(nm), mat2str(round(10 * cal.resid) / 10), min(nm), max(nm));
    otherwise
        error('psm:cal', 'Unknown calibration type "%s".', type);
end
end
