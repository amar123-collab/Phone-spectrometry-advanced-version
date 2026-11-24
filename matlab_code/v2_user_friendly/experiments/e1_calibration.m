%% e1_calibration.m  -  Experiment E1: wavelength calibration and resolution
% Edit the SETTINGS block and press Run (MATLAB R2020b+ or GNU Octave).
%
% Needs, all taken WITHOUT moving the phone or the screen:
%   * a reference photo (white light, solvent cuvette) - defines the strip
%   * one or more photos of calibration lines with the lamp blocked:
%     laser spots (one photo per laser, or all in one photo) or a compact
%     fluorescent lamp (CFL) at the entrance of the spectrometer
%   * optional: photos of ONE line at different iris openings (resolution)
%
% Output (in the current folder):
%   E1_calibration.mat   calibration + strip, used by e3_dilution.m
%   E1_results.txt       numbers for the paper (residuals, dispersion, FWHM)
%   E1_lines.csv         line positions and fit residuals
%   E1_calibration.png/.pdf   check figure (draft of Figure 4)

clear; close all;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'core'));  addpath(here);

%% ======================= SETTINGS =======================
folder    = fullfile(here, 'E1_photos');
refFile   = 'ref_white.jpg';                  % white light through solvent cuvette
linearise = true;

% Calibration photos and the wavelengths (nm) of the lines in each photo.
% Laser pointers:  lineFiles = {'laser405.jpg','laser532.jpg','laser650.jpg'};
%                  lineNm    = {405, 532, 650};
% One CFL photo:   lineFiles = {'cfl.jpg'};  lineNm = {[435.8 487 546.1 611]};
%   (Hg 404.7 / 435.8 / 546.1 nm; Tb 487, 543 nm; Eu about 611 nm). List the
%   lines from blue to red; the n brightest peaks are paired with them in
%   order, so CHECK the figure. If a wrong peak is picked, read the right
%   pixel positions off panel (a) and type them into linePx.
lineFiles = {'laser405.jpg', 'laser532.jpg', 'laser650.jpg'};
lineNm    = {405, 532, 650};
linePx    = [];                               % [] = automatic; else pixel positions, blue -> red

% Resolution: photos of one line at different iris openings ({} = skip).
% For a laser, put a diffuser (tracing paper, ground glass) just before the
% iris so that the laser fills the iris like the lamp does; otherwise the
% spot width shows the laser beam, not the instrument.
resFiles  = {};                               % e.g. {'res532_iris8mm.jpg', 'res532_iris3mm.jpg'}
resLabels = {};                               % e.g. {'iris 8 mm (working)', 'iris 3 mm'}
resNm     = 532;                              % wavelength of that line

minPeakFrac = 0.05;                           % weakest accepted line, fraction of strongest
%% ========================================================

ref = psm_loadPhoto(fullfile(folder, refFile), 'full', linearise);
% Both strips run over the WHOLE photo along the rainbow, so that lines at
% the dim ends (e.g. 405 nm) are not cut off; the dim ends are masked later.
strip  = ex_fullLength(psm_findStrip(ref), size(ref.I));        % for absorbance (e3_dilution)
stripW = ex_fullLength(psm_findStrip(ref, 1.0), size(ref.I));   % full rainbow height: catches spots slightly off axis
Iref   = psm_profile(ref, stripW, [0 0], 1);
n      = numel(Iref);  px = 1:n;
% lamp-blocked photos have no blue end to align on: the phone must not move
fprintf('Reference %s: strip %d px long.\n', ref.name, n);

%% ---- calibration lines
pkAll = [];  nmAll = [];  prLines = zeros(numel(lineFiles), n);  satWarn = {};
for k = 1:numel(lineFiles)
    L  = psm_loadPhoto(fullfile(folder, lineFiles{k}), 'full', linearise);
    pr = psm_profile(L, stripW, [0 0], 1);
    prLines(k, :) = pr;
    nm = sort(lineNm{k}(:)');
    pk = psm_findPeaks(pr, numel(nm), 10, minPeakFrac);
    if numel(pk) < numel(nm)
        error('%s: found %d line(s), expected %d. Lower minPeakFrac or use linePx.', L.name, numel(pk), numel(nm));
    end
    pk = ex_refinePeak(pr, pk);
    s = ex_satInStrip(L, stripW, round(pk));
    if any(s), satWarn{end+1} = L.name; end %#ok<SAGROW>
    pkAll = [pkAll, pk];  nmAll = [nmAll, nm]; %#ok<AGROW>
end
[nmAll, o] = sort(nmAll);  pkAll = pkAll(o);
if ~isempty(linePx)
    if numel(linePx) ~= numel(nmAll), error('linePx must have %d values.', numel(nmAll)); end
    pkAll = linePx(:)';
end
if ~issorted(pkAll) && ~issorted(fliplr(pkAll))
    warning('Line positions are not in wavelength order - check pairing in panel (a).');
end

cal = psm_calibration('prism', pkAll, nmAll);
fprintf('%s\n', cal.info);
calPoly = [];
if numel(pkAll) >= 4
    calPoly = psm_calibration('laser', pkAll, nmAll, 2);   % independent check of the model
end
lam  = psm_wavelength(cal, px);
disp_nm = abs(gradient(lam));                 % nm per pixel along the strip

%% ---- resolution
res = struct('label', {}, 'fwhmPx', {}, 'fwhmNm', {}, 'pk', {});
for k = 1:numel(resFiles)
    R  = psm_loadPhoto(fullfile(folder, resFiles{k}), 'full', linearise);
    pr = psm_profile(R, stripW, [0 0], 1);
    [~, i0] = min(abs(lam - resNm));
    pk = psm_findPeaks(pr, 3, 10, 0.2);
    [~, j] = min(abs(pk - i0));  pk = pk(j);
    w  = ex_fwhm(pr, pk);
    res(k).label = resLabels{k};  res(k).pk = pk;  res(k).fwhmPx = w;
    res(k).fwhmNm = w * disp_nm(pk);  res(k).profile = pr;
    if any(ex_satInStrip(R, stripW, pk)), satWarn{end+1} = R.name; end %#ok<SAGROW>
end
% widths of the calibration lines themselves (laser = instrument function)
lineFWHM = nan(size(pkAll));
for j = 1:numel(pkAll)
    k = find(cellfun(@(v) any(abs(v - nmAll(j)) < 0.01), lineNm), 1);
    lineFWHM(j) = ex_fwhm(prLines(k, :), round(pkAll(j))) * disp_nm(round(pkAll(j)));
end

%% ---- report
fid = fopen('E1_results.txt', 'w');
out = @(varargin) cellfun(@(f) fprintf(f, varargin{:}), {1, fid});
out('E1 WAVELENGTH CALIBRATION AND RESOLUTION   (%s)\n', datestr(now));
out('Reference: %s   linearised: %d   strip length: %d px\n', ref.name, linearise, n);
out('Model: pixel = a + b / lambda^2   a = %.2f px, b = %.4g px nm^2\n', cal.p(1), cal.p(2));
out('\n  line (nm)   pixel     residual (nm)   FWHM (nm)\n');
for j = 1:numel(pkAll)
    out('  %8.1f   %7.1f   %+8.2f        %6.1f\n', nmAll(j), pkAll(j), cal.resid(j), lineFWHM(j));
end
if numel(pkAll) > 2
    out('  RMS residual: %.2f nm (%d lines, %d free parameters)\n', sqrt(mean(cal.resid.^2)), numel(pkAll), 2);
else
    out('  Only 2 lines: exact fit, no accuracy check. Add a third line.\n');
end
if ~isempty(calPoly)
    out('  Quadratic check: residuals %s nm\n', mat2str(round(100 * calPoly.resid) / 100));
end
use = Iref > 0.3 * max(Iref);
out('\nUsable range (reference > 30 %% of max): %.0f - %.0f nm; lines span %.0f - %.0f nm\n', min(lam(use)), max(lam(use)), min(nmAll), max(nmAll));
out('Dispersion (nm per pixel):\n');
for wl = [420 450 500 550 600 650 700]
    [d, i] = min(abs(lam - wl));
    if d < 5, out('  at %3d nm: %.3f nm/px\n', wl, disp_nm(i)); end
end
if ~isempty(res)
    out('\nResolution at %g nm (FWHM of the line):\n', resNm);
    for k = 1:numel(res)
        out('  %-28s %6.1f px = %5.1f nm\n', res(k).label, res(k).fwhmPx, res(k).fwhmNm);
    end
end
if ~isempty(satWarn)
    out('\nWARNING: saturated pixels on the line in %s -> peak too wide; use a shorter shutter time.\n', strjoin(unique(satWarn), ', '));
end
fclose(fid);

ex_writecsv('E1_lines.csv', {'lambda_nm', 'pixel', 'residual_nm', 'fwhm_nm'}, [nmAll(:), pkAll(:), cal.resid(:), lineFWHM(:)]);
refBlue = ref.blue;  refFileUsed = ref.file; %#ok<NASGU>
save('E1_calibration.mat', 'cal', 'strip', 'stripW', 'refBlue', 'refFileUsed', 'nmAll', 'pkAll', 'res', 'linearise');

%% ---- figure (draft Figure 4)
f = figure('Color', 'w', 'Position', [100 100 1100 380]);
subplot(1, 3, 1);  hold on;
plot(px, Iref / max(Iref), 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
for k = 1:numel(lineFiles), plot(px, prLines(k, :) / max(prLines(k, :)), 'LineWidth', 1.3); end
for j = 1:numel(pkAll), plot(pkAll([j j]), [0 1.05], 'k:'); end
xlabel('Position along strip (px, blue \rightarrow red)');  ylabel('Normalised intensity');
title('(a) Calibration lines');  box off;
subplot(1, 3, 2);  hold on;
pp = linspace(1, n, 300);
plot(psm_wavelength(cal, pp), pp, 'k-', 'LineWidth', 1.2);
plot(nmAll, pkAll, 'o', 'MarkerSize', 7, 'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerEdgeColor', 'w');
xlabel('Wavelength (nm)');  ylabel('Pixel');  title('(b) Prism model');  box off;
subplot(1, 3, 3);  hold on;
if ~isempty(res)
    for k = 1:numel(res)
        pr = res(k).profile;  i = max(1, res(k).pk - 60):min(n, res(k).pk + 60);
        plot(lam(i), (pr(i) - min(pr(i))) / (max(pr(i)) - min(pr(i))), 'LineWidth', 1.5);
    end
    legend({res.label}, 'Location', 'southoutside', 'Box', 'off');
    xlabel('Wavelength (nm)');  ylabel('Normalised intensity');  title('(c) Line width');
else
    stem(nmAll, cal.resid, 'filled');  plot([min(nmAll) max(nmAll)], [0 0], 'k:');
    xlabel('Wavelength (nm)');  ylabel('Residual (nm)');  title('(c) Residuals');
end
box off;
ex_savefig(f, 'E1_calibration');
fprintf('Saved E1_calibration.mat, E1_results.txt, E1_lines.csv and the figure in %s\n', pwd);
