%% e3_dilution.m  -  Experiment E3: linearity (dilution series), with and without linearisation
% Edit the SETTINGS block and press Run (MATLAB R2020b+ or GNU Octave).
%
% Every photo is processed twice: with the sRGB gamma undone (linearised)
% and with the encoded pixel values (as many school phone spectrometers do).
% Phone absorbance at the analytical wavelength is compared with the UV-Vis
% absorbance of the SAME cuvette.
%
% Photos: fixed (Pro/manual) exposure, phone and screen not moved. Run
% e1_calibration.m on calibration photos from the SAME session first: its
% strip and wavelength axis are used here (required).
%
% Output (in the current folder):
%   E3_results.txt        slopes, intercepts, R^2, departure from linearity, LOD
%   E3_points.csv         phone A (both modes, every photo) and UV-Vis A
%   E3_spectra.csv        mean phone spectra per concentration (both modes)
%   E3_linearity.png/.pdf draft of Figure 6

clear; close all;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'core'));  addpath(here);

%% ======================= SETTINGS =======================
folder   = fullfile(here, 'E3_photos');
e1File   = fullfile(here, 'E1_calibration.mat');   % from e1_calibration.m (required)
dyeName  = 'DCM in ethanol';
concUnit = 'µM';

darkFile   = 'dark.jpg';                           % '' = none
refFiles   = {'ref_01.jpg', 'ref_02.jpg', 'ref_03.jpg', 'ref_04.jpg'};   % solvent cuvette, interleaved
blankFiles = {'blank_1.jpg', 'blank_2.jpg', 'blank_3.jpg', 'blank_4.jpg', 'blank_5.jpg'};  % solvent cuvette as a "sample"

conc        = [2.5 6 12 18 24 36];               % concentrations, low -> high (DCM in ethanol: A(~470) ~ 0.1-1.5)
sampleFiles = { {'c1_1.jpg', 'c1_2.jpg', 'c1_3.jpg'}, ...
                {'c2_1.jpg', 'c2_2.jpg', 'c2_3.jpg'}, ...
                {'c3_1.jpg', 'c3_2.jpg', 'c3_3.jpg'}, ...
                {'c4_1.jpg', 'c4_2.jpg', 'c4_3.jpg'}, ...
                {'c5_1.jpg', 'c5_2.jpg', 'c5_3.jpg'}, ...
                {'c6_1.jpg', 'c6_2.jpg', 'c6_3.jpg'} };
uvFolder    = fullfile(here, 'E3_uvvis');
uvFiles     = {'c1.csv', 'c2.csv', 'c3.csv', 'c4.csv', 'c5.csv', 'c6.csv'};
uvBaseline  = [750 800];                           % nm, where the dye does not absorb ([] = none)
uvBlurNm    = 'auto';  % blur UV-Vis to the phone's resolution before comparing:
                       % 'auto' = FWHM measured in E1 (first resolution photo), a number (nm),
                       % or [] = no blur. Matters for bands with sharp structure; small for broad bands (DCM).

anaNm      = [];       % analytical wavelength (nm); [] = UV-Vis maximum of the strongest sample
halfWinNm  = 5;        % average +/- this around anaNm
fitMaxA    = 1.0;      % fit only points with UV-Vis A <= this
devLimit   = 0.10;     % "departs from linearity": |phone - fit| / fit above this
pathCm     = 1;        % cuvette path length (cm), for epsilon

alignOnBlue = false;   % follow small movements using the blue end of the rainbow. Keep FALSE
                       % for dyes that absorb blue light (DCM, yellow/orange dyes): they
                       % darken the blue end and the alignment goes wrong. The phone is
                       % clamped and E1 already requires that nothing moves.
minRefFrac  = 0.3;
smoothPx    = 5;
%% ========================================================

modes = {'linearised', 'encoded'};  lin = [true false];

%% ---- strip and wavelength axis
if exist(e1File, 'file') ~= 2
    error('Run e1_calibration.m first (same session, phone not moved): %s not found.', e1File);
end
E1 = load(e1File);
strip = E1.strip;  cal = E1.cal;  anchor.blue = E1.refBlue;
fprintf('Using calibration and strip from %s\n  %s\n', e1File, cal.info);
if strip.vertical, n = diff(round(strip.roi(1:2))) + 1; else, n = diff(round(strip.roi(3:4))) + 1; end
px = 1:n;  x = psm_wavelength(cal, px);

%% ---- UV-Vis
nC = numel(conc);  Auv = nan(1, nC);  uvSpec = cell(1, nC);
for c = 1:nC
    [wl, A] = psm_readSpectrumCSV(fullfile(uvFolder, uvFiles{c}));
    if ~isempty(uvBaseline)
        b = wl >= uvBaseline(1) & wl <= uvBaseline(2);
        if any(b), A = A - mean(A(b)); end
    end
    uvSpec{c} = [wl(:), A(:)];
end
% blur the UV-Vis spectra to the phone spectrometer's resolution
blurInfo = 'none';
if ischar(uvBlurNm) && strcmpi(uvBlurNm, 'auto')
    if isfield(E1, 'res') && ~isempty(E1.res) && ~isnan(E1.res(1).fwhmNm)
        uvBlurNm = E1.res(1).fwhmNm;  blurInfo = sprintf('%.1f nm FWHM (E1, %s)', uvBlurNm, E1.res(1).label);
    else
        uvBlurNm = [];  blurInfo = 'none (no resolution photo in E1)';
        warning('No resolution measured in E1: UV-Vis not blurred. Sharp bands will read lower on the phone.');
    end
elseif ~isempty(uvBlurNm)
    blurInfo = sprintf('%.1f nm FWHM (set by hand)', uvBlurNm);
end
uvSharp = uvSpec;
if ~isempty(uvBlurNm)
    for c = 1:nC, uvSpec{c} = ex_blurSpectrum(uvSpec{c}, uvBlurNm); end
end

if isempty(anaNm)
    w = uvSpec{end}(:, 1);  a = uvSpec{end}(:, 2);
    ok = w >= max(min(x), 400) & w <= min(max(x), 750);
    [~, i] = max(a .* ok);  anaNm = w(i);
end
for c = 1:nC
    w = uvSpec{c}(:, 1);  a = uvSpec{c}(:, 2);
    Auv(c) = mean(a(abs(w - anaNm) <= halfWinNm));
    w = uvSharp{c}(:, 1);  a = uvSharp{c}(:, 2);
    AuvSharp(c) = mean(a(abs(w - anaNm) <= halfWinNm)); %#ok<SAGROW>
end
win = abs(x - anaNm) <= halfWinNm;
if ~any(win), error('anaNm = %g nm is outside the calibrated strip.', anaNm); end

%% ---- phone photos, both modes
nR = max(cellfun(@numel, sampleFiles));
Aph   = nan(2, nC, nR);        % mode x concentration x repeat
Ablk  = nan(2, numel(blankFiles));
specM = zeros(2, nC, n);
warn  = {};
for m = 1:2
    % dark and references
    Id = 0;
    if ~isempty(darkFile)
        D  = psm_loadPhoto(fullfile(folder, darkFile), 'full', lin(m));
        Id = psm_profile(D, strip, [0 0], smoothPx);       % dark: nothing to align on
    end
    Ir = zeros(numel(refFiles), n);  refExp = nan(1, numel(refFiles));
    for r = 1:numel(refFiles)
        R = psm_loadPhoto(fullfile(folder, refFiles{r}), 'full', lin(m));
        Ir(r, :) = psm_profile(R, strip, psm_shift(R, anchor, alignOnBlue), smoothPx) - Id;
        refExp(r) = R.exposure;
    end
    Iref = mean(Ir, 1);
    if m == 1
        drift = (Ir(end, :) - Ir(1, :)) ./ Iref;
        refDrift = 100 * mean(drift(win));
    end
    % blanks
    for b = 1:numel(blankFiles)
        A = ex_phoneA(fullfile(folder, blankFiles{b}), lin(m), strip, anchor, alignOnBlue, smoothPx, Id, Iref, minRefFrac);
        Ablk(m, b) = mean(A(win), 'omitnan');
    end
    % samples
    for c = 1:nC
        S = zeros(numel(sampleFiles{c}), n);
        for r = 1:numel(sampleFiles{c})
            [A, P] = ex_phoneA(fullfile(folder, sampleFiles{c}{r}), lin(m), strip, anchor, alignOnBlue, smoothPx, Id, Iref, minRefFrac);
            S(r, :) = A;
            if m == 1 && ~isnan(P.exposure) && ~isnan(refExp(1)) && abs(P.exposure / refExp(1) - 1) > 0.02
                warn{end+1} = sprintf('%s: exposure %.4g s differs from reference %.4g s', P.name, P.exposure, refExp(1)); %#ok<SAGROW>
            end
            if m == 1 && any(P.sat(:)) && any(any(P.sat(round(strip.roi(1)):round(strip.roi(2)), round(strip.roi(3)):round(strip.roi(4)))))
                warn{end+1} = sprintf('%s: saturated pixels in the strip', P.name); %#ok<SAGROW>
            end
        end
        specM(m, c, :) = mean(S, 1);
        Aph(m, c, 1:size(S, 1)) = mean(S(:, win), 2, 'omitnan');
    end
end
%% ---- statistics
mA = mean(Aph, 3, 'omitnan');  sA = ex_nanstd(Aph, 3);
for m = 1:2
    use = Auv <= fitMaxA;
    F(m) = ex_linfit(Auv(use), mA(m, use)); %#ok<SAGROW>
    pred = F(m).m * Auv + F(m).c;
    dev  = (mA(m, :) - pred) ./ pred;
    k = find(abs(dev) > devLimit & Auv > 0.2, 1);
    if isempty(k), depart(m) = NaN; else, depart(m) = Auv(k); end %#ok<SAGROW>
    devAll(m, :) = dev; %#ok<SAGROW>
end
FuvC = ex_linfit(conc, AuvSharp);               % Beer-Lambert of the lab instrument (unblurred)
useS = AuvSharp <= fitMaxA;
Fsharp = ex_linfit(AuvSharp(useS), mA(1, useS));  % linearised phone vs UNblurred UV-Vis
sdBlank = ex_nanstd(Ablk(1, :));
LOD_A   = 3 * sdBlank / F(1).m;                 % in UV-Vis absorbance units
LOD_c   = LOD_A / FuvC.m;

%% ---- report
fid = fopen('E3_results.txt', 'w');
out = @(varargin) cellfun(@(f) fprintf(f, varargin{:}), {1, fid});
out('E3 LINEARITY (DILUTION SERIES)   %s   (%s)\n', dyeName, datestr(now));
out('%s\n', cal.info);
out('Analytical wavelength %.0f nm, +/- %g nm; fit for UV-Vis A <= %.2f\n', anaNm, halfWinNm, fitMaxA);
out('UV-Vis blurred to phone resolution: %s\n', blurInfo);
out('Reference drift first -> last photo: %+.1f %%\n\n', refDrift);
out('  conc (%s)   A UV-Vis   A phone lin (mean +/- SD)   A phone enc (mean +/- SD)   n\n', concUnit);
for c = 1:nC
    out('  %8.3g    %7.3f     %6.3f +/- %5.3f            %6.3f +/- %5.3f          %d\n', conc(c), Auv(c), ...
        mA(1, c), sA(1, c), mA(2, c), sA(2, c), sum(~isnan(Aph(1, c, :))));
end
out('\n  Blank (n = %d): lin %.4f +/- %.4f, enc %.4f +/- %.4f\n', sum(~isnan(Ablk(1, :))), ...
    mean(Ablk(1, :), 'omitnan'), sdBlank, mean(Ablk(2, :), 'omitnan'), ex_nanstd(Ablk(2, :)));
out('\nPhone vs UV-Vis, fit A_phone = m A_uv + c (points with A_uv <= %.2f):\n', fitMaxA);
for m = 1:2
    out('  %-10s  m = %.3f +/- %.3f   c = %+.3f +/- %.3f   R^2 = %.4f   n = %d\n', modes{m}, ...
        F(m).m, F(m).sm, F(m).c, F(m).sc, F(m).R2, F(m).n);
end
out('  slope ratio encoded / linearised = %.2f\n', F(2).m / F(1).m);
if ~isempty(uvBlurNm)
    out('  linearised vs UNblurred UV-Vis: m = %.3f  (the difference to the slope above is the resolution effect)\n', Fsharp.m);
end
for m = 1:2
    if isnan(depart(m)), out('  %-10s  within %.0f %% of the fit at all A_uv tested (max %.2f)\n', modes{m}, 100 * devLimit, max(Auv));
    else, out('  %-10s  departs by > %.0f %% from A_uv = %.2f\n', modes{m}, 100 * devLimit, depart(m)); end
end
out('\nUV-Vis Beer-Lambert: A = %.4f per %s (R^2 = %.4f)', FuvC.m, concUnit, FuvC.R2);
if strcmp(concUnit, 'µM') || strcmp(concUnit, 'uM')
    out('  ->  epsilon = %.0f M^-1 cm^-1 at %.0f nm', FuvC.m * 1e6 / pathCm, anaNm);
end
out('\nLOD (3 x SD blank / slope, linearised): A = %.4f  =  %.3g %s\n', LOD_A, LOD_c, concUnit);
if numel(blankFiles) < 5, out('  (fewer than 5 blank photos: LOD is a rough estimate)\n'); end
if ~isempty(warn), out('\nWARNINGS:\n'); out('  %s\n', warn{:}); end
fclose(fid);

hdr = {'conc', 'A_uv_blurred', 'A_uv'};  M = [conc(:), Auv(:), AuvSharp(:)];
for m = 1:2
    for r = 1:nR, hdr{end+1} = sprintf('A_%s_%d', modes{m}(1:3), r); M(:, end+1) = squeeze(Aph(m, :, r))'; end %#ok<SAGROW>
end
ex_writecsv('E3_points.csv', hdr, M);
hdr = {'pixel', 'x'};  M = [px(:), x(:)];
for m = 1:2
    for c = 1:nC, hdr{end+1} = sprintf('A_%s_c%g', modes{m}(1:3), conc(c)); M(:, end+1) = squeeze(specM(m, c, :)); end %#ok<SAGROW>
end
ex_writecsv('E3_spectra.csv', hdr, M);

%% ---- figure (draft Figure 6); all text outside the plot areas
cols = [linspace(0.90, 0.40, nC)', linspace(0.75, 0.05, nC)', linspace(0.90, 0.45, nC)'];  % light -> dark purple
f = figure('Color', 'w', 'Position', [100 100 1150 400]);
ax1 = subplot(1, 3, 1);  hold on;
ok = ~isnan(squeeze(specM(1, end, :)));
for c = 1:nC, plot(x(ok), squeeze(specM(1, c, ok)), 'Color', cols(c, :), 'LineWidth', 1.6); end
plot([anaNm anaNm], ylim, 'k:');
xlabel('Wavelength (nm)');  ylabel('Absorbance (phone, linearised)');  title('(a)');  box off;
legend(arrayfun(@(v) sprintf('%g %s', v, concUnit), conc, 'UniformOutput', false), ...
       'Location', 'eastoutside', 'Box', 'off');
subplot(1, 3, 2);  hold on;
lim = [0, 1.05 * max([Auv, mA(:)'])];
h1 = plot(lim, lim, '--', 'Color', [0.5 0.5 0.5]);
mk = {'o', 's'};  fc = {[0.10 0.30 0.65], 'w'};
for m = 1:2
    hE(m) = errorbar(Auv, mA(m, :), sA(m, :), mk{m});
    set(hE(m), 'Color', [0.10 0.30 0.65], 'MarkerFaceColor', fc{m}, 'MarkerSize', 7);
    xx = [0, fitMaxA];  plot(xx, F(m).m * xx + F(m).c, '-', 'Color', [0.10 0.30 0.65], 'LineWidth', 1 + (m == 1));
end
xlim(lim);  ylim(lim);  axis square;
xlabel('Absorbance, UV-Vis (at phone resolution)');  ylabel('Absorbance, phone');  title('(b)');  box off;
legend([h1, hE(1), hE(2)], {'1:1', sprintf('linearised (m = %.2f)', F(1).m), sprintf('encoded (m = %.2f)', F(2).m)}, ...
       'Location', 'southoutside', 'Box', 'off');
subplot(1, 3, 3);  hold on;
for m = 1:2
    hC(m) = plot(conc, mA(m, :), mk{m}, 'Color', [0.10 0.30 0.65], 'MarkerFaceColor', fc{m}, 'MarkerSize', 7);
end
hC(3) = plot(conc, AuvSharp, 'k^', 'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerSize', 6);
xlabel(sprintf('Concentration (%s)', concUnit));  ylabel('Absorbance');  title('(c)');  box off;
legend(hC, {'phone, linearised', 'phone, encoded', 'UV-Vis'}, 'Location', 'southoutside', 'Box', 'off');
ex_savefig(f, 'E3_linearity');
fprintf('Saved E3_results.txt, E3_points.csv, E3_spectra.csv and the figure in %s\n', pwd);
