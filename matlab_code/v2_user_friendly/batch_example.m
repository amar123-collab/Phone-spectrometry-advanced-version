%% batch_example.m  -  phone spectra from a fixed list of photos (no clicks)
% Edit the SETTINGS block and press Run. Works out of the box on the
% example photos in ../example_data/photos.
%
% Result: absorbance spectra of each sample photo relative to the reference
% photo, plotted and saved as CSV + PNG in the current folder.

clear; close all;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'core'));

%% ======================= SETTINGS =======================
folder      = fullfile(here, '..', 'example_data', 'photos');
refFile     = '01_reference_ethanol.jpg';   % ethanol in the cuvette (or no sample)
darkFile    = '00_dark.jpg';                % lamp blocked, same camera settings; '' = none
sampleFiles = {'02_DCM_49uM.jpg', '03_sample1.jpg', '04_sample2.jpg', '05_sample3.jpg', '06_sample4.jpg', '07_sample5.jpg'};
sampleNames = {'DCM 49.5 uM', 'Sample 1 (violet)', 'Sample 2 (cyan)', 'Sample 3 (magenta)', 'Sample 4 (pink)', 'Sample 5 (lime)'};

laserFile   = '';                           % photo of laser spots for calibration; '' = pixel axis
laserNm     = [405 532 650];                % wavelengths of those spots (any order)
dispersion  = 'prism';                      % 'prism' (curved fit) or 'grating' (straight line)

offsetMode  = 'none';                       % 'none' for fixed (Pro-mode) exposure,
                                            % 'auto' for automatic-exposure photos
alignOnBlue = false;                        % follow small movements between photos (false: phone clamped)
minRefFrac  = 0.3;                          % ignore dim ends (< 30 % of reference maximum)
smoothPx    = 5;                            % smoothing along the rainbow (pixels)
linearise   = true;                         % undo the JPEG gamma (keep true; false only
                                            % to show the effect of skipping it)
outName     = 'phone_spectra';              % output file prefix
%% ========================================================

ref   = psm_loadPhoto(fullfile(folder, refFile), 'full', linearise);
strip = psm_findStrip(ref);
Iref  = psm_profile(ref, strip, [0 0], smoothPx);
Idark = 0;
if ~isempty(darkFile)
    dark  = psm_loadPhoto(fullfile(folder, darkFile), 'full', linearise);
    Idark = psm_profile(dark, strip, psm_shift(dark, ref, alignOnBlue), smoothPx);
    Iref  = Iref - Idark;
end
px = 1:numel(Iref);

% wavelength axis
if ~isempty(laserFile)
    L   = psm_loadPhoto(fullfile(folder, laserFile), 'full', linearise);
    pk  = psm_findPeaks(psm_profile(L, strip, psm_shift(L, ref, alignOnBlue), smoothPx), numel(laserNm));
    if numel(pk) < numel(laserNm), error('Found only %d laser spots, expected %d.', numel(pk), numel(laserNm)); end
    cal = psm_calibration('laser', pk, sort(laserNm), 1 + strcmpi(dispersion, 'prism'));
else
    cal = psm_calibration('none');
end
fprintf('%s\n', cal.info);
if ~linearise, fprintf('NOTE: gamma NOT undone (linearise = false) - absorbance values are distorted.\n'); end
x = psm_wavelength(cal, px);
if strcmp(cal.type, 'none'), xlab = 'Position along the rainbow (pixels, blue \rightarrow red)';
else,                        xlab = 'Wavelength (nm)'; end

% samples
nS = numel(sampleFiles);
Is = zeros(nS, numel(px));  A = Is;  Araw = Is;
for k = 1:nS
    P = psm_loadPhoto(fullfile(folder, sampleFiles{k}), 'full', linearise);
    [sh, ok] = psm_shift(P, ref, alignOnBlue);
    if ~ok, warning('%s: no blue spot found, not aligned.', P.name); end
    if ~isnan(P.exposure) && ~isnan(ref.exposure) && abs(P.exposure / ref.exposure - 1) > 0.05
        fprintf('  note: %s exposure %.4g s differs from reference %.4g s\n', P.name, P.exposure, ref.exposure);
    end
    Is(k, :)   = psm_profile(P, strip, sh, smoothPx) - Idark;
    Araw(k, :) = psm_absorbance(Is(k, :), Iref, minRefFrac);
    A(k, :)    = psm_offset(Araw(k, :), offsetMode);
    [am, i] = max(A(k, :));
    fprintf('%-22s peak at %6.0f %s   A = %.2f\n', sampleNames{k}, x(i), ...
            ifelse(strcmp(cal.type, 'none'), 'px', 'nm'), am);
end

% figures
f1 = figure('Color', 'w', 'Name', 'Light along the rainbow');  hold on;
plot(x, Iref, 'k', 'LineWidth', 2, 'DisplayName', 'Reference');
for k = 1:nS, plot(x, Is(k, :), 'LineWidth', 1.4, 'DisplayName', sampleNames{k}); end
xlabel(xlab);  ylabel('Intensity (a.u.)');  legend('Box', 'off');  box off;

f2 = figure('Color', 'w', 'Name', 'Absorbance');  hold on;
for k = 1:nS, plot(x, A(k, :), 'LineWidth', 2, 'DisplayName', sampleNames{k}); end
yline(0, ':', 'HandleVisibility', 'off');
xlabel(xlab);  ylabel('Absorbance');  legend('Box', 'off', 'Location', 'best');  box off;

% save
out = table(px(:), x(:), Iref(:), 'VariableNames', {'pixel', 'x', 'I_ref'});
for k = 1:nS
    out.(sprintf('I_S%d', k))    = Is(k, :)';
    out.(sprintf('Araw_S%d', k)) = Araw(k, :)';
    out.(sprintf('A_S%d', k))    = A(k, :)';
end
writetable(out, [outName '.csv']);
exportgraphics(f1, [outName '_intensity.png'], 'Resolution', 300);
exportgraphics(f2, [outName '_absorbance.png'], 'Resolution', 300);
fprintf('Saved %s.csv and PNG figures in %s\n', outName, pwd);

function s = ifelse(c, a, b)
if c, s = a; else, s = b; end
end
