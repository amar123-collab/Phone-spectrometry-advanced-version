function P = psm_loadPhoto(file, mode, linearise)
%PSM_LOADPHOTO  Load a phone photo for spectrum analysis.
%   P = psm_loadPhoto(file)                  full analysis load
%   P = psm_loadPhoto(file, 'preview')       thumbnail only (fast)
%   P = psm_loadPhoto(file, 'full', false)   full load WITHOUT undoing the gamma
%
%   Applies the phone's rotation flag (imread ignores it), undoes the JPEG
%   (sRGB) gamma so intensities are approximately linear, and records
%   exposure time and ISO from the photo metadata.
%
%   linearise (default true): false keeps the encoded pixel values (0-1).
%   Use it only to show what the gamma does to absorbance (e.g. a dilution
%   series processed both ways); real measurements should be linearised.
%
%   Fields: name, file, thumb, thumbStep, exposure, iso, and (full load only)
%   I (intensity, sum of R+G+B), linearised (true/false), sat (saturated
%   pixels), blue / red (centroids [row col] of the blue and red ends of the
%   rainbow), medianI.

if nargin < 2 || isempty(mode), mode = 'full'; end
if nargin < 3 || isempty(linearise), linearise = true; end
thumbStep = 4;

raw = imread(file);
[~, n, e] = fileparts(file);
P.name = [n e];
P.file = file;
P.exposure = NaN;
P.iso = NaN;

ws = warning('off', 'all');                 % imfinfo warns about some phone tags
try
    info = imfinfo(file);  info = info(1);
    if isfield(info, 'Orientation')
        switch info.Orientation
            case 3, raw = rot90(raw, 2);
            case 6, raw = rot90(raw, -1);
            case 8, raw = rot90(raw, 1);
        end
    end
    if isfield(info, 'DigitalCamera')
        dc = info.DigitalCamera;
        if isfield(dc, 'ExposureTime'),    P.exposure = double(dc.ExposureTime); end
        if isfield(dc, 'ISOSpeedRatings'), P.iso = double(dc.ISOSpeedRatings(1)); end
    end
catch
end
warning(ws);

if size(raw, 3) == 1, raw = repmat(raw, 1, 1, 3); end
raw = raw(:, :, 1:3);
if isa(raw, 'uint16'), sc = 65535; else, sc = 255; end

t = raw(1:thumbStep:end, 1:thumbStep:end, :);
if ~isa(t, 'uint8'), t = uint8(double(t) / sc * 255); end
P.thumb = t;
P.thumbStep = thumbStep;
if strcmpi(mode, 'preview'), return; end

P.sat = any(raw >= 0.98 * sc, 3);
v = single(raw) / sc;
clear raw
lin = (v <= 0.04045) .* (v / 12.92) + (v > 0.04045) .* (((v + 0.055) / 1.055) .^ 2.4);
P.linearised = logical(linearise);
if P.linearised, P.I = sum(lin, 3); else, P.I = sum(v, 3); end
clear v
[P.blue, P.red] = colourCentroids(lin);        % landmarks always from linear values
s = P.I(1:8:end, 1:8:end);
P.medianI = median(s(:));
end

% -------------------------------------------------------------------------
function [b, r] = colourCentroids(lin)
R = lin(:, :, 1);  G = lin(:, :, 2);  B = lin(:, :, 3);
b = weightedCentroid(B .* (B > 1.6*R & B > 1.3*G & B > 0.5*highValue(B)));
r = weightedCentroid(R .* (R > 1.6*B & R > 1.3*G & R > 0.5*highValue(R)));
end

function v = highValue(X)
s = sort(X(1:4:end));                       % 99.9th percentile, no toolbox needed
v = s(max(1, round(0.999 * numel(s))));
end

function c = weightedCentroid(w)
w = double(w);  t = sum(w(:));
if t <= 0, c = [NaN NaN]; return; end
c = [sum((1:size(w, 1))' .* sum(w, 2)) / t, sum((1:size(w, 2)) .* sum(w, 1)) / t];
end
