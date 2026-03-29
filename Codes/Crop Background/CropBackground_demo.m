% CropBackground_demo. Generate realistic NIR hyperspectral cube from a plastic sample image.
%
% Reads demo_plastics.png (RGB photo of plastic pieces on a dark rubber mat),
% classifies each pixel by material, and assigns physically realistic NIR
% reflectance spectra modulated by the original pixel intensity. All spatial
% texture (background mat grain, cap highlights, plastic translucency) is
% preserved in the resulting hyperspectral cube.
%
% Materials:
%   Blue cap      — PP (polypropylene): C-H overtones at 1195, 1395, 1720 nm
%   Green cap     — HDPE (high-density polyethylene): C-H at 1210, 1400, 1730 nm
%   Clear plastic — PET (polyethylene terephthalate): aromatic C-H, ester bands
%   Black plastic — ABS + carbon black: broad absorption, very low reflectance
%   Background    — Dark rubber mat: low reflectance with full surface texture
%
% Output variables:
%   cube       — [rows x cols x 60] NIR reflectance cube (900–1750 nm)
%   wavelength — [1 x 60] wavelength axis in nm
%
% Usage:
%   CropBackground_demo       % creates 'cube' and 'wavelength' in workspace
%   app = CropBackground(cube);
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-29
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0

%% ---- Read image ----
imgPath = fullfile(fileparts(mfilename('fullpath')), 'demo_plastics.png');
if ~isfile(imgPath)
    error('CropBackground_demo:fileNotFound', ...
        'demo_plastics.png not found. Place it next to this script.');
end

rgb = im2double(imread(imgPath));

% Downsample for manageable cube size
scaleFactor = 0.3;
rgb = imresize(rgb, scaleFactor, 'bilinear');
[nRows, nCols, ~] = size(rgb);

fprintf('Image size: %d x %d\n', nRows, nCols);

%% ---- Spectral axis ----
nBands = 60;
wavelength = linspace(900, 1750, nBands);  % nm

% Gaussian peak helper
gauss = @(x, mu, A, sigma) A .* exp(-((x - mu).^2) ./ (2 * sigma.^2));

%% ---- Define NIR spectral SHAPES (normalised templates) ----
% These are the spectral fingerprints. Each pixel's actual intensity comes
% from the original RGB image, preserving all texture and shading.

% PP (polypropylene) — blue cap
% C-H 2nd overtone ~1195, combination ~1395, 1st overtone ~1720 nm
tPP = 1.0 + gauss(wavelength, 1195, -0.35, 25) ...
          + gauss(wavelength, 1395, -0.40, 20) ...
          + gauss(wavelength, 1720, -0.55, 28) ...
          + gauss(wavelength, 1160, -0.15, 15);

% HDPE (high-density polyethylene) — green cap
tHDPE = 1.0 + gauss(wavelength, 1210, -0.30, 22) ...
            + gauss(wavelength, 1400, -0.38, 22) ...
            + gauss(wavelength, 1730, -0.50, 30) ...
            + gauss(wavelength, 1175, -0.12, 18);

% PET (polyethylene terephthalate) — clear plastic
tPET = 1.0 + gauss(wavelength, 1130, -0.22, 20) ...
           + gauss(wavelength, 1410, -0.28, 18) ...
           + gauss(wavelength, 1680, -0.45, 25) ...
           + gauss(wavelength, 1340, -0.15, 22) ...
           + gauss(wavelength, 1510, -0.10, 15);

% ABS + carbon black — broad featureless absorption
tBlack = 1.0 + gauss(wavelength, 1180, -0.05, 40) ...
             + gauss(wavelength, 1690, -0.07, 45) ...
             - 0.08 * (wavelength - 900) / 850;

% Background rubber mat — nearly flat, slight rubber signature
tBG = 1.0 + gauss(wavelength, 1150, -0.06, 30) ...
          + gauss(wavelength, 1650, -0.04, 40) ...
          + 0.03 * sin(wavelength * 0.008);

% Normalise templates to [0, 1] range
tPP    = max(0.05, tPP / max(tPP));
tHDPE  = max(0.05, tHDPE / max(tHDPE));
tPET   = max(0.05, tPET / max(tPET));
tBlack = max(0.05, tBlack / max(tBlack));
tBG    = max(0.05, tBG / max(tBG));

%% ---- Soft material classification per pixel ----
% Use HSV to compute per-pixel membership weights for each material.
% Weights are soft (continuous), so boundaries blend naturally.
hsv = rgb2hsv(rgb);
H = hsv(:,:,1);
S = hsv(:,:,2);
V = hsv(:,:,3);
gray = rgb2gray(rgb);  % luminance for intensity modulation

% Weight functions (soft Gaussians in colour space)
softH = @(h, mu, sigma) exp(-min((h-mu).^2, min((h-mu+1).^2, (h-mu-1).^2)) / (2*sigma^2));

wBlue  = softH(H, 0.60, 0.06) .* (S > 0.20) .* (V > 0.25);
wGreen = softH(H, 0.35, 0.06) .* (S > 0.20) .* (V > 0.20);
wClear = max(0, 1 - S/0.15) .* max(0, (V - 0.40)/0.3);  % low S, high V
wBlack = max(0, 1 - V/0.22) .* max(0, 1 - S/0.20);       % very dark, low S

% Background: everything not strongly assigned to an object
wObj   = wBlue + wGreen + wClear + wBlack;
wBG    = max(0, 1 - wObj);

% Normalise so weights sum to 1 per pixel
wTotal = wBlue + wGreen + wClear + wBlack + wBG + 1e-10;
wBlue  = wBlue  ./ wTotal;
wGreen = wGreen ./ wTotal;
wClear = wClear ./ wTotal;
wBlack = wBlack ./ wTotal;
wBG    = wBG    ./ wTotal;

%% ---- Build hyperspectral cube ----
% For each band, blend the spectral templates weighted by material membership,
% then modulate by the original pixel intensity (grayscale) to preserve texture.
cube = zeros(nRows, nCols, nBands);

% Base intensity from the original image (preserves ALL texture)
% Scale so bright pixels are high reflectance, dark pixels are low
baseIntensity = gray;  % 0..1

for b = 1:nBands
    % Blended spectral shape per pixel
    spectralShape = wBlue  * tPP(b) ...
                  + wGreen * tHDPE(b) ...
                  + wClear * tPET(b) ...
                  + wBlack * tBlack(b) ...
                  + wBG    * tBG(b);

    % Modulate by original intensity — this preserves all spatial texture:
    % highlights on caps, grain on rubber mat, translucency in PET, etc.
    bandData = baseIntensity .* spectralShape;

    % Realistic noise model
    shotNoise     = sqrt(max(0.001, bandData)) .* randn(nRows, nCols) * 0.008;
    detectorNoise = 0.003 * randn(nRows, nCols);
    cube(:,:,b)   = max(0, bandData + shotNoise + detectorNoise);
end

%% ---- Output ----
clearvars -except cube wavelength

fprintf('Created: cube (%dx%dx%d)\n', size(cube,1), size(cube,2), size(cube,3));
fprintf('  Wavelength range: %.0f – %.0f nm (NIR)\n', wavelength(1), wavelength(end));
fprintf('  Materials: PP (blue), HDPE (green), PET (clear), ABS (black), rubber (background)\n');
fprintf('  Spatial texture fully preserved from original photograph\n');
fprintf('  Background (dark mat) has lowest NIR signal → ideal for cropping\n');
fprintf('\nRun:  app = CropBackground(cube);\n');
