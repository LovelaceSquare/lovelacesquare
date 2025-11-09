% test_cropBackground.m
%
% This script generates a synthetic hyperspectral cube containing a bright
% triangular region on a noisy background, then uses cropBackground to
% remove the background pixels based on total intensity thresholds.
%
% Author:  Adrián Gómez-Sánchez
% Date:    2024-12-14
% License: MIT
% Version: v1.0

%% 1) Setup parameters
rows      = 200;    % image height
cols      = 200;    % image width
nBands    = 30;     % number of spectral bands
noiseLevelBackground = 0.05;  % background noise std
noiseLevelSignal     = 0.01;  % signal noise std

% Wavelength axis (for plotting)
wavelengths = linspace(400, 900, nBands);  % nm

%% 2) Create triangular mask
% Define triangle vertices (x, y)
vx = [50, 150, 100];
vy = [150, 150, 50];
triangleMask = poly2mask(vx, vy, rows, cols);  % logical [rows x cols]

%% 3) Define spectral signature for the triangle
% Simulate a Gaussian-shaped reflectance spectrum
lambda0 = 650;                % peak at 650 nm
sigma   = 80/2.355;           % FWHM ~80 nm
spectralSig = exp( -((wavelengths - lambda0).^2) / (2*sigma^2) );
spectralSig = spectralSig / max(spectralSig);  % normalize to [0,1]

%% 4) Build the hyperspectral cube
hypCube = zeros(rows, cols, nBands);
for b = 1:nBands
    % Background noise everywhere
    bg = noiseLevelBackground * randn(rows, cols);
    % Signal only inside triangle, scaled by spectral signature + small noise
    sig = spectralSig(b) * ones(rows, cols);
    sig = sig + noiseLevelSignal * randn(rows, cols);
    sig(~triangleMask) = 0;  % zero outside triangle
    % Combine
    hypCube(:,:,b) = bg + sig;
end

%% 5) Compute global intensity map and pick thresholds
globalIntensity = sum(hypCube, 3);  % [rows x cols]
% Visualize intensity histogram to choose threshold
figure; imhist(globalIntensity(:), 100);
title('Histogram of Total Intensity');
xlabel('Sum Intensity'); ylabel('Pixel Count');

% Set thresholds to remove most background while keeping triangle
minT = 0.3 * max(globalIntensity(:));  
maxT = max(globalIntensity(:));

%% 6) Apply cropBackground
[croppedMatrix, keptIdx, removedIdx] = cropBackground(hypCube, minT, maxT);

%% 7) Reconstruct and visualize mask
mask = false(rows*cols,1);
mask(keptIdx) = true;
maskImage = reshape(mask, rows, cols);

figure;
subplot(2,2,1);
imshow(globalIntensity, []);
title('Global Intensity Map');

subplot(2,2,2);
imshow(maskImage);
title('Retained-Pixels Mask');

subplot(2,2,3);
% Show one representative band (e.g., band at index ~ mid-spectrum)
bShow = round(nBands/2);
imagesc(hypCube(:,:,bShow)); axis image off;
title(sprintf('Band %d (%.0f nm)', bShow, wavelengths(bShow)));

subplot(2,2,4);
imagesc(maskImage .* hypCube(:,:,bShow)); axis image off;
title('Masked Band (Signal Only)');

%% 8) Quantitative summary
nTotal    = rows*cols;
nKept     = numel(keptIdx);
nRemoved  = numel(removedIdx);
fprintf('Total pixels:   %d\n', nTotal);
fprintf('Kept pixels:    %d (%.2f%%)\n', nKept, 100*nKept/nTotal);
fprintf('Removed pixels: %d (%.2f%%)\n', nRemoved, 100*nRemoved/nTotal);

disp('Test of cropBackground completed. Review figures and summary.');
