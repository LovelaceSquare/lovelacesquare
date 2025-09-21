%% Test Script for MSC Function
% This script tests the functionality of the MSC function
%   • Mean‐spectrum reference
%   • Single-row reference (index)
%   • Explicit reference spectrum (row vector)
%
% Author:  Adrián Gómez-Sánchez
% Date Created: 2024-12-14
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v1.3   % (2025-07-14) – add explicit-spectrum test case

%% 1) Parameters for Simulated NIR Spectral Data
nSamples      = 20;               % Number of spectra (rows)
nWavelengths  = 1000;             % Number of wavelengths (columns)
wavelengthMin = 800;              % Minimum wavelength in nm
wavelengthMax = 2500;             % Maximum wavelength in nm
wavelengths   = linspace(wavelengthMin, wavelengthMax, nWavelengths);

%% 2) Generate Realistic "True" NIR Spectra
peakPositions  = [950, 1200, 1400, 1900, 2100, 2300];  % Typical NIR regions
peakAmplitudes = [1, 0.8, 1.2, 0.9, 1.1, 0.7];
peakWidths     = [50, 60, 55, 70, 65, 60];

trueSpectra = zeros(nSamples, nWavelengths);
for i = 1:nSamples
    spectrum = zeros(1, nWavelengths);
    for p = 1:length(peakPositions)
        ampVar  = peakAmplitudes(p)*(1 + 0.05*randn());         % ±5 %
        spectrum = spectrum + ampVar .* ...
            exp(-0.5*((wavelengths - peakPositions(p))/peakWidths(p)).^2);
    end
    trueSpectra(i, :) = spectrum;
end

%% 3) Simulate Multiplicative and Additive Effects
multiplicativeFactors = 0.85 + (rand(nSamples,1)*2);  % 0.85–2.85
additiveOffsets       = rand(nSamples,1)*3;           % 0–3

rawSpectra = (trueSpectra .* multiplicativeFactors) + additiveOffsets;

%% 5) Apply MSC
% 5.1) Mean spectrum reference
[correctedMean, referenceMean] = MSC(rawSpectra, 'Mean Spectrum', 1);

% 5.2) Specific sample reference (index)
refIndex = 10;
[correctedRefIdx, referenceSpecRef] = MSC(rawSpectra, 'Reference Index', refIndex);

% 5.3) Explicit external reference spectrum -----------------------------
customRef = trueSpectra(7, :);            % here we use sample #7 as “gold-standard”
[correctedCustom, referenceSpecCustom] = ...
        MSC(rawSpectra, 'Reference Spectrum', customRef);

%% 6) Verification of MSC Correction – Mean Squared Error (MSE)
mseMean   = mean((correctedMean   - trueSpectra).^2, 'all');
mseRefIdx = mean((correctedRefIdx - trueSpectra).^2, 'all');
mseCustom = mean((correctedCustom - trueSpectra).^2, 'all');

fprintf('MSE  Mean Spectrum reference:        %.4e\n', mseMean);
fprintf('MSE  Reference Index  %2d:           %.4e\n', refIndex, mseRefIdx);
fprintf('MSE  Explicit Reference Spectrum:    %.4e\n', mseCustom);

%% 7) Visualize Results
% (a) True Spectra
subplot(3,1,1);
plot(wavelengths, trueSpectra', 'LineWidth', 1.2);
title('True Spectra (Underlying Signal)');
xlabel('Wavelength (nm)'); ylabel('Intensity (a.u.)');
xlim([wavelengthMin, wavelengthMax]); grid on;

% (b) Raw Spectra with Additive and Multiplicative Effects
subplot(3,1,2);
plot(wavelengths, rawSpectra', 'LineWidth', 1.2);
title('Raw Spectra with Additive and Multiplicative Effects');
xlabel('Wavelength (nm)'); ylabel('Intensity (a.u.)');
xlim([wavelengthMin, wavelengthMax]); grid on;

% (c) Corrected Spectra After MSC – Mean Reference
subplot(3,1,3);
plot(wavelengths, correctedMean', 'LineWidth', 1.2);
title('Corrected Spectra After MSC (Mean Reference)');
xlabel('Wavelength (nm)'); ylabel('Intensity (a.u.)');
xlim([wavelengthMin, wavelengthMax]); grid on;

%% Overlays for one random sample
sampleIdx = randi(nSamples);   % pick a spectrum to inspect

figure;
plot(wavelengths, trueSpectra(sampleIdx,:),    'k-', 'LineWidth', 2); hold on;
plot(wavelengths, rawSpectra(sampleIdx,:),     'c--','LineWidth', 1.6);
plot(wavelengths, correctedMean(sampleIdx,:),  'r-.','LineWidth', 1.6);
plot(wavelengths, correctedRefIdx(sampleIdx,:), 'm:','LineWidth', 1.6);
plot(wavelengths, correctedCustom(sampleIdx,:), 'g-','LineWidth', 1.2);
hold off;

title(sprintf('Sample #%d – True vs. Raw vs. MSC Variants', sampleIdx));
xlabel('Wavelength (nm)'); ylabel('Intensity (a.u.)');
legend({'True','Raw','MSC Mean','MSC RefIdx', 'MSC Custom'}, ...
       'Location','Best');
xlim([wavelengthMin, wavelengthMax]); grid on;
