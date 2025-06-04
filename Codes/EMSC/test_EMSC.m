% test_EMSC.m
% This script tests the functionality of the EMSC function by simulating
% spectral data with known additive and multiplicative effects.

% Author: Adrián Gómez-Sánchez
% Date Created: 2024-12-14
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0

%% Parameters for Simulated Spectral Data
nSamples = 20;               % Number of spectra
nWavelengths = 500;          % Number of wavelengths per spectrum

trueSpectra = linspace(0, 1, nWavelengths);  % Underlying spectrum (linear)
trueSpectra = repmat(trueSpectra, nSamples, 1); % Replicate for all samples

% Add a sine wave component for more realistic spectra
for i = 1:nSamples
    trueSpectra(i, :) = trueSpectra(i, :) + sin(2 * pi * (1:nWavelengths) / nWavelengths) * (0.1 * i);
end

% Simulate multiplicative and additive effects
multiplicativeFactors = rand(nSamples, 1) * 0.5; % Between 0.75 and 1.25
additiveOffsets = abs(randn(nSamples, 1)) * 0.5;             % Small random offsets

% Apply effects to generate raw spectral data
rawSpectra = trueSpectra .* multiplicativeFactors + additiveOffsets;

%% Apply EMSC
refType = 'Mean';   % Use the mean spectrum as the reference
polyOrder = 1;      % Linear baseline modeling
[correctedSpectra, referenceSpec] = EMSC(rawSpectra, refType, polyOrder);

%% Visualize Results
figure;
subplot(3, 1, 1);
plot(1:nWavelengths, rawSpectra');
title('Raw Spectra with Additive and Multiplicative Effects');
xlabel('Wavelength'); ylabel('Intensity');
grid on;

subplot(3, 1, 2);
plot(1:nWavelengths, referenceSpec, 'LineWidth', 2);
title('Reference Spectrum (Mean or Median)');
xlabel('Wavelength'); ylabel('Intensity');
grid on;

subplot(3, 1, 3);
plot(1:nWavelengths, correctedSpectra');
title('Corrected Spectra After EMSC');
xlabel('Wavelength'); ylabel('Intensity');
grid on;
