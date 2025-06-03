% testWienerFiltering.m
%
% This script generates synthetic spectral data, adds noise, applies
% Wiener filtering, and compares the original, noisy, and filtered data
% for analysis.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 2024-12-14
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0
%

%% Parameters Setup
% Spectral data dimensions
nSamples = 10;         % Number of spectra (rows)
nWavelengths = 500;    % Number of spectral points (columns)

% Wiener filtering parameters
noiseVar = 0.01;       % Estimated noise variance
segmentLength = 128;   % Segment length for Welch's PSD
overlap = 64;          % Overlap between segments
nfft = 512;            % FFT length (must be >= nWavelengths)

% Set a random seed for reproducibility
rng(42);

%% Generate Synthetic Spectral Data
% Generate sine and cosine-based spectral data with varying frequencies
t = linspace(0, 2 * pi, nWavelengths);
data = zeros(nSamples, nWavelengths);
for i = 1:nSamples
    data(i, :) = sin((i + 1) * t) + cos((i + 2) * t);
end

% Add Gaussian noise to the data
noisyData = data + sqrt(noiseVar) * randn(size(data));

%% Apply Wiener Filtering
% Call the WienerFiltering function
filteredData = WienerFiltering(noisyData, noiseVar, segmentLength, overlap, nfft);

%% Visualization
% Plot the original, noisy, and filtered data for a single spectrum
spectrumIndex = 1;  % Choose a spectrum to visualize

figure;

% Original spectrum
subplot(3, 1, 1);
plot(t, data(spectrumIndex, :), 'LineWidth', 1.5);
title(sprintf('Original Spectrum (Sample %d)', spectrumIndex));
xlabel('Wavelength (arbitrary units)');
ylabel('Amplitude');

% Noisy spectrum
subplot(3, 1, 2);
plot(t, noisyData(spectrumIndex, :), 'LineWidth', 1.5);
title(sprintf('Noisy Spectrum (Sample %d)', spectrumIndex));
xlabel('Wavelength (arbitrary units)');
ylabel('Amplitude');

% Filtered spectrum
subplot(3, 1, 3);
plot(t, filteredData(spectrumIndex, :), 'LineWidth', 1.5);
title(sprintf('Filtered Spectrum (Sample %d)', spectrumIndex));
xlabel('Wavelength (arbitrary units)');
ylabel('Amplitude');

% Display results
disp('Wiener filtering applied successfully. Results plotted for inspection.');

%% Statistical Comparison
% Compute mean squared error (MSE) before and after filtering
originalMSE = mean((noisyData(:) - data(:)).^2);
filteredMSE = mean((filteredData(:) - data(:)).^2);

% Display MSE results
fprintf('Mean Squared Error (Noisy vs. Original): %.6f\n', originalMSE);
fprintf('Mean Squared Error (Filtered vs. Original): %.6f\n', filteredMSE);
