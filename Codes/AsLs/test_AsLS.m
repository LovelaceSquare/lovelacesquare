% test_AsLS.m
%
% This script generates synthetic spectral data with known baselines,
% applies AsLS baseline correction, and evaluates the results visually
% and statistically.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 2024-12-16
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: 1.0

%% Setup Parameters
% Dimensions of synthetic data
nRows = 5;            % Number of spectra (rows)
nCols = 500;          % Number of data points per spectrum (columns)

% AsLS parameters
lambda = 8e5;         % Smoothing parameter for baseline smoothness
p = 0.06;            % Asymmetry parameter (typical range: 0.001 - 0.01)

% Noise and baseline parameters
baselineAmplitude = 0.4; % Amplitude of the baseline
noiseLevel = 0.05;       % Standard deviation of Gaussian noise

% Seed for reproducibility
rng(42);

%% Generate Synthetic Spectral Data
% Create x-axis values (e.g., wavelengths)
x = linspace(0, 10, nCols);

% Preallocate matrices for data and true baselines
trueBaseline = zeros(nRows, nCols);
data = zeros(nRows, nCols);

for i = 1:nRows
    % Generate synthetic spectrum with Gaussian peaks
    spectrum = ...
        exp(-((x - 3).^2) / (2 * 0.3^2)) + ...  % Narrow peak at x=3
        0.8 * exp(-((x - 6).^2) / (2 * 0.5^2)) + ... % Broader peak at x=6
        0.5 * exp(-((x - 8).^2) / (2 * 0.2^2));  % Narrow peak at x=8

    % Generate a baseline using a sinusoidal and linear trend
    baseline = baselineAmplitude * (0.5 * sin(2 * pi * x / 20) + 0.1 * x);

    % Add noise and baseline to the spectrum
    noisySpectrum = spectrum + baseline + noiseLevel * randn(1, nCols);

    % Store results
    trueBaseline(i, :) = baseline;
    data(i, :) = noisySpectrum;
end

%% Apply AsLS Baseline Correction
[correctedData, estimatedBaseline] = AsLS(data, lambda, p);

%% Visualization
% Select one spectrum for detailed visualization
spectrumIndex = 1;

figure;

% Original spectrum with true baseline
subplot(3, 1, 1);
plot(x, data(spectrumIndex, :), 'LineWidth', 1.5); hold on;
plot(x, trueBaseline(spectrumIndex, :), '--', 'LineWidth', 1.5);
title('Original Spectrum with True Baseline');
xlabel('Wavelength (arbitrary units)');
ylabel('Intensity');
legend('Noisy Spectrum', 'True Baseline');

% Estimated baseline
subplot(3, 1, 2);
plot(x, data(spectrumIndex, :), 'LineWidth', 1.5); hold on;
plot(x, estimatedBaseline(spectrumIndex, :), '--', 'LineWidth', 1.5);
title('Estimated Baseline');
xlabel('Wavelength (arbitrary units)');
ylabel('Intensity');
legend('Noisy Spectrum', 'Estimated Baseline');

% Corrected spectrum
subplot(3, 1, 3);
plot(x, correctedData(spectrumIndex, :), 'LineWidth', 1.5);
title('Baseline-Corrected Spectrum');
xlabel('Wavelength (arbitrary units)');
ylabel('Intensity');
legend('Corrected Spectrum');

%% Statistical Evaluation
% Compute Mean Squared Error (MSE) for baseline estimation
mseBaseline = mean((trueBaseline(:) - estimatedBaseline(:)).^2);

% Display MSE result
fprintf('Mean Squared Error for Estimated Baseline: %.6e\n', mseBaseline);

% Completion message
disp('AsLS baseline correction test completed. Review plots and results.');
