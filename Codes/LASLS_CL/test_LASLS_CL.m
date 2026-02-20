% test_LASLS_CL.m
%
% This script generates synthetic spectral data with known baselines,
% applies LASLS_CL baseline correction with local multi-interval parameters,
% and evaluates the results visually and statistically.
%
% REFERENCES:
%   Eilers, Paul H.C., and Hans F.M. Boelens.
%   "Baseline correction with asymmetric least squares smoothing."
%   Leiden University Medical Centre Report 1.1 (2005): 5.
%
% Authors: Adrián Gómez-Sánchez, Berta Torres-Cobos, Rodrigo Rocha de Oliveira
% Date Created: 2024-12-16
% License: MIT
% Repository: https://github.com/LovelaceSquare/lovelacesquare
% Reviewed by Lovelace's Square: Yes
% Version: 1.1
%
% The script creates synthetic spectral data containing multiple Gaussian peaks
% with an added baseline and Gaussian noise. It then applies the Local Asymmetric
% Least Squares (LASLS_CL) baseline correction algorithm to a selected spectrum using
% specified intervals with unique asymmetry and smoothing parameters. The resulting
% estimated baseline is subtracted from the original spectrum to obtain a baseline-
% corrected spectrum. Visual plots and a statistical evaluation (Mean Squared Error)
% are provided for assessment.

%% Setup Parameters
% Dimensions of synthetic data
nRows = 5;            % Number of spectra (rows)
nCols = 500;          % Number of data points per spectrum (columns)

% LASLS_CL parameters for a single spectrum (test on one row)
% Define intervals around expected peaks (convert x-values to indices):
% For peaks near x=3, x=6, and x=8.
intervals = [140, 160; 280, 320; 375, 425];

pVals = [0.0001; 0.0001; 0.0005];       % Asymmetry parameters for each interval
lambdasAsym = [1e5; 1e5; 1e5];        % Local smoothing penalties for each interval
lambdaWhit = 10;                   % Smoothing penalty outside intervals
mu = 10;                           % Global first-derivative penalty
maxIter = 50;                       % Maximum iterations for IRLS
tol = 1e-6;                         % Convergence tolerance

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
        exp(-((x - 3).^2) / (0.1^2)) + ...         % Narrow peak at x=3
        0.8 * exp(-((x - 6).^2) / (0.2^2)) + ...     % Broader peak at x=6
        0.5 * exp(-((x - 8).^2) / (0.2^2));           % Narrow peak at x=8

    % Generate a baseline using a sinusoidal and linear trend
    baseline = baselineAmplitude * (0.5 * sin(2 * pi * x / 20) + 0.1 * x);

    % Add noise and baseline to the spectrum
    noisySpectrum = spectrum + baseline + noiseLevel * randn(1, nCols);

    % Store results
    trueBaseline(i, :) = baseline;
    data(i, :) = noisySpectrum;
end

%% Apply LASLS_CL Baseline Correction on a Single Spectrum
% Select one spectrum for detailed analysis
spectrumIndex = 1;
y = data(spectrumIndex, :)';

% Apply LASLS_CL baseline correction
[estimatedBaseline, weights] = LASLS_CL(y, intervals, pVals, lambdasAsym, lambdaWhit, mu, maxIter, tol);

% Compute baseline-corrected spectrum
correctedSpectrum = y - estimatedBaseline;

%% Visualization
figure;

% Original spectrum with true baseline
subplot(3, 1, 1);
plot(x, data(spectrumIndex, :), 'LineWidth', 1.5); hold on;
plot(x, trueBaseline(spectrumIndex, :), '--', 'LineWidth', 1.5);
title('Original Spectrum with True Baseline');
xlabel('Wavelength (arbitrary units)');
ylabel('Intensity');
legend('Noisy Spectrum', 'True Baseline');

% Estimated baseline from LASLS_CL
subplot(3, 1, 2);
plot(x, data(spectrumIndex, :), 'LineWidth', 1.5); hold on;
plot(x, estimatedBaseline, '--', 'LineWidth', 1.5);
title('Estimated Baseline (LASLS_CL)');
xlabel('Wavelength (arbitrary units)');
ylabel('Intensity');
legend('Noisy Spectrum', 'Estimated Baseline');

% Baseline-corrected spectrum
subplot(3, 1, 3);
plot(x, correctedSpectrum, 'LineWidth', 1.5);
title('Baseline-Corrected Spectrum');
xlabel('Wavelength (arbitrary units)');
ylabel('Intensity');
legend('Corrected Spectrum');

%% Statistical Evaluation
% Compute Mean Squared Error (MSE) for baseline estimation
mseBaseline = mean((trueBaseline(spectrumIndex, :)' - estimatedBaseline).^2);

% Display MSE result
fprintf('Mean Squared Error for Estimated Baseline: %.6e\n', mseBaseline);

% Completion message
disp('LASLS_CL baseline correction test completed. Review plots and results.');
