% test_whittakerSmoother.m

% This script generates synthetic spectral data with noise, applies the
% Whittaker smoother, and evaluates the results visually and statistically.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 2024-12-14
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0

%% Setup Parameters
% Define the dimensions of the synthetic data
nRows = 5;             % Number of spectra (rows)
nCols = 500;           % Number of data points per spectrum (columns)

% Smoothing parameters
lambda = 50;         % Smoothing parameter (adjust for desired smoothness)
d = 1;                 % Order of the finite difference operator

% Noise level
noiseLevel = 0.1;     % Standard deviation of added Gaussian noise

% Seed for reproducibility
rng(42);

%% Generate Synthetic Spectral Data
% Create synthetic spectral peaks using Gaussian functions
x = linspace(0, 10, nCols);
inputMatrix = zeros(nRows, nCols);

for i = 1:nRows
    % Generate a spectrum with multiple Gaussian peaks
    inputMatrix(i, :) = ...
        exp(-((x - 3).^2) / (2 * 0.2^2)) + ... % Peak at x = 3
        0.5 * exp(-((x - 6).^2) / (2 * 0.5^2)) + ... % Broader peak at x = 6
        0.3 * exp(-((x - 8).^2) / (2 * 0.3^2)); % Narrow peak at x = 8
    
    % Add oscillatory baseline
    inputMatrix(i, :) = inputMatrix(i, :) + 0.1 * sin(2 * pi * x / 10);
end

% Add Gaussian noise to the spectra
noisyMatrix = inputMatrix + noiseLevel * randn(size(inputMatrix));

%% Apply Whittaker Smoother
% Call the WhittakerSmoother function
smoothedMatrix = whittakerSmoother(noisyMatrix, lambda, d);

%% Visualization
% Compare one example spectrum before and after smoothing
spectrumIndex = 1;  % Index of the spectrum to visualize

figure;

% Original spectrum
subplot(3, 1, 1);
plot(x, inputMatrix(spectrumIndex, :), 'LineWidth', 1.5);
title('Original Spectrum');
xlabel('Wavelength (arbitrary units)');
ylabel('Intensity');

% Noisy spectrum
subplot(3, 1, 2);
plot(x, noisyMatrix(spectrumIndex, :), 'LineWidth', 1.5);
title('Noisy Spectrum');
xlabel('Wavelength (arbitrary units)');
ylabel('Intensity');

% Smoothed spectrum
subplot(3, 1, 3);
plot(x, smoothedMatrix(spectrumIndex, :), 'LineWidth', 1.5);
title('Smoothed Spectrum');
xlabel('Wavelength (arbitrary units)');
ylabel('Intensity');

%% Statistical Evaluation
% Compute mean squared error (MSE) before and after smoothing
mseNoisy = mean((noisyMatrix(:) - inputMatrix(:)).^2);
mseSmoothed = mean((smoothedMatrix(:) - inputMatrix(:)).^2);

% Display MSE results
fprintf('Mean Squared Error (Noisy vs. Original): %.6f\n', mseNoisy);
fprintf('Mean Squared Error (Smoothed vs. Original): %.6f\n', mseSmoothed);

% Completion Message
disp('Whittaker smoothing applied successfully. Review the results above.');
