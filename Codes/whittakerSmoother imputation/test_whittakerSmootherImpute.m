% test_whittakerSmootherImpute.m
%
% This script generates synthetic spectral data with noise and missing values,
% applies the Whittaker smoother with imputation, and evaluates the results.

% Author: Adrián Gómez-Sánchez
% Date Created: 2024-12-14
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0

%% Setup Parameters
% Define the dimensions of the synthetic data
nRows = 50;         % Number of spectra (rows)
nCols = 500;       % Number of data points per spectrum (columns)

% Smoothing parameters
lambda = 500;       % Smoothing parameter (adjust for desired smoothness)
d = 2;             % Order of the finite difference operator (common = 2)

% Noise level
noiseLevel = 0.05;  % Standard deviation of added Gaussian noise

% Iteration parameters
maxIter = 150;      % Maximum iterations for iterative imputation
tol = 1e-10;        % Convergence tolerance


%% Generate Synthetic Spectral Data
x = linspace(0, 10, nCols);
inputMatrix = zeros(nRows, nCols);

for i = 1:nRows
    % Generate a spectrum with multiple Gaussian peaks
    inputMatrix(i, :) = ...
        exp(-((x - 3).^2) / (2 * 0.2^2)) + ...       % Peak at x = 3
        0.5 * exp(-((x - 6).^2) / (2 * 0.5^2)) + ... % Broader peak at x = 6
        0.3 * exp(-((x - 8).^2) / (2 * 0.3^2));      % Narrow peak at x = 8

    % Add oscillatory baseline
    inputMatrix(i, :) = inputMatrix(i, :) + 0.1 * sin(2 * pi * x / 10);
end

% Add Gaussian noise
noisyMatrix = inputMatrix + noiseLevel * randn(size(inputMatrix));

%% Introduce Missing Data (NaNs)
missingFraction = 0.1; % ~10% missing
totalPoints = numel(noisyMatrix);
numMissing = round(missingFraction * totalPoints);
missingIndices = randperm(totalPoints, numMissing);
noisyMatrix(missingIndices) = NaN;

%% Apply Whittaker Smoother with Imputation
[smoothedMatrix, finalImputedMatrix] = whittakerSmootherImpute(noisyMatrix, ...
    lambda, d, maxIter, tol);

%% Visualization
% Compare one example spectrum before and after
spectrumIndex = 1;

figure('Name','Whittaker Smoother with Imputation','NumberTitle','off');

% Original spectrum
subplot(3,1,1);
plot(x, inputMatrix(spectrumIndex, :), 'LineWidth',1.5);
title('Original Spectrum');
xlabel('Wavelength');
ylabel('Intensity');

% Noisy spectrum with NaNs
subplot(3,1,2);
plot(x, noisyMatrix(spectrumIndex, :), 'LineWidth',1.5);
title('Noisy Spectrum (with NaNs)');
xlabel('Wavelength');
ylabel('Intensity');

% Final smoothed (imputed) spectrum
subplot(3,1,3);
plot(x, smoothedMatrix(spectrumIndex, :), 'LineWidth',1.5);
hold on;
% Mark the originally missing points
nanPositions = isnan(noisyMatrix(spectrumIndex, :));
plot(x(nanPositions), smoothedMatrix(spectrumIndex, nanPositions), 'ro', ...
    'MarkerSize',5, 'DisplayName','Imputed Points');
hold off;
title('Whittaker-Smoothed (Imputed) Spectrum');
xlabel('Wavelength');
ylabel('Intensity');
legend('Smoothed','Imputed Values','Location','best');

%% Statistical Evaluation
% Compare only valid points in the original input
validMask = ~isnan(noisyMatrix);
mseNoisy = mean((noisyMatrix(validMask) - inputMatrix(validMask)).^2);
mseSmoothed = mean((smoothedMatrix(validMask) - inputMatrix(validMask)).^2);

fprintf('Mean Squared Error (Noisy vs. Original):    %.6f\n', mseNoisy);
fprintf('Mean Squared Error (Smoothed vs. Original): %.6f\n', mseSmoothed);

disp('Whittaker smoothing with imputation completed successfully.');
