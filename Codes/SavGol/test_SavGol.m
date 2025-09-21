% test_SavGol.m
%
% This script generates synthetic spectral data with noise, applies the
% Savitzky-Golay filter, and evaluates the results visually and statistically.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 2024-12-14
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0

%% Setup Parameters
% Dimensions of the synthetic data
nSamples = 5;           % Number of spectra (rows)
nPoints = 500;          % Number of data points per spectrum (columns)

% Savitzky-Golay parameters
windowSize = 31;        % Must be odd and > polyOrder
polyOrder = 1;          % Polynomial order (e.g., linear)
derivOrder = 0;         % Derivative order (0 = smoothing)
edgeMethod = 'Extrapolation'; % Edge handling strategy: 'None', 'Reflection', etc.

% Noise level
noiseLevel = 0.1;      % Standard deviation of Gaussian noise

% Random seed for reproducibility
rng(42);

%% Generate Synthetic Spectral Data
% Create synthetic spectral peaks with Gaussian profiles
x = linspace(0, 10, nPoints);
data = zeros(nSamples, nPoints);

for i = 1:nSamples
    % Generate a spectrum with multiple Gaussian peaks
    data(i, :) = ...
        exp(-((x - 3).^2) / (2 * 0.3^2)) + ... % Peak at x = 3
        0.8 * exp(-((x - 6).^2) / (2 * 0.4^2)) + ... % Broader peak at x = 6
        0.5 * exp(-((x - 8).^2) / (2 * 0.2^2)); % Narrow peak at x = 8
    
    % Add oscillatory baseline
    data(i, :) = data(i, :) + 0.05 * sin(2 * pi * x / 2);
end

% Add Gaussian noise
noisyData = data + noiseLevel * randn(size(data));

%% Apply Savitzky-Golay Filter
% Call the SavGol function
filteredData = SavGol(noisyData, windowSize, polyOrder, derivOrder, edgeMethod);

%% Visualization
% Compare one example spectrum before and after filtering
spectrumIndex = 1;  % Index of the spectrum to visualize

figure;

% Original spectrum
subplot(3, 1, 1);
plot(x, data(spectrumIndex, :), 'LineWidth', 1.5);
title('Original Spectrum');
xlabel('Wavelength (arbitrary units)');
ylabel('Intensity');

% Noisy spectrum
subplot(3, 1, 2);
plot(x, noisyData(spectrumIndex, :), 'LineWidth', 1.5);
title('Noisy Spectrum');
xlabel('Wavelength (arbitrary units)');
ylabel('Intensity');

% Filtered spectrum
subplot(3, 1, 3);
plot(x, filteredData(spectrumIndex, :), 'LineWidth', 1.5);
title('Filtered Spectrum');
xlabel('Wavelength (arbitrary units)');
ylabel('Intensity');

%% Statistical Evaluation
% Compute mean squared error (MSE) before and after filtering
mseNoisy = mean((noisyData(:) - data(:)).^2);
mseFiltered = mean((filteredData(:) - data(:)).^2);

% Display MSE results
fprintf('Mean Squared Error (Noisy vs. Original): %.6f\n', mseNoisy);
fprintf('Mean Squared Error (Filtered vs. Original): %.6f\n', mseFiltered);

% Completion Message
disp('Savitzky-Golay filtering applied successfully. Review the results above.');
