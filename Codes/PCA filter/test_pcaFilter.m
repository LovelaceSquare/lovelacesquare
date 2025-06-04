% test_pcaFilter.m
%
% This script generates synthetic Raman spectra as a linear combination
% of three underlying components with added noise. Each component
% contains unique Raman peaks. PCA filtering is applied to denoise the spectra,
% and a comparison is made between the original noisy and PCA-denoised spectra.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 2024-12-14
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0

%% Setup Parameters
nSamples = 150;           % Number of spectra (samples)
nPoints = 1000;          % Number of spectral points (e.g., Raman shift in cm^-1)
numComponents = 3;       % Number of principal components to retain for denoising
noiseLevel = 0.05;       % Standard deviation of Gaussian noise

% Define the x-axis (e.g., Raman shift)
x = linspace(0, 1000, nPoints);

%% Define Underlying Components
% Component 1: Peaks at 200, 500, 800 cm^-1
component1 = exp(-((x - 200).^2) / (2*30^2)) + ...
             exp(-((x - 500).^2) / (2*40^2)) + ...
             exp(-((x - 800).^2) / (2*25^2));

% Component 2: Peaks at 150, 600 cm^-1
component2 = exp(-((x - 150).^2) / (2*35^2)) + ...
             exp(-((x - 600).^2) / (2*50^2));

% Component 3: Peaks at 300, 700, 900 cm^-1
component3 = exp(-((x - 300).^2) / (2*20^2)) + ...
             exp(-((x - 700).^2) / (2*30^2)) + ...
             exp(-((x - 900).^2) / (2*40^2));

% Combine components into a matrix (S)
S = [component1; component2; component3];

%% Generate Linear Combinations of Components
% Randomly generate mixing coefficients (C) for each sample
C = rand(nSamples, 3); % Mixing matrix

% Generate dataset: D = C * S
data = C * S;

% Add Gaussian noise to the data
noisyData = data + noiseLevel * randn(size(data));

%% Apply PCA Filtering (Denoising)
filteredData = pcaFilter(noisyData, numComponents);

%% Visualization: Plot Every 10th Spectrum
sampleIndices = 1:50:nSamples; % Select every 10th spectrum
figure;
hold on;
for idx = sampleIndices
    % Plot original noisy spectrum in blue
    plot(x, noisyData(idx, :), 'b-', 'LineWidth', 1);
    % Plot corresponding PCA-denoised spectrum in red
    plot(x, filteredData(idx, :), 'r-', 'LineWidth', 1);
end
hold off;
xlabel('Raman Shift (cm^{-1})');
ylabel('Intensity (a.u.)');
title('Comparison: Original vs. PCA-Denoised Raman Spectra (Every 10th Sample)');
legend({'Original Spectrum', 'Denoised Spectrum'}, 'Location', 'best');
grid on;

%% Completion
disp('PCA filtering of synthetic Raman spectra completed. Review the plot for comparison.');
