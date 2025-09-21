% test_normMatrix.m
%
% This script generates synthetic Gaussian spectra with varying intensities,
% applies row-wise Euclidean normalization using normMatrix, and visualizes 
% the results before and after normalization for a clear comparison.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 2025-07-20
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v1.0

%% Generate Synthetic Gaussian Spectra
% Define channel axis
x = linspace(-5, 5, 200);

% Create three Gaussian spectra with different centers, widths, and amplitudes
spectrum1 = 1.0 * exp(-((x - 0.0).^2) / (2 * 0.5^2));   % Center 0.0, sigma=0.5, amp=1.0
spectrum2 = 3.0 * exp(-((x + 1.5).^2) / (2 * 1.0^2));   % Center -1.5, sigma=1.0, amp=3.0
spectrum3 = 0.5 * exp(-((x - 2.0).^2) / (2 * 0.3^2));   % Center 2.0, sigma=0.3, amp=0.5

% Combine into a matrix: each row is one spectrum
data = [spectrum1; spectrum2; spectrum3];

%% Normalize Spectra
% Perform row-wise Euclidean normalization
normalizedData = normMatrix(data, 'euclidean', 'row');

%% Plot Before and After Normalization
figure('Position',[100 100 900 350]);

% Original spectra
subplot(1,2,1);
plot(x, data, 'LineWidth',1.5);
title('Original Gaussian Spectra');
xlabel('Channel');
ylabel('Intensity');
legend('Spec 1','Spec 2','Spec 3','Location','best');
grid on;

% Normalized spectra
subplot(1,2,2);
plot(x, normalizedData, 'LineWidth',1.5);
title('Row-wise Euclidean Normalized Spectra');
xlabel('Channel');
ylabel('Normalized Intensity');
legend('Spec 1','Spec 2','Spec 3','Location','best');
grid on;

%% Completion Message
disp('normMatrix test completed: Gaussian spectra normalized and plotted.');
