% =========================================================================
% test_NIPALS.m
% =========================================================================
%
% Test script for NIPALS (Nonlinear Iterative Partial Least Squares)
% Compares NIPALS decomposition on incomplete data (with NaNs) against
% SVD on complete data as the reference solution.
%
% Author:   Adrián Gómez-Sánchez
% Date:     2025-01-02
% License:  MIT
% Reviewed: Lovelace's Square
% Version:  v 1.0
%
% -------------------------------------------------------------------------
% REFERENCES:
%   • Wold, H. (1966).
%     "Nonlinear estimation by iterative least square procedures."
%     Research Papers in Statistics: Festschrift for J. Neyman, 411-444.
%
%   • Wold, S., Esbensen, K., & Geladi, P. (1987).
%     "Principal component analysis." Chemometrics and intelligent
%     laboratory systems, 2(1-3), 37-52.
%
% =========================================================================
clear; clc; close all;

%% ========================================================================
%% 1. Generate Synthetic Bilinear Spectral Data
%% ========================================================================
fprintf('Generating synthetic Gaussian spectral data...\n\n');

% Reproducibility
rng(42, 'twister');

% Data dimensions
nSamples     = 60;    % Number of samples (e.g., time points)
nVariables   = 150;   % Number of variables (e.g., wavelengths)
nComponents  = 3;     % Number of components

% Variable axis (e.g., wavelength range)
varIndex = linspace(400, 700, nVariables);  % 400-700 nm

% --- Define true loadings with distinct Gaussian peaks ---
loadingsTrue = zeros(nComponents, nVariables);
loadingsTrue(1,:) = exp(-((varIndex - 470)/50).^2);   % Component 1: 470 nm
loadingsTrue(2,:) = exp(-((varIndex - 550)/40).^2);   % Component 2: 550 nm
loadingsTrue(3,:) = exp(-((varIndex - 630)/50).^2);   % Component 3: 630 nm

% --- Generate random scores for each sample ---
scoresTrue = rand(nSamples, nComponents);

% --- Create noiseless bilinear data ---
D_noiseless = scoresTrue * loadingsTrue;

%% ========================================================================
%% 2. Add Low-Level Gaussian Noise
%% ========================================================================
noiseLevel = 0.02;  % 2% noise
noiseMatrix = noiseLevel * randn(nSamples, nVariables);
D_complete  = D_noiseless + noiseMatrix;

%% ========================================================================
%% 3. Introduce Missing Values
%% ========================================================================
missingPercent = 10;  % 10% missing data
missingMask    = rand(nSamples, nVariables) < (missingPercent / 100);
D_incomplete   = D_complete;
D_incomplete(missingMask) = NaN;

actualMissing = 100 * sum(isnan(D_incomplete(:))) / numel(D_incomplete);

fprintf('Synthetic spectral data generated.\n');
fprintf(' - Number of Samples    : %d\n', nSamples);
fprintf(' - Number of Variables  : %d\n', nVariables);
fprintf(' - Number of Components : %d\n', nComponents);
fprintf(' - Noise Level          : %.2f%%\n', noiseLevel*100);
fprintf(' - Missing Data         : %.2f%%\n\n', actualMissing);

%% ========================================================================
%% 4. Analyze Complete Data using SVD
%% ========================================================================
fprintf('Running SVD on complete data (reference solution)...\n');

% Perform SVD
[U_svd, S_svd, V_svd] = svd(D_complete, 'econ');

% Extract scores and loadings for the desired number of components
U_svd_trunc = U_svd(:, 1:nComponents);
S_svd_trunc = S_svd(1:nComponents, 1:nComponents);
V_svd_trunc = V_svd(:, 1:nComponents);

Scores_svd   = U_svd_trunc * S_svd_trunc;    % (nSamples x nComponents)
Loadings_svd = V_svd_trunc';                 % (nComponents x nVariables)

% Reconstruct the data from SVD
D_svd_reconstructed = Scores_svd * Loadings_svd;

fprintf('SVD completed.\n\n');

%% ========================================================================
%% 5. Analyze Incomplete Data using NIPALS
%% ========================================================================
fprintf('Running NIPALS on incomplete data (with missing values)...\n');
fprintf('========================================\n\n');

maxIter = 100;

% Ensure function is available
assert(exist('NIPALS','file')==2, 'NIPALS.m not found on path.');

[Dr_nipals, T_nipals, P_nipals] = NIPALS(D_incomplete, maxIter, nComponents);

% Calculate quality metrics
[r2_nipals, lof_nipals] = lofNaN(D_incomplete, T_nipals, P_nipals);

fprintf('NIPALS completed.\n');
fprintf('  Final R²  : %.2f%%\n', r2_nipals);
fprintf('  Final LOF : %.4f%%\n\n', lof_nipals);

%% ========================================================================
%% 6. Align NIPALS Components to SVD Reference
%% ========================================================================
%% 6. Align NIPALS Component Signs to SVD Reference (no reordering)
%% ========================================================================
fprintf('Aligning NIPALS component signs to SVD reference (no reordering)...\n\n');

% Initialize "aligned" variables as direct copies
Scores_nipals_aligned   = T_nipals;
Loadings_nipals_aligned = P_nipals;

for comp = 1:nComponents
    % Correlation between corresponding SVD and NIPALS loadings
    corrVal = corr(Loadings_svd(comp, :)', P_nipals(comp, :)', 'Rows', 'complete');

    % If correlation is negative, flip sign of both scores and loadings
    if corrVal < 0
        Scores_nipals_aligned(:, comp)   = -Scores_nipals_aligned(:, comp);
        Loadings_nipals_aligned(comp, :) = -Loadings_nipals_aligned(comp, :);
        fprintf('  Component %d sign flipped (corr = %.4f)\n', comp, corrVal);
    else
        fprintf('  Component %d sign OK (corr = %.4f)\n', comp, corrVal);
    end
end

fprintf('\n');


%% ========================================================================
%% 7. Calculate Mean Squared Errors (MSE)
%% ========================================================================
validMask = ~missingMask;  % Logical mask for non-missing data

% Calculate MSE for NIPALS (on non-missing values)
mseNipals = mean((D_complete(validMask) - Dr_nipals(validMask)).^2);

% Calculate MSE for SVD
mseSvd = mean((D_complete(validMask) - D_svd_reconstructed(validMask)).^2);

fprintf('Mean Squared Error (MSE) on non-missing data:\n');
fprintf(' - NIPALS Reconstruction vs Original: %.6e\n', mseNipals);
fprintf(' - SVD Reconstruction vs Original:    %.6e\n\n', mseSvd);

%% ========================================================================
%% 8. Visualization of Scores and Loadings (SVD vs NIPALS)
%% ========================================================================
fprintf('========================================\n');
fprintf('Generating comparison plots...\n\n');

figure('Name', 'SVD vs NIPALS Comparison', 'Position', [100, 100, 1400, 800]);

% --- Plot Scores Comparison ---
for i = 1:nComponents
    subplot(2, nComponents, i);
    hold on;
    plot(Scores_svd(:,i), '-o', 'LineWidth', 2, 'MarkerSize', 4, ...
        'DisplayName', sprintf('SVD Comp %d', i));
    plot(Scores_nipals_aligned(:,i), '--x', 'LineWidth', 2, 'MarkerSize', 4, ...
        'DisplayName', sprintf('NIPALS Comp %d', i));
    hold off;
    xlabel('Sample Index', 'FontSize', 10);
    ylabel('Score Value', 'FontSize', 10);
    title(sprintf('Scores - Component %d', i), 'FontSize', 11, 'FontWeight', 'bold');
    legend('Location', 'best');
    grid on;
end

% --- Plot Loadings Comparison ---
for i = 1:nComponents
    subplot(2, nComponents, nComponents + i);
    hold on;
    plot(varIndex, Loadings_svd(i,:), '-', 'LineWidth', 2, ...
        'DisplayName', sprintf('SVD Comp %d', i));
    plot(varIndex, Loadings_nipals_aligned(i,:), '--', 'LineWidth', 2, ...
        'DisplayName', sprintf('NIPALS Comp %d', i));
    hold off;
    xlabel('Wavelength (nm)', 'FontSize', 10);
    ylabel('Loading Value', 'FontSize', 10);
    title(sprintf('Loadings - Component %d', i), 'FontSize', 11, 'FontWeight', 'bold');
    legend('Location', 'best');
    grid on;
end

sgtitle('SVD (Complete Data) vs NIPALS (Incomplete Data) Comparison', ...
    'FontSize', 14, 'FontWeight', 'bold');

%% ========================================================================
%% 9. Reconstruction Quality Visualization
%% ========================================================================
figure('Name', 'NIPALS Reconstruction Quality', 'Position', [150, 150, 1400, 420]);

subplot(1,4,1);
imagesc(varIndex, 1:nSamples, D_complete); colorbar;
title('Complete Data (with noise)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Wavelength (nm)', 'FontSize', 10);
ylabel('Sample', 'FontSize', 10);

subplot(1,4,2);
imagesc(varIndex, 1:nSamples, D_incomplete); colorbar;
title('Incomplete Data (with NaNs)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Wavelength (nm)', 'FontSize', 10);
ylabel('Sample', 'FontSize', 10);

subplot(1,4,3);
imagesc(varIndex, 1:nSamples, Dr_nipals); colorbar;
title('NIPALS Reconstruction', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Wavelength (nm)', 'FontSize', 10);
ylabel('Sample', 'FontSize', 10);

subplot(1,4,4);
residual_nipals = D_complete - Dr_nipals;
imagesc(varIndex, 1:nSamples, residual_nipals); colorbar;
title('Residual (Complete - NIPALS)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Wavelength (nm)', 'FontSize', 10);
ylabel('Sample', 'FontSize', 10);

sgtitle(sprintf('NIPALS Reconstruction (R² = %.2f%%, LOF = %.4f%%)', r2_nipals, lof_nipals), ...
    'FontSize', 14, 'FontWeight', 'bold');

%% ========================================================================
%% 10. Missing Value Imputation Quality
%% ========================================================================
fprintf('========================================\n');
fprintf('Evaluating missing value imputation quality...\n\n');

% Extract imputed values for missing positions
imputedValues = Dr_nipals(missingMask);
trueValues = D_complete(missingMask);

% Calculate imputation error
imputError_rmse = sqrt(mean((trueValues - imputedValues).^2));
imputError_r2 = 100 * (1 - sum((trueValues - imputedValues).^2) / ...
    sum((trueValues - mean(trueValues)).^2));

figure('Name', 'Missing Value Imputation Quality', 'Position', [200, 200, 1000, 400]);

subplot(1,2,1);
scatter(trueValues, imputedValues, 30, 'filled', 'MarkerFaceAlpha', 0.6);
hold on;
plot([min(trueValues), max(trueValues)], [min(trueValues), max(trueValues)], ...
    'r--', 'LineWidth', 2);
hold off;
xlabel('True Values', 'FontSize', 11);
ylabel('Imputed Values', 'FontSize', 11);
title('Imputed vs True Values', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
axis equal tight;

subplot(1,2,2);
errorHist = trueValues - imputedValues;
histogram(errorHist, 30, 'FaceColor', [0.3 0.7 0.9], 'EdgeColor', 'k');
xlabel('Imputation Error', 'FontSize', 11);
ylabel('Frequency', 'FontSize', 11);
title(sprintf('Error Distribution (RMSE=%.4f, R²=%.2f%%)', imputError_rmse, imputError_r2), ...
    'FontSize', 12, 'FontWeight', 'bold');
grid on;

%% ========================================================================
%% Final Summary
%% ========================================================================
fprintf('========================================\n');
fprintf('Test completed successfully!\n');
fprintf('========================================\n');
fprintf('Data Quality:\n');
fprintf('  Noise level           : %.2f%%\n', noiseLevel*100);
fprintf('  Missing data          : %.2f%%\n\n', actualMissing);
fprintf('NIPALS Performance:\n');
fprintf('  Final R²              : %.2f%%\n', r2_nipals);
fprintf('  Final LOF             : %.4f%%\n', lof_nipals);
fprintf('  Reconstruction MSE    : %.6e\n', mseNipals);
fprintf('  Imputation RMSE       : %.4f\n', imputError_rmse);
fprintf('  Imputation R²         : %.2f%%\n\n', imputError_r2);
fprintf('Reference (SVD):\n');
fprintf('  Reconstruction MSE    : %.6e\n', mseSvd);
fprintf('========================================\n');
