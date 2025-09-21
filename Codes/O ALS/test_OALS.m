% testOALS_Spectra.m
%
% This script generates synthetic spectral data resembling Raman spectra composed
% of three distinct components, adds low-level Gaussian noise, introduces a specified
% percentage of missing values, and compares the decomposition results using
% Singular Value Decomposition (SVD) on complete data and Orthogonalized Alternating
% Least Squares (OALS) on incomplete data.
%
% Author:   Adrián Gómez-Sánchez
% Date:     2025-02-01
% License:  MIT
% Reviewed: Lovelace's Square
% Version:  v 1.0
%

clear; clc; close all;

%% 1) Generate Synthetic Bilinear Spectral Data
% -------------------------------------------------------------------------
% Simulate 3 underlying spectral components (loadings) with distinct peak shapes
% and combine them with random scores to form bilinear data.

nSamples     = 40;     % Number of samples
nVariables   = 500;    % Number of spectral points (e.g., wavenumbers)
nComponents  = 3;      % Number of true components

% --- Wavenumber axis for spectral data ---
wavenumber = linspace(500, 3500, nVariables);  % e.g., 500 to 3500 cm^{-1}

% --- Define true loadings with distinct Gaussian peaks ---
loadingsTrue = zeros(nComponents, nVariables);

% Component 1: Broad Gaussian peak around 1000 cm^{-1}
peak1_center = 1000;
peak1_width  = 80;
loadingsTrue(1,:) = exp(-0.5 * ((wavenumber - peak1_center) / peak1_width).^2);

% Component 2: Narrow Gaussian peak around 2000 cm^{-1}
peak2_center = 2000;
peak2_width  = 40;
loadingsTrue(2,:) = exp(-0.5 * ((wavenumber - peak2_center) / peak2_width).^2);

% Component 3: Medium Gaussian peak around 3000 cm^{-1}
peak3_center = 3000;
peak3_width  = 60;
loadingsTrue(3,:) = exp(-0.5 * ((wavenumber - peak3_center) / peak3_width).^2);

% --- Generate random scores for each sample ---
rng(456);  % For reproducibility
scoresTrue = rand(nSamples, nComponents);

% --- Create noiseless bilinear data ---
D_noiseless = scoresTrue * loadingsTrue;

%% 2) Add Low-Level Gaussian Noise
% -------------------------------------------------------------------------
% Introduce Gaussian noise to simulate measurement imperfections.

noiseLevel = 0.01;  % 1% noise
noiseMatrix = noiseLevel * randn(nSamples, nVariables);
D_complete  = D_noiseless + noiseMatrix;

%% 3) Introduce Missing Values
% -------------------------------------------------------------------------
% Randomly assign a specified percentage of data points as missing (NaN).

missingPercent = 50;  % 50% missing data
missingMask    = rand(nSamples, nVariables) < (missingPercent / 100);
D_incomplete   = D_complete;
D_incomplete(missingMask) = NaN;

fprintf('Synthetic spectral data generated.\n');
fprintf(' - Number of Samples    : %d\n', nSamples);
fprintf(' - Number of Variables  : %d\n', nVariables);
fprintf(' - Number of Components : %d\n', nComponents);
fprintf(' - Noise Level          : %.2f%%\n', noiseLevel*100);
fprintf(' - Missing Data         : %.2f%%\n\n', missingPercent);

%% 4) Analyze Complete Data using SVD
% -------------------------------------------------------------------------
% Perform Singular Value Decomposition on the complete (no missing data)
% dataset to establish a baseline for comparison.

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

%% 5) Analyze Incomplete Data using OALS
% -------------------------------------------------------------------------
% Use the OALS function to decompose the incomplete dataset, handling missing values.

iter = 200;  % Maximum number of iterations

% OALS requires the number of components as a mandatory argument
[Dr_oals, T_oals, P_oals, r2_oals, lofc_oals] = OALS(D_incomplete, iter, nComponents);

%% 6) Compare SVD and OALS Results
% -------------------------------------------------------------------------
% Align signs and order of OALS components to match those from SVD based on maximum correlation.

% Initialize variables to store aligned OALS components
Scores_oals_aligned = zeros(nSamples, nComponents);
Loadings_oals_aligned = zeros(nComponents, nVariables);

% Keep track of which SVD components have been matched
matchedSvdComponents = false(1, nComponents);

% Compute the full correlation matrix between SVD loadings and OALS loadings
corrMatrix = corr(Loadings_svd', P_oals', 'Rows', 'complete');  % Size: (nComponents x nComponents)

for oalsComp = 1:nComponents
    % Find the SVD component with the highest absolute correlation to the current OALS component
    [maxCorrVal, bestSvdComp] = max(abs(corrMatrix(:, oalsComp)));
    
    % Check if this SVD component has already been matched
    if matchedSvdComponents(bestSvdComp)
        % Find the next best match
        [maxCorrVal, bestSvdComp] = max(abs(corrMatrix(:, oalsComp)));
        if matchedSvdComponents(bestSvdComp)
            error('Unable to match OALS components to SVD components.');
        end
    end
    
    % Assign the OALS component to the matched SVD component
    matchedSvdComponents(bestSvdComp) = true;
    Loadings_oals_aligned(bestSvdComp, :) = P_oals(oalsComp, :);
    Scores_oals_aligned(:, bestSvdComp) = T_oals(:, oalsComp);
    
    % Flip sign if the actual correlation is negative
    actualCorr = corr(Loadings_svd(bestSvdComp,:)', P_oals(oalsComp,:)','Rows','complete');
    if actualCorr < 0
        Loadings_oals_aligned(bestSvdComp, :) = -Loadings_oals_aligned(bestSvdComp, :);
        Scores_oals_aligned(:, bestSvdComp) = -Scores_oals_aligned(:, bestSvdComp);
    end
    
    % Compute scaling factor to best match SVD scores
    scale = (Scores_oals_aligned(:, bestSvdComp)' * Scores_svd(:, bestSvdComp)) / ...
            (Scores_oals_aligned(:, bestSvdComp)' * Scores_oals_aligned(:, bestSvdComp));
    
    % Apply scaling to OALS scores and loadings
    Scores_oals_aligned(:, bestSvdComp) = Scores_oals_aligned(:, bestSvdComp) * scale;
    Loadings_oals_aligned(bestSvdComp, :) = Loadings_oals_aligned(bestSvdComp, :) / scale;
end

%% 7) Calculate Mean Squared Errors (MSE)
% -------------------------------------------------------------------------
% Compute MSE between reconstructed data and original complete data only on non-missing entries.

validMask = ~missingMask;  % Logical mask for non-missing data

% Calculate MSE for OALS
mseOals = mean( (D_complete(validMask) - Dr_oals(validMask)).^2 );

% Calculate MSE for SVD
mseSvd = mean( (D_complete(validMask) - D_svd_reconstructed(validMask)).^2 );

fprintf('Mean Squared Error (MSE):\n');
fprintf(' - OALS Reconstruction vs Original (non-missing): %.6e\n', mseOals);
fprintf(' - SVD Reconstruction vs Original (non-missing):  %.6e\n\n', mseSvd);

%% 8) Visualization of Scores and Loadings
% -------------------------------------------------------------------------
% Create subplots to compare scores and loadings from SVD and OALS.

figure('Name','SVD vs OALS Comparison','Color','w');

% --- Subplot 1: Scores Comparison ---
subplot(1,2,1); hold on;
% Plot SVD scores
for compIdx = 1:nComponents
    plot(Scores_svd(:, compIdx), '-o','LineWidth',1.5, 'MarkerSize',6, 'DisplayName', sprintf('SVD Comp %d', compIdx));
end
% Plot OALS scores (aligned)
for compIdx = 1:nComponents
    plot(Scores_oals_aligned(:, compIdx), '--x','LineWidth',1.5, 'MarkerSize',6, 'DisplayName', sprintf('OALS Comp %d', compIdx));
end
hold off;
xlabel('Sample Index');
ylabel('Score Value');
title('Scores Comparison');
legend('Location','best');
grid on;

% --- Subplot 2: Loadings Comparison ---
subplot(1,2,2); hold on;
% Plot SVD loadings
for compIdx = 1:nComponents
    plot(wavenumber, Loadings_svd(compIdx,:), '-','LineWidth',1.5, 'DisplayName', sprintf('SVD Comp %d', compIdx));
end
% Plot OALS loadings (aligned)
for compIdx = 1:nComponents
    plot(wavenumber, Loadings_oals_aligned(compIdx,:), '--','LineWidth',1.5, 'DisplayName', sprintf('OALS Comp %d', compIdx));
end
hold off;
xlabel('Wavenumber (cm^{-1})');
ylabel('Loading Intensity');
title('Loadings Comparison');
legend('Location','best');
grid on;

sgtitle('Comparison of SVD (Complete Data) vs OALS (Missing Data)');

%% 9) Display Final Metrics
% -------------------------------------------------------------------------
% Show final lack of fit and MSE values.

fprintf('Final Metrics:\n');
fprintf(' - OALS Final Lack of Fit (LOF): %.6f\n', lofc_oals(end));
fprintf(' - OALS Mean Squared Error (MSE)  : %.6e\n', mseOals);
fprintf(' - SVD Mean Squared Error (MSE)   : %.6e\n\n', mseSvd);

disp('Test completed. Review figures for Scores and Loadings comparison.');
