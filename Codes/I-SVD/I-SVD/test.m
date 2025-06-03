% testI_SVD_Spectra.m
%
% This script generates synthetic spectral data resembling Raman spectra composed
% of three distinct components, adds low-level Gaussian noise, introduces a specified
% percentage of missing values, and evaluates the performance of the Iterative
% SVD-based PCA Imputation (I_SVD) algorithm by comparing the imputed scores and
% loadings against those obtained from Singular Value Decomposition (SVD) on
% complete data.
%
% Author:   Adrián Gómez-Sánchez
% Date:     2025-02-10
% License:  MIT
% Reviewed: Lovelace's Square
% Version:  v1.0
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
rng(123);  % For reproducibility
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

%% 5) Impute Incomplete Data using I_SVD
% -------------------------------------------------------------------------
% Use the I_SVD function to impute missing values in the incomplete dataset.

nComp = nComponents;    % Number of components
maxIter = 100;          % Maximum number of iterations
tol = 1e-10;            % Tolerance for convergence

% Perform iterative SVD-based PCA imputation
[Dimp_isvd, T_isvd, P_isvd, r2_isvd, lofc_isvd] = I_SVD(D_incomplete, nComp, maxIter, tol);

%% 6) Compare I_SVD Results with SVD Baseline
% -------------------------------------------------------------------------
% Align and compare the scores and loadings obtained from I_SVD with those from SVD.

% Initialize variables to store aligned I_SVD components
Scores_isvd_aligned = zeros(nSamples, nComponents);
Loadings_isvd_aligned = zeros(nComponents, nVariables);

% Keep track of which SVD components have been matched
matchedSvdComponents = false(1, nComponents);

% Compute the correlation matrix between SVD loadings and I_SVD loadings
corrMatrix = corr(Loadings_svd', P_isvd', 'Rows', 'complete');  % Size: (nComponents x nComponents)

for isvdComp = 1:nComponents
    % Find the SVD component with the highest absolute correlation to the current I_SVD component
    [maxCorrVal, bestSvdComp] = max(abs(corrMatrix(:, isvdComp)));

    % Check if this SVD component has already been matched
    if matchedSvdComponents(bestSvdComp)
        % Find the next best match
        [sortedCorr, sortedIdx] = sort(abs(corrMatrix(:, isvdComp)), 'descend');
        for idx = 1:length(sortedCorr)
            candidateComp = sortedIdx(idx);
            if ~matchedSvdComponents(candidateComp)
                bestSvdComp = candidateComp;
                maxCorrVal = sortedCorr(idx);
                break;
            end
        end
        if matchedSvdComponents(bestSvdComp)
            error('Unable to match ISVD components to SVD components.');
        end
    end

    % Assign the I_SVD component to the matched SVD component
    matchedSvdComponents(bestSvdComp) = true;
    Loadings_isvd_aligned(bestSvdComp, :) = P_isvd(isvdComp, :);
    Scores_isvd_aligned(:, bestSvdComp) = T_isvd(:, isvdComp);

    % Flip sign if the actual correlation is negative
    actualCorr = corr(Loadings_svd(bestSvdComp,:)', P_isvd(isvdComp,:)','Rows','complete');
    if actualCorr < 0
        Loadings_isvd_aligned(bestSvdComp, :) = -Loadings_isvd_aligned(bestSvdComp, :);
        Scores_isvd_aligned(:, bestSvdComp) = -Scores_isvd_aligned(:, bestSvdComp);
    end

    % Compute scaling factor to best match SVD scores
    scale = (Scores_isvd_aligned(:, bestSvdComp)' * Scores_svd(:, bestSvdComp)) / ...
            (Scores_isvd_aligned(:, bestSvdComp)' * Scores_isvd_aligned(:, bestSvdComp));

    % Apply scaling to I_SVD scores and loadings
    Scores_isvd_aligned(:, bestSvdComp) = Scores_isvd_aligned(:, bestSvdComp) * scale;
    Loadings_isvd_aligned(bestSvdComp, :) = Loadings_isvd_aligned(bestSvdComp, :) / scale;
end

%% 7) Calculate Similarity Metrics for Scores and Loadings
% -------------------------------------------------------------------------
% Compute correlation coefficients and Mean Squared Errors (MSE) between
% the aligned I_SVD and SVD scores and loadings.

% Initialize metrics
corrScores = zeros(1, nComponents);
corrLoadings = zeros(1, nComponents);
mseScores = zeros(1, nComponents);
mseLoadings = zeros(1, nComponents);

for compIdx = 1:nComponents
    % Correlation for Scores
    corrScores(compIdx) = corr(Scores_svd(:, compIdx), Scores_isvd_aligned(:, compIdx));

    % Correlation for Loadings
    corrLoadings(compIdx) = corr(Loadings_svd(compIdx, :)', Loadings_isvd_aligned(compIdx, :)');

    % MSE for Scores
    mseScores(compIdx) = mean( (Scores_svd(:, compIdx) - Scores_isvd_aligned(:, compIdx)).^2 );

    % MSE for Loadings
    mseLoadings(compIdx) = mean( (Loadings_svd(compIdx, :)' - Loadings_isvd_aligned(compIdx, :)').^2 );
end

%% 8) Visualization of Scores and Loadings Comparison
% -------------------------------------------------------------------------
% Create subplots to compare scores and loadings from SVD and I_SVD.

figure('Name','SVD vs ISVD Scores and Loadings Comparison','Color','w');

% --- Subplot 1: Scores Comparison ---
subplot(2,1,1); hold on;
for compIdx = 1:nComponents
    plot(Scores_svd(:, compIdx), '-o','LineWidth',1.5, 'MarkerSize',6, 'DisplayName', sprintf('SVD Comp %d', compIdx));
    plot(Scores_isvd_aligned(:, compIdx), '--x','LineWidth',1.5, 'MarkerSize',6, 'DisplayName', sprintf('ISVD Comp %d', compIdx));
end
hold off;
xlabel('Sample Index');
ylabel('Score Value');
title('Scores Comparison: SVD vs I\_SVD');
legend('Location','best');
grid on;

% --- Subplot 2: Loadings Comparison ---
subplot(2,1,2); hold on;
for compIdx = 1:nComponents
    plot(wavenumber, Loadings_svd(compIdx,:), '-','LineWidth',1.5, 'DisplayName', sprintf('SVD Comp %d', compIdx));
    plot(wavenumber, Loadings_isvd_aligned(compIdx,:), '--','LineWidth',1.5, 'DisplayName', sprintf('ISVD Comp %d', compIdx));
end
hold off;
xlabel('Wavenumber (cm^{-1})');
ylabel('Loading Intensity');
title('Loadings Comparison: SVD vs I\_SVD');
legend('Location','best');
grid on;

sgtitle('Comparison of SVD (Complete Data) vs I\_SVD (Imputed Data)');

%% 9) Display Similarity Metrics
% -------------------------------------------------------------------------
% Show correlation coefficients and MSE values for scores and loadings.

fprintf('Similarity Metrics between SVD and ISVD Results:\n\n');

for compIdx = 1:nComponents
    fprintf('Component %d:\n', compIdx);
    fprintf(' - Score Correlation: %.4f\n', corrScores(compIdx));
    fprintf(' - Score MSE        : %.6e\n', mseScores(compIdx));
    fprintf(' - Loading Correlation: %.4f\n', corrLoadings(compIdx));
    fprintf(' - Loading MSE        : %.6e\n\n', mseLoadings(compIdx));
end

%% 10) Plot Convergence Metrics from I_SVD
% -------------------------------------------------------------------------
% Display the convergence of explained variance and lack of fit over iterations.

figure('Name','ISVD Convergence Metrics','Color','w');

% Plot Explained Variance (r2)
subplot(2,1,1); hold on;
plot(1:length(r2_isvd), r2_isvd, '-bo', 'LineWidth',1.5, 'MarkerSize',6);
xlabel('Iteration');
ylabel('Explained Variance (R^2)');
title('Convergence of Explained Variance');
grid on; hold off;

% Plot Lack of Fit (lofc)
subplot(2,1,2); hold on;
plot(1:length(lofc_isvd), lofc_isvd, '-ro', 'LineWidth',1.5, 'MarkerSize',6);
xlabel('Iteration');
ylabel('Lack of Fit (LOF)');
title('Convergence of Lack of Fit');
grid on; hold off;

sgtitle('I\_SVD Convergence Metrics');

%% 11) Display Final Metrics
% -------------------------------------------------------------------------
% Show final lack of fit and MSE values.

fprintf('Final Metrics for ISVD Imputation:\n');
fprintf(' - Final Lack of Fit (LOF) : %.6f\n', lofc_isvd(end));
fprintf(' - Mean Squared Error (MSE) : %.6e\n\n', mean(mseScores));

disp('Test completed. Review figures for Scores and Loadings comparison and Convergence metrics.');
