% =========================================================================
% test_MCR_ALS_Lite.m
% =========================================================================
%
% Test script for MCR_ALS_Lite (Multivariate Curve Resolution - Alternating Least Squares)
% Demonstrates the algorithm on synthetic spectral data generated from
% Gaussian-shaped concentration and spectral profiles.
%
% Author:   Adrián Gómez-Sánchez
% Date:     2025-10-30
% License:  MIT
% Reviewed: Lovelace's Square
% Version:  v 1.0
%
% -------------------------------------------------------------------------
% REFERENCES:
%   • Lawton, W. H., & Sylvestre, E. A. (1971).
%     "Self modeling curve resolution." Technometrics, 13(3), 617–633.
%
%   • de Juan, A., & Tauler, R. (2021).
%     "Multivariate Curve Resolution: 50 years addressing the mixture
%      analysis problem – A review."
%      Analytica Chimica Acta, 1145, 59–78. Elsevier.
%
% =========================================================================
clear; clc; close all;

%% ========================================================================
%% 1. Generate Synthetic Gaussian Data
%% ========================================================================
fprintf('Generating synthetic Gaussian spectral data...\n\n');

% Reproducibility
rng(42, 'twister');

% Data dimensions
nSamples     = 50;    % Number of samples (e.g., time points)
nVariables   = 200;   % Number of variables (e.g., wavelengths)
nComponents  = 3;     % Number of chemical components

sampleIndex  = linspace(0, 1, nSamples);
varIndex     = linspace(400, 700, nVariables);  % e.g., wavelength range (nm)

% --- Create Gaussian concentration profiles (C_true) ---
C_true = zeros(nSamples, nComponents);
C_true(:,1) = exp(-((sampleIndex - 0.2)/0.2).^2);  % Peak early
C_true(:,2) = exp(-((sampleIndex - 0.5)/0.2).^2);  % Peak mid
C_true(:,3) = exp(-((sampleIndex - 0.8)/0.2).^2);  % Peak late

% --- Create Gaussian spectral profiles (S_true) ---
S_true = zeros(nComponents, nVariables);
S_true(1,:) = exp(-((varIndex - 450)/60).^2);   % Component 1: 450 nm
S_true(2,:) = exp(-((varIndex - 550)/40).^2);   % Component 2: 550 nm
S_true(3,:) = exp(-((varIndex - 650)/60).^2);   % Component 3: 650 nm

% --- Generate mixture data: D = C * S + noise ---
D_true = C_true * S_true;

% Add Gaussian noise (SNR ≈ 100)
noiseLevel = 0.025;             % fraction of std noise added to D
noise = noiseLevel * randn(nSamples, nVariables);
D = D_true + noise;

% Optional SNR estimate (Frobenius norm based)
if noiseLevel > 0
    snr_db = 20*log10( norm(D_true,'fro') / max(norm(noise,'fro'), eps) );
    fprintf('Estimated SNR: %.2f dB\n', snr_db);
end

fprintf('Data generated successfully.\n');
fprintf('  Samples     : %d\n', nSamples);
fprintf('  Variables   : %d\n', nVariables);
fprintf('  Components  : %d\n', nComponents);
fprintf('  Noise level : %.2f%%\n\n', noiseLevel * 100);

%% ========================================================================
%% 2. Initialize MCR-ALS Lite
%% ========================================================================
fprintf('Initializing MCR-ALS Lite...\n\n');

% Random nonnegative initialization
C_init = abs(rand(nSamples, nComponents));

%% ========================================================================
%% 3. Run MCR-ALS Lite
%% ========================================================================
fprintf('Running MCR-ALS Lite...\n');
fprintf('========================================\n\n');

maxIter = 200;
tol     = 1e-10;

% Ensure function is available
assert(exist('MCR_ALS_Lite','file')==2, 'MCR_ALS_Lite.m not found on path.');

[C, S, lof] = MCR_ALS_Lite(D, C_init, [], maxIter, tol);

%% ========================================================================
%% 3b. Normalize BOTH true and recovered by their OWN L2 norms
%%     (component-wise; recovered versions kept separate for reconstruction)
%% ========================================================================
fprintf('Normalizing BOTH true and recovered by component-wise L2 norms...\n\n');

% Keep originals for reconstruction (C_rec * S_rec ≈ D)
C_rec = C;
S_rec = S;

% ---- True solution L2-normalization ----
% C_true: normalize columns (components)
normC_true = sqrt(sum(C_true.^2, 1));            % 1 x k
normC_true = max(normC_true, eps);
C_true_n = C_true ./ normC_true;                 % R2016b+ implicit expansion
% S_true: normalize rows (components)
normS_true = sqrt(sum(S_true.^2, 2));            % k x 1
normS_true = max(normS_true, eps);
S_true_n = S_true ./ normS_true;

% ---- Recovered solution L2-normalization ----
% C: normalize columns (components)
normC = sqrt(sum(C.^2, 1));                      % 1 x k
normC = max(normC, eps);
Cn = C ./ normC;
% S: normalize rows (components)
normS = sqrt(sum(S.^2, 2));                      % k x 1
normS = max(normS, eps);
Sn = S ./ normS;

% NOTE: After L2-normalization, Cn*Sn ~= D. Use C_rec, S_rec for reconstruction.

%% ========================================================================
%% 3c. Determine proper component order using corrcoef (every pair)
%%     We compute correlations on BOTH spectra and concentrations, then
%%     average them to form a matching score. No toolboxes; manual greedy.
%% ========================================================================
fprintf('Computing pairwise correlations and matching components...\n\n');

corrS = zeros(nComponents, nComponents); % rows: true, cols: recovered
corrC = zeros(nComponents, nComponents);
for i = 1:nComponents
    for j = 1:nComponents
        % spectra correlation (rows)
        corrS(i,j) = safeCorr(S_true_n(i,:), Sn(j,:));
        % concentration correlation (columns)
        corrC(i,j) = safeCorr(C_true_n(:,i), Cn(:,j));
    end
end

% Combine (average absolute correlations) to be robust
score = (abs(corrS) + abs(corrC)) / 2;

% Manual greedy assignment maximizing score with unique columns
reorderIdx = zeros(1, nComponents);
available = true(1, nComponents);
for i = 1:nComponents
    [~, jbest] = max(score(i,:));
    while ~available(jbest)
        score(i, jbest) = -Inf;           % mark as unusable and retry
        [~, jbest] = max(score(i,:));
    end
    reorderIdx(i) = jbest;
    available(jbest) = false;
end
assert(numel(unique(reorderIdx))==nComponents, 'Matching did not produce unique assignment.');

% Reorder recovered (normalized AND original) according to the match
Cn_match   = Cn(:, reorderIdx);
Sn_match   = Sn(reorderIdx, :);
C_rec_match = C_rec(:, reorderIdx);
S_rec_match = S_rec(reorderIdx, :);

% Report correlations for the chosen assignment
fprintf('Matching (true i  ->  recovered j):\n');
for i = 1:nComponents
    j = reorderIdx(i);
    fprintf('  i=%d -> j=%d   |  corrS=%.4f, corrC=%.4f\n', ...
        i, j, corrS(i,j), corrC(i,j));
end
fprintf('\n');

% Optional: visualize the score matrix (comment out if not needed)
figure('Name', 'Matching score | Average abs corr (C & S)', 'Position', [120, 120, 520, 420]);
imagesc(score, [0 1]); colorbar; axis image;
xlabel('Recovered component j'); ylabel('True component i');
title('Matching score = mean(|corr_S|, |corr_C|)');
set(gca,'XTick',1:nComponents,'YTick',1:nComponents);

%% ========================================================================
%% 4. Compare Results with True Profiles (L2-normalized & matched)
%% ========================================================================
fprintf('========================================\n');
fprintf('Comparing recovered (L2-normalized & matched) vs. true profiles...\n\n');

figure('Name', 'MCR-ALS Lite Results vs True Profiles (L2-normalized & matched)', ...
       'Position', [100, 100, 1400, 800]);

% --- Plot Concentration Profiles ---
for i = 1:nComponents
    subplot(2, nComponents, i);
    hold on;
    plot(C_true_n(:,i), 'k-', 'LineWidth', 2, 'DisplayName', 'True (L2)');
    plot(Cn_match(:,i), 'r--',  'LineWidth', 2, 'DisplayName', 'Recovered (L2)');
    hold off;
    xlabel('Sample Index', 'FontSize', 10);
    ylabel('Concentration (L2 norm = 1)', 'FontSize', 10);
    title(sprintf('C - Component %d', i), 'FontSize', 11, 'FontWeight', 'bold');
    legend('Location', 'best');
    grid on;
end

% --- Plot Spectral Profiles ---
for i = 1:nComponents
    subplot(2, nComponents, nComponents + i);
    hold on;
    plot(varIndex, S_true_n(i,:), 'k-', 'LineWidth', 2, 'DisplayName', 'True (L2)');
    plot(varIndex, Sn_match(i,:), 'b--',  'LineWidth', 2, 'DisplayName', 'Recovered (L2)');
    hold off;
    xlabel('Wavelength (nm)', 'FontSize', 10);
    ylabel('Intensity (L2 norm = 1)', 'FontSize', 10);
    title(sprintf('S - Component %d', i), 'FontSize', 11, 'FontWeight', 'bold');
    legend('Location', 'best');
    grid on;
end

sgtitle('MCR-ALS Lite: L2-normalized & Matched Recovered vs True Gaussian Profiles', ...
    'FontSize', 14, 'FontWeight', 'bold');

%% ========================================================================
%% 5. Reconstruction Quality (use ORIGINAL C_rec and S_rec, optionally matched)
%% ========================================================================
D_reconstructed = C_rec * S_rec;
residual = D - D_reconstructed;

figure('Name', 'MCR-ALS Lite Reconstruction Quality', 'Position', [150, 150, 1200, 420]);

subplot(1,3,1);
imagesc(D); colorbar;
title('Original Data (D)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Variable', 'FontSize', 10);
ylabel('Sample', 'FontSize', 10);

subplot(1,3,2);
imagesc(D_reconstructed); colorbar;
title('Reconstructed Data (C_{orig} * S_{orig})', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Variable', 'FontSize', 10);
ylabel('Sample', 'FontSize', 10);

subplot(1,3,3);
imagesc(residual); colorbar;
title('Residual (D - C_{orig} * S_{orig})', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Variable', 'FontSize', 10);
ylabel('Sample', 'FontSize', 10);

sgtitle(sprintf('MCR-ALS Lite Reconstruction (Final LOF = %.4f%%)', lof(end)), ...
    'FontSize', 14, 'FontWeight', 'bold');

fprintf('========================================\n');
fprintf('Test completed successfully!\n');
fprintf('Final LOF: %.4f%%\n', lof(end));
fprintf('========================================\n');

%% ========================================================================
%% Local helpers (robust correlation)
%% ========================================================================
function r = safeCorr(a, b)
    % Returns scalar correlation between vectors a and b.
    % If either vector has near-zero variance, returns 0 to avoid NaN.
    a = a(:); b = b(:);
    sa = std(a); sb = std(b);
    if sa < eps || sb < eps
        r = 0;
        return
    end
    C = corrcoef(a, b);
    r = C(1,2);
end
