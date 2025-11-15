% =========================================================================
% test_PARAFAC_ALS_Lite.m
% =========================================================================
%
% Test script for PARAFAC_ALS_Lite (Parallel Factor Analysis - Alternating Least Squares)
% Demonstrates the algorithm on synthetic 3-way data generated from
% Gaussian-shaped factor profiles.
%
% Author:   Adrián Gómez-Sánchez
% Date:     2025-11-01
% License:  MIT
% Reviewed: Lovelace's Square
% Version:  v 1.0
%
% -------------------------------------------------------------------------
% REFERENCES:
%   • Harshman, R. A. (1970).
%     "Foundations of the PARAFAC procedure: Models and conditions for an
%      'explanatory' multimodal factor analysis." UCLA Working Papers in Phonetics.
%
%   • Bro, R. (1997).
%     "PARAFAC. Tutorial and applications." Chemometrics and Intelligent
%     Laboratory Systems, 38(2), 149–171.
%
% =========================================================================
clear; clc; close all;

%% ========================================================================
%% 1. Generate Synthetic Gaussian 3-Way Data
%% ========================================================================
fprintf('Generating synthetic Gaussian 3-way tensor data...\n\n');

% Reproducibility
rng(42, 'twister');

% Tensor dimensions
I = 50;           % Mode 1 dimension (e.g., samples/time points)
J = 40;           % Mode 2 dimension (e.g., wavelengths)
K = 30;           % Mode 3 dimension (e.g., elution time)
R = 3;            % Number of components

mode1Index = linspace(0, 1, I);
mode2Index = linspace(400, 700, J);   % e.g., wavelength range (nm)
mode3Index = linspace(0, 1, K);

% --- Create Gaussian factor profiles (A_true, B_true, C_true) ---
A_true = zeros(I, R);
A_true(:,1) = exp(-((mode1Index - 0.2)/0.6).^2);  % Peak early
A_true(:,2) = exp(-((mode1Index - 0.5)/0.4).^2);  % Peak mid
A_true(:,3) = exp(-((mode1Index - 0.8)/0.1).^2);  % Peak late

B_true = zeros(J, R);
B_true(:,1) = exp(-((mode2Index - 450)/60).^2);   % Component 1: 450 nm
B_true(:,2) = exp(-((mode2Index - 550)/200).^2);   % Component 2: 550 nm
B_true(:,3) = exp(-((mode2Index - 650)/10).^2);   % Component 3: 650 nm

C_true = zeros(K, R);
C_true(:,1) = exp(-((mode3Index - 0.3)/0.5).^2);  % Component 1
C_true(:,2) = exp(-((mode3Index - 0.5)/0.15).^2); % Component 2
C_true(:,3) = exp(-((mode3Index - 0.7)/1).^2);  % Component 3

% --- Generate tensor: X = sum_r a_r ⊗ b_r ⊗ c_r + noise ---
X_true = zeros(I, J, K);
for r = 1:R
    for i = 1:I
        for j = 1:J
            for k = 1:K
                X_true(i,j,k) = X_true(i,j,k) + A_true(i,r) * B_true(j,r) * C_true(k,r);
            end
        end
    end
end

% Add Gaussian noise
noiseLevel = 0.025;             % fraction of std noise added to X
noise = noiseLevel * randn(I, J, K);
X = X_true + noise;

% Optional SNR estimate (Frobenius norm based)
if noiseLevel > 0
    snr_db = 20*log10( norm(X_true(:)) / max(norm(noise(:)), eps) );
    fprintf('Estimated SNR: %.2f dB\n', snr_db);
end

fprintf('Tensor generated successfully.\n');
fprintf('  Mode 1 (I) : %d\n', I);
fprintf('  Mode 2 (J) : %d\n', J);
fprintf('  Mode 3 (K) : %d\n', K);
fprintf('  Components : %d\n', R);
fprintf('  Noise level: %.2f%%\n\n', noiseLevel * 100);

%% ========================================================================
%% 2. Initialize PARAFAC-ALS Lite
%% ========================================================================
fprintf('Initializing PARAFAC-ALS Lite...\n\n');

% Random nonnegative initialization for mode-1
B_init = abs(rand(J, R));
C_init = abs(rand(K, R));

%% ========================================================================
%% 3. Run PARAFAC-ALS Lite
%% ========================================================================
fprintf('Running PARAFAC-ALS Lite...\n');
fprintf('========================================\n\n');

maxIter = 200;
tol     = 1e-10;

% Ensure function is available
assert(exist('PARAFAC_ALS_Lite','file')==2, 'PARAFAC_ALS_Lite.m not found on path.');

[A, B, C, lof] = PARAFAC_ALS_Lite(X, [], B_init, C_init, maxIter, tol);

%% ========================================================================
%% 3b. Normalize BOTH true and recovered by their OWN L2 norms
%%     (component-wise; recovered versions kept separate for reconstruction)
%% ========================================================================
fprintf('Normalizing BOTH true and recovered by component-wise L2 norms...\n\n');

% Keep originals for reconstruction
A_rec = A;
B_rec = B;
C_rec = C;

% ---- True solution L2-normalization ----
% A_true: normalize columns (components)
normA_true = sqrt(sum(A_true.^2, 1));            % 1 x R
normA_true = max(normA_true, eps);
A_true_n = A_true ./ normA_true;

% B_true: normalize columns (components)
normB_true = sqrt(sum(B_true.^2, 1));            % 1 x R
normB_true = max(normB_true, eps);
B_true_n = B_true ./ normB_true;

% C_true: normalize columns (components)
normC_true = sqrt(sum(C_true.^2, 1));            % 1 x R
normC_true = max(normC_true, eps);
C_true_n = C_true ./ normC_true;

% ---- Recovered solution L2-normalization ----
% A: normalize columns (components)
normA = sqrt(sum(A.^2, 1));                      % 1 x R
normA = max(normA, eps);
An = A ./ normA;

% B: normalize columns (components)
normB = sqrt(sum(B.^2, 1));                      % 1 x R
normB = max(normB, eps);
Bn = B ./ normB;

% C: normalize columns (components)
normC = sqrt(sum(C.^2, 1));                      % 1 x R
normC = max(normC, eps);
Cn = C ./ normC;

%% ========================================================================
%% 3c. Determine proper component order using corrcoef (all three modes)
%%     We compute correlations on A, B, and C, then average them to form
%%     a matching score. No toolboxes; manual greedy.
%% ========================================================================
fprintf('Computing pairwise correlations and matching components...\n\n');

corrA = zeros(R, R); % rows: true, cols: recovered
corrB = zeros(R, R);
corrC = zeros(R, R);
for i = 1:R
    for j = 1:R
        % mode-1 correlation (columns of A)
        corrA(i,j) = safeCorr(A_true_n(:,i), An(:,j));
        % mode-2 correlation (columns of B)
        corrB(i,j) = safeCorr(B_true_n(:,i), Bn(:,j));
        % mode-3 correlation (columns of C)
        corrC(i,j) = safeCorr(C_true_n(:,i), Cn(:,j));
    end
end

% Combine (average absolute correlations) to be robust
score = (abs(corrA) + abs(corrB) + abs(corrC)) / 3;

% Manual greedy assignment maximizing score with unique columns
reorderIdx = zeros(1, R);
available = true(1, R);
for i = 1:R
    [~, jbest] = max(score(i,:));
    while ~available(jbest)
        score(i, jbest) = -Inf;           % mark as unusable and retry
        [~, jbest] = max(score(i,:));
    end
    reorderIdx(i) = jbest;
    available(jbest) = false;
end
assert(numel(unique(reorderIdx))==R, 'Matching did not produce unique assignment.');

% Reorder recovered (normalized AND original) according to the match
An_match = An(:, reorderIdx);
Bn_match = Bn(:, reorderIdx);
Cn_match = Cn(:, reorderIdx);
A_rec_match = A_rec(:, reorderIdx);
B_rec_match = B_rec(:, reorderIdx);
C_rec_match = C_rec(:, reorderIdx);

% Report correlations for the chosen assignment
fprintf('Matching (true i  ->  recovered j):\n');
for i = 1:R
    j = reorderIdx(i);
    fprintf('  i=%d -> j=%d   |  corrA=%.4f, corrB=%.4f, corrC=%.4f\n', ...
        i, j, corrA(i,j), corrB(i,j), corrC(i,j));
end
fprintf('\n');

% Optional: visualize the score matrix (comment out if not needed)
figure('Name', 'Matching score | Average abs corr (A, B, C)', 'Position', [120, 120, 520, 420]);
imagesc(score, [0 1]); colorbar; axis image;
xlabel('Recovered component j'); ylabel('True component i');
title('Matching score = mean(|corr_A|, |corr_B|, |corr_C|)');
set(gca,'XTick',1:R,'YTick',1:R);

%% ========================================================================
%% 4. Compare Results with True Profiles (L2-normalized & matched)
%% ========================================================================
fprintf('========================================\n');
fprintf('Comparing recovered (L2-normalized & matched) vs. true profiles...\n\n');

figure('Name', 'PARAFAC-ALS Lite Results vs True Profiles (L2-normalized & matched)', ...
       'Position', [100, 100, 1400, 900]);

% --- Plot Factor A (Mode 1) ---
for i = 1:R
    subplot(3, R, i);
    hold on;
    plot(A_true_n(:,i), 'k-', 'LineWidth', 2, 'DisplayName', 'True (L2)');
    plot(An_match(:,i), 'r--',  'LineWidth', 2, 'DisplayName', 'Recovered (L2)');
    hold off;
    xlabel('Mode-1 Index', 'FontSize', 10);
    ylabel('Loading (L2 norm = 1)', 'FontSize', 10);
    title(sprintf('A - Component %d', i), 'FontSize', 11, 'FontWeight', 'bold');
    legend('Location', 'best');
    grid on;
end

% --- Plot Factor B (Mode 2) ---
for i = 1:R
    subplot(3, R, R + i);
    hold on;
    plot(mode2Index, B_true_n(:,i), 'k-', 'LineWidth', 2, 'DisplayName', 'True (L2)');
    plot(mode2Index, Bn_match(:,i), 'b--',  'LineWidth', 2, 'DisplayName', 'Recovered (L2)');
    hold off;
    xlabel('Wavelength (nm)', 'FontSize', 10);
    ylabel('Intensity (L2 norm = 1)', 'FontSize', 10);
    title(sprintf('B - Component %d', i), 'FontSize', 11, 'FontWeight', 'bold');
    legend('Location', 'best');
    grid on;
end

% --- Plot Factor C (Mode 3) ---
for i = 1:R
    subplot(3, R, 2*R + i);
    hold on;
    plot(C_true_n(:,i), 'k-', 'LineWidth', 2, 'DisplayName', 'True (L2)');
    plot(Cn_match(:,i), 'g--',  'LineWidth', 2, 'DisplayName', 'Recovered (L2)');
    hold off;
    xlabel('Mode-3 Index', 'FontSize', 10);
    ylabel('Loading (L2 norm = 1)', 'FontSize', 10);
    title(sprintf('C - Component %d', i), 'FontSize', 11, 'FontWeight', 'bold');
    legend('Location', 'best');
    grid on;
end

sgtitle('PARAFAC-ALS Lite: L2-normalized & Matched Recovered vs True Gaussian Profiles', ...
    'FontSize', 14, 'FontWeight', 'bold');

%% ========================================================================
%% 5. Reconstruction Quality (use ORIGINAL A_rec, B_rec, C_rec, optionally matched)
%% ========================================================================
X_reconstructed = reconstruct_tensor(A_rec_match, B_rec_match, C_rec_match);
residual = X - X_reconstructed;

% Visualize slices
figure('Name', 'PARAFAC-ALS Lite Reconstruction Quality', 'Position', [150, 150, 1400, 420]);

% Show middle slice along mode-3
slice_k = round(K/2);

subplot(1,3,1);
imagesc(squeeze(X(:,:,slice_k))); colorbar;
title(sprintf('Original X (slice k=%d)', slice_k), 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Mode-2', 'FontSize', 10);
ylabel('Mode-1', 'FontSize', 10);

subplot(1,3,2);
imagesc(squeeze(X_reconstructed(:,:,slice_k))); colorbar;
title(sprintf('Reconstructed X (slice k=%d)', slice_k), 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Mode-2', 'FontSize', 10);
ylabel('Mode-1', 'FontSize', 10);

subplot(1,3,3);
imagesc(squeeze(residual(:,:,slice_k))); colorbar;
title(sprintf('Residual (slice k=%d)', slice_k), 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Mode-2', 'FontSize', 10);
ylabel('Mode-1', 'FontSize', 10);

sgtitle(sprintf('PARAFAC-ALS Lite Reconstruction (Final LOF = %.4f%%)', lof(end)), ...
    'FontSize', 14, 'FontWeight', 'bold');

fprintf('========================================\n');
fprintf('Test completed successfully!\n');
fprintf('Final LOF: %.4f%%\n', lof(end));
fprintf('========================================\n');

%% ========================================================================
%% Local helpers
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

function X_hat = reconstruct_tensor(A, B, C)
    % Reconstructs tensor from PARAFAC factors
    % X ≈ sum_{r=1}^R a_r ⊗ b_r ⊗ c_r
    [I, R] = size(A);
    [J, ~] = size(B);
    [K, ~] = size(C);

    X_hat = zeros(I, J, K);
    for r = 1:R
        % Outer product: a_r ⊗ b_r ⊗ c_r
        for i = 1:I
            for j = 1:J
                for k = 1:K
                    X_hat(i,j,k) = X_hat(i,j,k) + A(i,r) * B(j,r) * C(k,r);
                end
            end
        end
    end
end
