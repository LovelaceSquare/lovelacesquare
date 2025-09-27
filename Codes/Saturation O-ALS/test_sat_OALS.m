% test_sat_OALS.m.m
%
% PURPOSE
%   Demonstrate recovery of saturated spectral peaks by treating saturated
%   readings as missing (NaN) and reconstructing with OALS on a low-rank model.
%
% INTUITION 
%   Consider each spectrum as a point in a space where each variable
%   (wavelength/channel) is an axis. With 10,000 variables, spectra live in a
%   10,000-D space. In many chemometric problems the data are bilinear: a
%   measured spectrum can be represented as a combination of a few “pure”
%   spectra (components) with corresponding coefficients. Geometrically, the
%   data lie close to a low-dimensional (affine) subspace. PCA identifies and
%   rotates this subspace to describe the data with fewer dimensions while
%   retaining essential structure. Crucial for missing/saturated values:
%   even if a spectrum have missing variables, the point
%   still belongs to the same low-dimensional subspace. The missing variables 
%   are estimated based on this subspace.
%
% MODELING CHOICES
%   • This demonstration does NOT center the data (no mean-centering).
%   • Synthetic spectra are bilinear (rank = 3) with low Gaussian noise.
%   • Detector saturation is simulated via hard clipping at ADC_max.
%
% PRACTICAL CONSIDERATIONS
%   • Bilinear model and multiple samples.
%     – OALS assumes a bilinear (low-rank) data matrix with several spectra
%       (samples). Single-spectrum recovery is ill-posed.
%
%   • Coverage of saturated regions (informative missingness).
%     – A variable (wavelength/channel) that is saturated must be observed
%       unsaturated in at least one other sample to be recoverable. If a column
%       is saturated (or missing) in all samples, that region cannot be recovered.
%     – Long contiguous saturated segments reduce stability. Recovery improves
%       when saturated blocks are short relative to the local spectral structure
%       and when neighboring samples provide diverse, unsaturated information.
%
%   • Rank selection (number of components).
%     – The number of components (nPC) must be known or estimated beforehand.
%       Too small or too large → biased reconstructions. Use cross-validation, 
%       scree profiles, or external knowledge to set nPC.
%
%   • Selective peak strategy.
%     – If the target peak is selective (dominant and localized), restrict the
%       analysis to that spectral window and use a single principal component
%       (nPC = 1) to simplify and stabilize the reconstruction in that region.
%
%   • Centering and baselines.
%     – This demonstration does not center the data; the latent structure is
%       therefore affine. If substantial baselines/offsets or fluorescence are
%       present in real data, consider baseline correction.
%
%   • Saturation detection and instrument limits.
%     – Accurate identification of saturated points is required (knowledge of
%       ADC ceiling and any firmware clipping). Quantization and digitization
%       near the ceiling may introduce bias even in “unsaturated” readings.
%
%   • Identifiability vs. local minima.
%     – OALS solves a nonconvex problem; results can depend on initialization.
%       Use multiple restarts and compare final lack-of-fit. Ensure the pattern
%       of missingness provides enough constraints (rows and columns with
%       sufficient valid entries).
%
% OUTPUT
%   A two-panel figure:
%     LEFT  – Saturated (clipped) spectrum in black.
%     RIGHT – True pre-ceiling spectrum in blue, the same spectrum with NaNs
%             in black (gaps at saturated points), and the OALS reconstruction
%             in dashed red.
%   Console metrics: MSE vs the true pre-ceiling signal (global and on
%   saturated regions).
%
% REQUIREMENTS
%   OALS.m and its dependencies: ScoresLS.m, LoadingsLS.m, lofNaN.m
%
% Author:   Adrián Gómez-Sánchez
% Date:     2025-09-27
% License:  MIT
% Reviewed: Lovelace's Square
% Version:  v 1.2

clear; clc; close all;

%% 1) Synthetic bilinear spectra (no centering)
rng(20250927);
nSamples     = 60;
nVariables   = 1200;
nComponents  = 3;
wavenumber   = linspace(500, 3500, nVariables);

% Loadings: Gaussian peaks near 1000/2000/3000 cm^-1
g = @(x,mu,s) exp(-0.5*((x-mu)/s).^2);
Ptrue = [ g(wavenumber,1000,90);   % broad
          g(wavenumber,2000,45);   % narrow
          g(wavenumber,3000,65) ]; % medium

% Scores scaled to induce clipping after an overall scale-up
Ttrue = 0.2 + rand(nSamples, nComponents) .* [1.3 1.6 1.1];

% Bilinear signal with low Gaussian noise (1%); no mean-centering is applied
D_noiseless = Ttrue * Ptrue;
noiseLevel  = 0.01;
D_complete  = D_noiseless + noiseLevel * randn(nSamples, nVariables);

%% 2) Simulate detector saturation (hard clipping at ADC_max)
scaleFactor = 1.35;                   % pushes intensities toward ceiling
D_scaled    = scaleFactor * D_complete;    % "true" pre-ceiling signal

ADC_max   = 1.0;                      % detector ceiling
D_clipped = min(D_scaled, ADC_max);   % measured (flat tops at saturation)
satMask   = D_scaled >= ADC_max;      % saturated entries

% For OALS: treat saturated entries as missing values
D_for_OALS          = D_clipped;
D_for_OALS(satMask) = NaN;

sat_ratio = 100 * nnz(satMask) / numel(satMask);
fprintf('Saturated points: %.2f %% of all entries\n\n', sat_ratio);

%% 3) OALS on NaN-masked data
iter = 300;
[Dr_oals, ~, ~, ~, ~] = OALS(D_for_OALS, iter, nComponents);

%% 4) Select one representative spectrum (most saturation)
[~, worstIdx] = max(sum(satMask,2));
fprintf('Sample selected for plotting: %d\n\n', worstIdx);

%% 5) Two-panel figure (specified layering)
figure('Name','Saturated peak recovery (comparison)','Color','w');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% Common y-limits
ymin = min([D_clipped(worstIdx,:), D_scaled(worstIdx,:), Dr_oals(worstIdx,:)], [], 'omitnan');
ymax = max([D_clipped(worstIdx,:), D_scaled(worstIdx,:), Dr_oals(worstIdx,:)], [], 'omitnan');
pad  = 0.05*(ymax - ymin + eps);
yl   = [ymin - pad, ymax + pad];

% LEFT: saturated (clipped) spectrum in black
nexttile; hold on;
plot(wavenumber, D_clipped(worstIdx,:), 'k-', 'LineWidth', 1.6);
xlabel('Wavenumber (cm^{-1})'); ylabel('Intensity');
title('Saturated (clipped)'); box on;
xlim([wavenumber(1) wavenumber(end)]); ylim(yl);

% RIGHT: true (blue) FIRST, then NaNs (black), then OALS (red dashed)
nexttile; hold on;
plot(wavenumber, D_scaled(worstIdx,:),   'b-',  'LineWidth', 1.6, 'DisplayName','True (pre-ceiling)');
plot(wavenumber, D_for_OALS(worstIdx,:), 'k-',  'LineWidth', 1.6, 'DisplayName','Saturated with NaNs');
plot(wavenumber, Dr_oals(worstIdx,:),    'r--', 'LineWidth', 1.8, 'DisplayName','OALS reconstruction');
xlabel('Wavenumber (cm^{-1})'); ylabel('Intensity');
title('True (blue), NaNs (black), OALS (red dashed)');
legend('Location','best'); box on;
xlim([wavenumber(1) wavenumber(end)]); ylim(yl);

%% 6) Quantitative evaluation vs true pre-ceiling signal
mse = @(A,B,mask) mean( (A(mask) - B(mask)).^2 );
validAll        = true(size(D_scaled));
validSaturated  = satMask;

mse_clipped_all = mse(D_scaled, D_clipped, validAll);
mse_oals_all    = mse(D_scaled, Dr_oals,   validAll);
mse_clipped_sat = mse(D_scaled, D_clipped, validSaturated);
mse_oals_sat    = mse(D_scaled, Dr_oals,   validSaturated);

fprintf('MSE vs TRUE (pre-ceiling)\n');
fprintf('  Global:    Clipped = %.3e | OALS = %.3e\n', mse_clipped_all, mse_oals_all);
fprintf('  Saturated: Clipped = %.3e | OALS = %.3e\n', mse_clipped_sat,  mse_oals_sat);
