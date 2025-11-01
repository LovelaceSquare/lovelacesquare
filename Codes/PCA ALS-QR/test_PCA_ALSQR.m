% test_PCA_ALSQR.m
%
% Demonstration script for PCA_ALSQR on synthetic spectral data with
% missing values. Compares reconstructed results to SVD on complete data.
%
% Author:         Adrián Gómez-Sánchez
% Date Created:   2025-03-01
% License:        MIT
% Reviewed by:    Lovelace's Square
% Version:        v 1.0
%

clear; clc; close all;

%% --- 1) Generate synthetic bilinear spectral data --------------------
nSamples   = 100;
nVariables = 1000;
nComp      = 3;

x = linspace(500, 3500, nVariables);
loadingsTrue = [
    exp(-0.5*((x-1000)/800).^2);
    exp(-0.5*((x-2000)/400).^2);
    exp(-0.5*((x-3000)/600).^2)
];
rng(123);
scoresTrue = rand(nSamples, nComp);
D_noiseless = scoresTrue * loadingsTrue;
D_complete  = D_noiseless + 0.01 * randn(nSamples, nVariables);

%% --- 2) Introduce missing data --------------------------------------
missingPercent = 50;
mask = rand(nSamples, nVariables) < (missingPercent/100);
D_incomplete = D_complete;
D_incomplete(mask) = NaN;

fprintf('Running PCA_ALSQR...\n');
[Dr, T, P, r2, lofc] = PCA_ALSQR(D_incomplete, 200, nComp);

fprintf('Final Lack of Fit: %.6f %%\n', lofc(end));

%% --- 3) Compare with SVD on complete data ----------------------------
[U, S, V] = svd(D_complete, 'econ');
Scores_svd   = U(:,1:nComp)*S(1:nComp,1:nComp);
Loadings_svd = V(:,1:nComp)';

%% --- 3b) Align signs of PCA_ALSQR with SVD ---------------------------
for k = 1:nComp
    corr_sign = sign( sum( Loadings_svd(k,:) .* P(k,:) ) );
    if corr_sign == 0, corr_sign = 1; end
    P(k,:) = P(k,:) * corr_sign;
    T(:,k) = T(:,k) * corr_sign;
end

fprintf('Final Lack of Fit: %.6f %%\n', lofc(end));

%% --- 4) Visualization -------------------------------------------------
figure('Name','PCA_ALSQR vs SVD','Color','w');
...


%% --- 4) Visualization -------------------------------------------------
figure('Name','PCA_ALSQR vs SVD','Color','w');

subplot(1,2,1);
plot(Scores_svd,'LineWidth',1.5); hold on;
plot(T,'--','LineWidth',1.2);
xlabel('Sample'); ylabel('Score Value');
title('Scores Comparison'); grid on;

subplot(1,2,2);
plot(x,Loadings_svd','LineWidth',1.5); hold on;
plot(x,P','--','LineWidth',1.2);
xlabel('Wavenumber (cm^{-1})'); ylabel('Loading Intensity');
title('Loadings Comparison'); grid on;

% --- Dynamic main title showing missing percentage --------------------
sgtitle(sprintf('SVD (Complete) vs PCA-ALS-QR (%.0f%% Missing)', missingPercent));

