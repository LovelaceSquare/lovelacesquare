% test_pure.m
%
% Deterministic showcase:
%   • Builds data matrix D with known pure *samples*
%   • Applies SIMPLISMA to D' (row direction) to find pure samples
%   • Plots D(samples,:)' (original left panel) and prints truth vs found


% Author: Adrián Gómez-Sánchez
% Date Created: 2025-08-05
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0
% -------------------------------------------------------------------------


clear, clc, close all

%% Parameters
nSamples   = 50;
nVars      = 400;
nComp      = 3;        % number of components = number of pure samples we expect
noisePct   = 1;        % SIMPLISMA noise threshold (% of max mean)
sigmaNoise = 0.05;
rng(0);

%% Synthetic components (three isolated peaks just to make mixtures realistic)
x = linspace(0,10,nVars);
pureSpectra = [
    exp(-((x-2.0).^2)/(2*1^2));
    exp(-((x-5.5).^2)/(2*3^2));
    exp(-((x-8.0).^2)/(2*0.7^2))
];
pureSpectra = pureSpectra ./ max(pureSpectra,[],2);

%% Concentrations: first nComp rows are pure samples (ground truth)
C = rand(nSamples, nComp);
C = C ./ sum(C,2);
C(1:nComp,:) = eye(nComp);     % rows 1..nComp are pure
truthSampleIdx = 1:nComp;

%% Data matrix
D = C * pureSpectra + sigmaNoise*randn(nSamples,nVars);

%% SIMPLISMA only in row direction (find pure samples)
if exist('Simplisma','file') == 2
    [~, sampIdx] = Simplisma(D.', nComp, noisePct);  % run on D'
else
    [~, sampIdx] = pure(D.', nComp, noisePct);       % fallback
end

%% Minimal display-only summary
disp('--- SIMPLISMA (row direction only) ---');
disp('Truth SAMPLE indices (rows):');  disp(truthSampleIdx);
disp('Found SAMPLE indices (rows):');  disp(sampIdx);
samps_ok = isequal(sort(sampIdx(:)), sort(truthSampleIdx(:)));
fprintf('Pure SAMPLES correctly identified: %s\n', string(samps_ok));

%% Plot (original left panel retained)
figure('Name','SIMPLISMA – pure samples (row direction only)', ...
       'Units','normalized','Position',[.18 .2 .64 .5]);
plot(1:nVars, D(sampIdx,:).', 'LineWidth',1.4);
title('Spectra of pure samples (found)');
xlabel('Variable (wavelength) index');
ylabel('Intensity');
legend("Sample "+sampIdx, 'Location','best');
grid on