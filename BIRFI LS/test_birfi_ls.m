% test_birfi_ls.m
%
% This script generates synthetic decay data with a known instrument 
% response function (IRF) with two shoulders using an exponential decay 
% convolved with the synthetic IRF, adds Poisson noise, applies birfi_ls 
% to estimate the IRF, and evaluates the results visually by comparing the 
% estimated IRF with the true IRF.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 2024-12-16
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: 1.0

%% Setup Parameters
N = 500;            % Length of the exponential decay signal
irf_size = 100;      % Fixed size for the IRF
tau = 50;            % Decay constant for the exponential decay
lambda = 3*10E6;    % Regularization parameter

% Generate the exponential decay signal
t = (1:N);
signal = 10000 * exp(-t/tau);  % Exponential decay with maximum ~10,000

%% Generate Synthetic True IRF
x_irf = 1:irf_size;
center = (irf_size + 1) / 2;
% Main central peak
main_peak = 0.6 * exp(-((x_irf - center).^2) / (5 * 2^2));
% Right shoulder (shifted right by 10 samples)
right_shoulder = 0.2 * exp(-((x_irf - (center + 10)).^2) / (5 * 2^2));
irf_true = main_peak + right_shoulder;
irf_true = irf_true / sum(irf_true);  % Normalize so that sum(irf_true) = 1

%% Generate Synthetic Decay Data
% Convolve the signal with the true IRF to simulate the decay
decay_clean = conv(signal, irf_true, 'full');

% Scale the decay so that its maximum equals 10,000 counts
scaling_factor = 10000 / max(decay_clean);
decay_clean = decay_clean * scaling_factor;

% Add Poisson noise to simulate counting statistics
decay_noisy = poissrnd(decay_clean);

%% Apply birfi_ls to Estimate IRF
% Note: The default lambda (1e-3) is used if not provided.
estimated_irf = birfi_ls(decay_noisy, irf_size, lambda);

%% Visualization
figure;

% Top subplot: Noisy decay data with true IRF overlaid
subplot(2,1,1);
x_decay = 1:length(decay_noisy);
plot(x_decay, decay_noisy, 'LineWidth', 1.5); hold on;
% Scale the true IRF to have a maximum of 10,000 counts for display
irf_display = irf_true / max(irf_true) * 10000;
plot(1:irf_size, irf_display, '--', 'LineWidth', 1.5);
title('Noisy Convolved Decay Data with True IRF Overlay');
xlabel('Sample Index');
ylabel('Counts');
legend('Noisy Decay', 'True IRF (Scaled)');
grid on;

% Bottom subplot: True IRF vs Estimated IRF (normalized)
subplot(2,1,2);
plot(1:irf_size, irf_true/norm(irf_true), 'LineWidth', 1.5); hold on;
plot(1:irf_size, estimated_irf/norm(estimated_irf), '--', 'LineWidth', 1.5);
title('True IRF vs Estimated IRF');
xlabel('Sample Index');
ylabel('Normalized Amplitude');
legend('True IRF', 'Estimated IRF');
grid on;

fprintf('Test completed. Review the plotted IRF comparison.\n');


