% test_kernelize.m
%
% Generate a **rank-3 mono-exponential decay** data matrix, apply the
% `kernelize` function, and display basic sanity checks.
%
% Author:  Adrián Gómez-Sánchez
% Date Created: 2025-07-14
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v1.0

%% 1) Synthetic rank-3 mono-exponential data --------------------------------
nSamples  = 15;                 % rows
nTimePts  = 120;                % columns
t         = (0:nTimePts-1)';    % time vector

% Three decay constants (samples share these "basis" exponentials)
tau       = [2 10 35];

% Build 3×T matrix of exponentials
E = exp(-t./tau);               % each column scaled automatically (T × 3)
E = E.';                        % 3 × T

% Random positive mixing coefficients for each sample (A ≥ 0)
A = rand(nSamples, 3);    % amplitudes ~ U(0.8,1.8)

% Rank-3 data matrix + small Gaussian noise
D  = A * E + 0.01*randn(nSamples, nTimePts);

%% 2) Kernelise -------------------------------------------------------------
numK   = 20;        % number of kernels (≥2)
kWidth = 40;       % kernel length (≤ nTimePts)

Dk = kernelize(D, numK, kWidth);   % size → [15 × 20 × 96]

%% 3) Quick checks ----------------------------------------------------------
fprintf('Original  D size: %s\n', mat2str(size(D)));
fprintf('Kernelised Dk size: %s\n\n', mat2str(size(Dk)));

% Plot the first sample before and after kernelisation (one kernel slice)
figure;
subplot(1,2,1);
plot(t, D(1,:), 'LineWidth',1.2);
xlabel('Time index'); ylabel('Signal'); grid on; title('Sample #1 – raw');

subplot(1,2,2);
plot(squeeze(Dk(1,2,:)));                 % kernel #2 as example
xlabel('Convolved index'); ylabel('Response');
grid on; title('Sample #1 – after kernel #2');

%% Completion message
disp('Kernelisation test completed — inspect the figure and console output.');
