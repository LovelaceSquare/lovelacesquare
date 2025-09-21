% test_fourierFilter.m
%
% This script generates synthetic spectral data that simulates a smooth 
% signal with superimposed high-frequency noise. The goal is to filter out 
% the noise by retaining only the low-frequency components (i.e., the smooth 
% spectral features) using a frequency-domain filtering approach.
%
% Author:  Adrián Gómez-Sánchez
% Date:    2024-12-14
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0

%% Setup Parameters
nSamples = 5;           % Number of spectral signals (rows)
nTimePoints = 1000;     % Number of data points per signal
dt = 0.001;             % Sampling interval in seconds (fs = 1/dt = 1000 Hz)
fs = 1/dt;              % Sampling frequency

% Define the time vector (e.g., wavelength or time)
t = (0:nTimePoints-1) * dt;

% Define the frequency interval to keep (low-pass filter to remove noise)
% Here, we choose to retain frequencies between 0 and 20 Hz.
freqIntervals = [0 20];

%% Generate Synthetic Spectral Data
% Create synthetic spectra that are smooth (e.g., a combination of sinusoids 
% representing broad spectral features) and add high-frequency noise.
data = zeros(nSamples, nTimePoints);
for i = 1:nSamples
    % Create a smooth baseline spectrum using a sum of low-frequency sinusoids.
    smoothSpectrum = sin(2*pi*5*t + 0.5*i) + 0.8*sin(2*pi*10*t + 0.3*i);
    
    % Add high-frequency noise (e.g., random fluctuations above 50 Hz)
    highFreqNoise = 0.5 * randn(1, nTimePoints);
    
    % Combine the smooth spectrum with the high-frequency noise
    data(i, :) = smoothSpectrum + highFreqNoise;
end

%% Apply Fourier Filter
% The fourierFilter function retains only the frequency components in the 
% specified interval (0–20 Hz) which should largely remove the high-frequency noise.
[filteredData, freqVector] = fourierFilter(data, freqIntervals, dt);

%% Visualization
% Select one spectrum for detailed visualization (change sampleIdx as needed)
sampleIdx = 1;

% Plot the original and filtered spectra in the time domain
figure;

subplot(3,1,1);
plot(t, data(sampleIdx, :), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Intensity');
title('Original Spectrum with High-Frequency Noise');
grid on;

subplot(3,1,2);
plot(t, filteredData(sampleIdx, :), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Intensity');
title('Filtered Spectrum (Low-Pass: 0-20 Hz)');
grid on;

% Plot the frequency spectrum (magnitude) of the original and filtered signals
nFFT = 2^nextpow2(nTimePoints);
origFFT = fft(data(sampleIdx, :), nFFT) / nTimePoints;
filteredFFT = fft(filteredData(sampleIdx, :), nFFT) / nTimePoints;
fAxis = fs/2 * linspace(0,1, nFFT/2+1);

subplot(3,1,3);
plot(fAxis, 2*abs(origFFT(1:nFFT/2+1)), 'b', 'LineWidth', 1.5); hold on;
plot(fAxis, 2*abs(filteredFFT(1:nFFT/2+1)), 'r', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('Original', 'Filtered');
title('Frequency Spectrum Comparison');
grid on;

%% Statistical Evaluation (Optional)
% Compute the energy (sum of squared magnitudes) in the pass-band for both
% the original and filtered spectra.
energyOrig = sum(abs(origFFT(1:nFFT/2+1)).^2);
energyFiltered = sum(abs(filteredFFT(1:nFFT/2+1)).^2);
fprintf('Total Energy (Original Spectrum in Passband): %.6f\n', energyOrig);
fprintf('Total Energy (Filtered Spectrum in Passband): %.6f\n', energyFiltered);

disp('Fourier filtering test completed. Review the plots and energy metrics.');
