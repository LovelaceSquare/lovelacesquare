% test.m
%
% This script generates several synthetic Raman spectra with varied features.
% For each spectrum, it introduces Gaussian noise and simulates cosmic peaks 
% by adding high-amplitude spikes. The CosmicPeakCorrection function is
% then applied to remove the cosmic peaks, and the results are visualized 
% using line plots for each individual spectrum.
%
% Authors: Adrián Gómez-Sánchez and Rodrigo Rocha de Oliveira
% Date Created: 2024-12-14
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0

%% Define Wavenumber Axis
% Typical Raman shift range (in cm^-1)
nWavenumbers = 1500;
wavenumber = linspace(500, 2000, nWavenumbers);

%% Generate Multiple Synthetic Raman Spectra with Distinct Features
nSamples = 3;  % Number of different spectra to generate
data = zeros(nSamples, nWavenumbers);

for sample = 1:nSamples
    % Randomly vary peak parameters for each spectrum:
    % Narrow peak parameters
    amp_narrow   = 90 + 20 * rand();          % Amplitude between 90 and 110
    center_narrow = 700 + 200 * rand();         % Center between 700 and 900 cm^-1
    sigma_narrow  = 8 + 4 * rand();             % Sigma between 8 and 12

    % Intermediate peak parameters
    amp_mid   = 40 + 20 * rand();               % Amplitude between 40 and 60
    center_mid = 1000 + 150 * rand();           % Center between 1000 and 1150 cm^-1
    sigma_mid  = 10 + 6 * rand();               % Sigma between 10 and 16

    % Broad peak parameters
    amp_wide   = 70 + 20 * rand();              % Amplitude between 70 and 90
    center_wide = 1300 + 400 * rand();           % Center between 1300 and 1700 cm^-1
    sigma_wide  = 20 + 15 * rand();             % Sigma between 20 and 35

    % Baseline offset (randomized)
    baseline  = 5 + 10 * rand();
    
    % Generate individual Gaussian peaks
    narrow_peak = amp_narrow * exp(-((wavenumber - center_narrow).^2) / (2*sigma_narrow^2));
    mid_peak    = amp_mid    * exp(-((wavenumber - center_mid).^2)    / (2*sigma_mid^2));
    wide_peak   = amp_wide   * exp(-((wavenumber - center_wide).^2)   / (2*sigma_wide^2));
    
    % Create the Raman spectrum by summing the peaks and adding the baseline.
    spectrum = narrow_peak + mid_peak + wide_peak + baseline;
    
    % Add Gaussian noise to simulate realistic measurement noise.
    noiseLevel = rand()*5;  % Noise level varying slightly between spectra.
    spectrum = spectrum + noiseLevel * randn(size(spectrum));
    
    % Save the generated spectrum.
    data(sample,:) = spectrum;
end

%% Simulate Cosmic Peaks
% For each sample, inject several random cosmic spikes.
numCosmicPeaks = 10;  % Number of cosmic peaks per spectrum.
rng(0);  % For reproducibility

for sample = 1:nSamples
    indices = randperm(nWavenumbers, numCosmicPeaks);
    for idx = indices
        % Add a large amplitude spike at the selected index.
        data(sample, idx) = data(sample, idx) + 200 * rand();
    end
end

%% Define Parameters for CosmicPeakCorrection
% Tune these parameters based on the synthetic data characteristics.
derivativeOrder  = 7;    % Compute first-order derivative.
channelsToRemove = 2;    % Remove 2 neighboring channels on each side of each detected spike.
threshold        = 500;   % Threshold for spike detection.

%% Apply Cosmic Peak Correction and Visualize Results
% For each sample, apply the correction and produce a figure with three subplots:
% 1. Original spectrum (with noise and cosmic peaks)
% 2. Detected cosmic peak mask (using a stem plot)
% 3. Corrected spectrum after cosmic peak removal
for sample = 1:nSamples
    [correctedData, peakMask] = CosmicPeakCorrection(data(sample,:), derivativeOrder, channelsToRemove, threshold);
    
    figure;
    
    % Plot Original Spectrum
    subplot(3,1,1);
    plot(wavenumber, data(sample,:), 'b-', 'LineWidth', 1.5);
    xlabel('Raman Shift (cm^{-1})');
    ylabel('Intensity');
    title(sprintf('Sample %d: Original Raman Spectrum with Noise & Cosmic Peaks', sample));
    grid on;
    
    % Plot Absolute Value of the Derivative with Threshold
    subplot(3,1,2);
    deriv = abs(diff(data(sample,:), derivativeOrder)); 
    x_deriv = wavenumber(1:length(deriv));  % Adjust x-axis for the derivative
    plot(x_deriv, deriv, 'LineWidth', 1.2);
    hold on;
    yline(threshold, 'r--', 'LineWidth', 1.5);  % Overlay the threshold line
    hold off;
    xlabel('Raman Shift (cm^{-1})');
    ylabel('Absolute Derivative');
    title(sprintf('Sample %d: Absolute Derivative with Threshold', sample));
    grid on;

    % Plot Corrected Spectrum
    subplot(3,1,3);
    plot(wavenumber, correctedData, 'g-', 'LineWidth', 1.5);
    xlabel('Raman Shift (cm^{-1})');
    ylabel('Intensity');
    title(sprintf('Sample %d: Corrected Raman Spectrum', sample));
    grid on;
end

%% Completion Message
disp('Raman spectrum cosmic peak correction test completed. Please review the generated plots for each sample.');