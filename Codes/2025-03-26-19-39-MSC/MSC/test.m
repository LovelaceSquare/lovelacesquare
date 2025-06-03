%% Test Script for MSC Function
% This script tests the functionality of the MSC function
%
% Author:  Adrián Gómez-Sánchez
% Date Created: 2024-12-14
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v1.0

%% 1) Parameters for Simulated NIR Spectral Data
nSamples = 20;          % Number of spectra (rows)
nWavelengths = 1000;    % Number of wavelengths (columns)
wavelengthMin = 800;    % Minimum wavelength in nm
wavelengthMax = 2500;   % Maximum wavelength in nm
wavelengths = linspace(wavelengthMin, wavelengthMax, nWavelengths); % Wavelength axis

%% 2) Generate Realistic "True" NIR Spectra
% Simulate NIR spectra using a combination of Gaussian peaks to mimic real absorption features.

% Define parameters for Gaussian peaks (positions in nm, amplitudes, and widths)
peakPositions = [950, 1200, 1400, 1900, 2100, 2300]; % Typical NIR absorption regions
peakAmplitudes = [1, 0.8, 1.2, 0.9, 1.1, 0.7];      % Amplitudes of the peaks
peakWidths = [50, 60, 55, 70, 65, 60];              % Standard deviations of the peaks

% Initialize the true spectra matrix
trueSpectra = zeros(nSamples, nWavelengths);

% Generate each sample's spectrum with slight variations in peak amplitudes
for i = 1:nSamples
    spectrum = zeros(1, nWavelengths);
    for p = 1:length(peakPositions)
        % Introduce slight random variations in amplitude for each sample
        ampVariation = 1 + 0.05 * randn(); % ±5% variation
        amplitude = peakAmplitudes(p) * ampVariation;

        % Generate Gaussian peak
        spectrum = spectrum + amplitude * exp(-0.5 * ((wavelengths - peakPositions(p))/peakWidths(p)).^2);
    end

    % Ensure non-negativity (already ensured by Gaussian peaks)
    trueSpectra(i, :) = spectrum;
end

%% 3) Simulate Multiplicative and Additive Effects
% Define known multiplicative factors and additive offsets

% Multiplicative factors: random values between 0.85 and 1.15
multiplicativeFactors = 0.85 + (rand(nSamples, 1) * 2);

% Additive offsets: non-negative small random values between 0 and 0.1
additiveOffsets = rand(nSamples, 1) * 3;
% Ensures that rawSpectra will not be negative since trueSpectra >= 0

% Generate raw spectra by applying multiplicative and additive effects
% Broadcasting multiplicativeFactors and additiveOffsets across wavelengths
rawSpectra = (trueSpectra .* multiplicativeFactors) + additiveOffsets;

%% 5) Apply MSC
% Use the mean of all spectra as the reference.
% (refIndex is irrelevant for 'Mean Spectrum', but must still be given)
[correctedSpectra, referenceSpec] = MSC(rawSpectra, 'Mean Spectrum', 1);

%% 6) Verification of MSC Correction
% Since MSC is intended to correct for multiplicative and additive scatter,
% the correctedSpectra should closely match the trueSpectra.

% Calculate the Mean Squared Error (MSE) between corrected and true spectra
mse = mean((correctedSpectra - trueSpectra).^2, 'all');



%% 7) Visualize Results
figure('Position', [100, 100, 1200, 800]);

% (a) True Spectra
subplot(3,1,1);
plot(wavelengths, trueSpectra', 'LineWidth', 1.2);
title('True Spectra (Underlying Signal)');
xlabel('Wavelength (nm)'); ylabel('Intensity (a.u.)');
xlim([wavelengthMin, wavelengthMax]);
grid on;

% (b) Raw Spectra with Additive and Multiplicative Effects
subplot(3,1,2);
plot(wavelengths, rawSpectra', 'LineWidth', 1.2);
title('Raw Spectra with Additive and Multiplicative Effects');
xlabel('Wavelength (nm)'); ylabel('Intensity (a.u.)');
xlim([wavelengthMin, wavelengthMax]);
grid on;

% (c) Corrected Spectra After MSC
subplot(3,1,3);
plot(wavelengths, correctedSpectra', 'LineWidth', 1.2);
title('Corrected Spectra After MSC');
xlabel('Wavelength (nm)'); ylabel('Intensity (a.u.)');
xlim([wavelengthMin, wavelengthMax]);
grid on;

% Optional: Overlay True and Corrected Spectra for a Specific Sample
figure('Position', [1500, 100, 800, 600]);
sampleIdx = randi(nSamples);  % Randomly select a sample to compare
plot(wavelengths, trueSpectra(sampleIdx, :), 'b-', 'LineWidth', 2);
hold on;
plot(wavelengths, correctedSpectra(sampleIdx, :), 'r--', 'LineWidth', 2);
title(sprintf('Comparison of True and Corrected Spectra for Sample #%d', sampleIdx));
xlabel('Wavelength (nm)'); ylabel('Intensity (a.u.)');
legend('True Spectrum', 'Corrected Spectrum', 'Location', 'Best');
xlim([wavelengthMin, wavelengthMax]);
grid on;
hold off;