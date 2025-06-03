%% Test Script for Autoscale Function
% This script tests the functionality of the autoscale function.
%
% Author:  Adrián Gómez-Sánchez
% Date Created: 2025-04-27
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v1.1

%% 1) Parameters for Simulated NIR Spectral Data
nSamples = 30;          % Number of spectra (rows)
nWavelengths = 1500;    % Number of wavelengths (columns)
wavelengthMin = 800;    % Minimum wavelength in nm
wavelengthMax = 2500;   % Maximum wavelength in nm
wavelengths = linspace(wavelengthMin, wavelengthMax, nWavelengths); % Wavelength axis

%% 2) Generate Realistic NIR Spectra with Gaussian Peaks
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
        spectrum = spectrum + amplitude * exp(-0.5 * ((wavelengths - peakPositions(p)) / peakWidths(p)).^2);
    end

    % Ensure non-negativity (already ensured by Gaussian peaks)
    trueSpectra(i, :) = spectrum;
end

%% 3) Simulate Multiplicative and Additive Effects
% Define known multiplicative factors and additive offsets

% Multiplicative factors: random values between 0.85 and 1.15
multiplicativeFactors = 0.85 + (rand(nSamples, 1) * 0.3); % 0.85 to 1.15

% Additive offsets: random values between -0.5 and 0.5
additiveOffsets = -0.5 + rand(nSamples, 1);

% Generate raw spectra by applying multiplicative and additive effects
% Broadcasting multiplicativeFactors and additiveOffsets across wavelengths
rawSpectra = (trueSpectra .* multiplicativeFactors) + additiveOffsets;

%% 4) Apply Autoscale
% Perform autoscaling in both 'column' and 'row' directions

% (a) Autoscale by columns (features)
[scaledSpectraColumn, scalingParamsColumn] = autoscale(rawSpectra, 'column');

% (b) Autoscale by rows (samples)
[scaledSpectraRow, scalingParamsRow] = autoscale(rawSpectra, 'row');

%% 5) Verification of Autoscaling
% Verify that the scaled data has zero mean and unit standard deviation

% Define a numerical tolerance
tolerance = 1e-10;

% (a) Verification for column-wise scaling
scaledMeansColumn = mean(scaledSpectraColumn, 1); % Mean of each column
scaledStdDevsColumn = std(scaledSpectraColumn, 0, 1);  % Std of each column

% Check that all column means are approximately zero
if all(abs(scaledMeansColumn) < tolerance)
    fprintf('Verification Passed: All column means are approximately zero.\n');
else
    numFailed = sum(abs(scaledMeansColumn) >= tolerance);
    warning('Verification Failed: %d/%d column means are not zero.', numFailed, nWavelengths);
end

% Check that all column standard deviations are approximately one
if all(abs(scaledStdDevsColumn - 1) < tolerance)
    fprintf('Verification Passed: All column standard deviations are approximately one.\n');
else
    numFailed = sum(abs(scaledStdDevsColumn - 1) >= tolerance);
    warning('Verification Failed: %d/%d column standard deviations are not one.', numFailed, nWavelengths);
end

% (b) Verification for row-wise scaling
scaledMeansRow = mean(scaledSpectraRow, 2); % Mean of each row
scaledStdDevsRow = std(scaledSpectraRow, 0, 2);  % Std of each row

% Check that all row means are approximately zero
if all(abs(scaledMeansRow) < tolerance)
    fprintf('Verification Passed: All row means are approximately zero.\n');
else
    numFailed = sum(abs(scaledMeansRow) >= tolerance);
    warning('Verification Failed: %d/%d row means are not zero.', numFailed, nSamples);
end

% Check that all row standard deviations are approximately one
if all(abs(scaledStdDevsRow - 1) < tolerance)
    fprintf('Verification Passed: All row standard deviations are approximately one.\n');
else
    numFailed = sum(abs(scaledStdDevsRow - 1) >= tolerance);
    warning('Verification Failed: %d/%d row standard deviations are not one.', numFailed, nSamples);
end

%% 6) Visualization of Results
% Plot original and scaled spectra for both scaling directions

% (a) Original Raw Spectra
figure('Position', [100, 100, 1600, 900]);

subplot(2,2,1);
plot(wavelengths, rawSpectra', 'Color', [0.7, 0.7, 0.7]);
title('Original Raw Spectra with Multiplicative and Additive Effects');
xlabel('Wavelength (nm)'); ylabel('Intensity (a.u.)');
xlim([wavelengthMin, wavelengthMax]);
grid on;

% (b) Scaled Spectra After Column-wise Autoscaling
subplot(2,2,2);
plot(wavelengths, scaledSpectraColumn', 'b-', 'LineWidth', 1.2);
title('Scaled Spectra After Column-wise Autoscaling');
xlabel('Wavelength (nm)'); ylabel('Scaled Intensity (a.u.)');
xlim([wavelengthMin, wavelengthMax]);
grid on;

% (c) Scaled Spectra After Row-wise Autoscaling
subplot(2,2,3);
plot(wavelengths, scaledSpectraRow', 'r-', 'LineWidth', 1.2);
title('Scaled Spectra After Row-wise Autoscaling');
xlabel('Wavelength (nm)'); ylabel('Scaled Intensity (a.u.)');
xlim([wavelengthMin, wavelengthMax]);
grid on;

% (d) Overlay of a Sample's Raw and Scaled Spectra (Both Directions)
subplot(2,2,4);
sampleIdx = randi(nSamples);  % Randomly select a sample to compare
plot(wavelengths, rawSpectra(sampleIdx, :), 'Color', [0.7, 0.7, 0.7], 'LineWidth', 1.2);
hold on;
plot(wavelengths, scaledSpectraColumn(sampleIdx, :), 'b-', 'LineWidth', 1.5);
plot(wavelengths, scaledSpectraRow(sampleIdx, :), 'r--', 'LineWidth', 1.5);
title(sprintf('Sample #%d: Raw vs. Column-wise vs. Row-wise Scaled', sampleIdx));
xlabel('Wavelength (nm)'); ylabel('Intensity (a.u.)');
legend('Raw Spectrum', 'Column-wise Scaled', 'Row-wise Scaled', 'Location', 'Best');
xlim([wavelengthMin, wavelengthMax]);
grid on;
hold off;

% (Optional) Additional Figure: Comparison of Means and Std Deviations
figure('Position', [1800, 100, 1200, 600]);

subplot(1,2,1);
plot(wavelengths, mean(rawSpectra, 1), 'Color', [0.7, 0.7, 0.7], 'LineWidth', 1.2);
hold on;
plot(wavelengths, scalingParamsColumn.mean, 'b-', 'LineWidth', 1.5);
title('Mean of Raw vs. Column-wise Scaled Data');
xlabel('Wavelength (nm)'); ylabel('Mean Intensity (a.u.)');
legend('Raw Mean', 'Scaled Mean', 'Location', 'Best');
xlim([wavelengthMin, wavelengthMax]);
grid on;
hold off;

subplot(1,2,2);
plot(wavelengths, std(rawSpectra, 0, 1), 'Color', [0.7, 0.7, 0.7], 'LineWidth', 1.2);
hold on;
plot(wavelengths, scalingParamsColumn.std, 'b-', 'LineWidth', 1.5);
title('Std Dev of Raw vs. Column-wise Scaled Data');
xlabel('Wavelength (nm)'); ylabel('Std Dev Intensity (a.u.)');
legend('Raw Std Dev', 'Scaled Std Dev', 'Location', 'Best');
xlim([wavelengthMin, wavelengthMax]);
grid on;
hold off;

%% 7) Display Summary of Scaling Parameters
fprintf('\n--- Autoscale Function Test Summary ---\n');
fprintf('Number of Samples: %d\n', nSamples);
fprintf('Number of Wavelengths: %d\n', nWavelengths);
fprintf('Column-wise Scaling:\n');
fprintf('  - Verification Passed: %d/%d column means are approximately zero.\n', sum(abs(scaledMeansColumn) < tolerance), nWavelengths);
fprintf('  - Verification Passed: %d/%d column standard deviations are approximately one.\n', sum(abs(scaledStdDevsColumn - 1) < tolerance), nWavelengths);
fprintf('Row-wise Scaling:\n');
fprintf('  - Verification Passed: %d/%d row means are approximately zero.\n', sum(abs(scaledMeansRow) < tolerance), nSamples);
fprintf('  - Verification Passed: %d/%d row standard deviations are approximately one.\n', sum(abs(scaledStdDevsRow - 1) < tolerance), nSamples);
fprintf('Autoscaling completed successfully.\n');
