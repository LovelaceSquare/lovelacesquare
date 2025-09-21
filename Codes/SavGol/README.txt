# SavGol: Savitzky–Golay Smoothing and Differentiation (MATLAB)

Applies Savitzky–Golay smoothing or differentiation to each row of a 2D data matrix with customizable window, polynomial, derivative order, and edge handling.

**Reference**:  
Savitzky, A., & Golay, M. J. E. (1964). *Smoothing and differentiation of data by simplified least squares procedures*.  
*Analytical Chemistry*, 36(8), 1627–1639.

---

## Overview

The `SavGol` function smooths or differentiates each row of a two-dimensional numeric matrix using the Savitzky–Golay filter—a moving-window least-squares polynomial fit method.

Key features include:

- **Window Size (`windowSize`)**: Odd-length span over which a polynomial is fit  
- **Polynomial Order (`polyOrder`)**: Degree of the fitted polynomial (must be less than `windowSize`)  
- **Derivative Order (`derivOrder`)**: Specifies zero for smoothing only or higher integer orders for numerical differentiation  
- **Edge Handling (`edgeMethod`)**: Strategies to manage boundary effects, including:
  - `'None'`: Zero padding  
  - `'Reflection'`: Reflect data at ends  
  - `'Replication'`: Replicate end values  
  - `'Extrapolation'`: Polynomial-based extrapolation  

`SavGol` constructs convolution coefficients via a pseudoinverse of a Vandermonde matrix for local polynomial fits, then applies these coefficients row-wise using convolution (`conv`). Edge regions are treated according to the selected `edgeMethod` to avoid artifacts.

---

## Inputs

- `data` (matrix): Input numeric matrix `[M × N]`  
- `windowSize` (odd integer): Length of the smoothing window  
- `polyOrder` (integer): Degree of polynomial (`polyOrder < windowSize`)  
- `derivOrder` (integer): Derivative order (`0` = smoothing)  
- `edgeMethod` (string): One of `'None'`, `'Reflection'`, `'Replication'`, `'Extrapolation'`

## Outputs

- `filteredData`: Matrix with smoothed or differentiated rows (same size as input)

---

## Usage Example

Paste into MATLAB:

```matlab
%% Setup Parameters
nSamples = 5;           % Number of spectra (rows)
nPoints = 500;          % Number of data points per spectrum (columns)
windowSize = 31;        % Must be odd and > polyOrder
polyOrder = 1;          % Polynomial order (e.g., linear)
derivOrder = 0;         % Derivative order (0 = smoothing)
edgeMethod = 'Extrapolation'; % Edge handling strategy: 'None', 'Reflection', etc.
noiseLevel = 0.1;       % Standard deviation of Gaussian noise

%% Generate Synthetic Spectral Data
x = linspace(0, 10, nPoints);
data = zeros(nSamples, nPoints);
for i = 1:nSamples
    % Generate a spectrum with multiple Gaussian peaks
    data(i, :) = ...
        exp(-((x - 3).^2) / (2 * 0.3^2)) + ... % Peak at x = 3
        0.8 * exp(-((x - 6).^2) / (2 * 0.4^2)) + ... % Broader peak at x = 6
        0.5 * exp(-((x - 8).^2) / (2 * 0.2^2)); % Narrow peak at x = 8
    data(i, :) = data(i, :) + 0.05 * sin(2 * pi * x / 2);
end

% Add Gaussian noise
noisyData = data + noiseLevel * randn(size(data));

%% Apply Savitzky-Golay Filter
filteredData = SavGol(noisyData, windowSize, polyOrder, derivOrder, edgeMethod);

%% Visualization
spectrumIndex = 1;  % Index of the spectrum to visualize
figure;
subplot(3, 1, 1);
plot(x, data(spectrumIndex, :), 'LineWidth', 1.5);
title('Original Spectrum');
xlabel('Wavelength (arbitrary units)');
ylabel('Intensity');
subplot(3, 1, 2);
plot(x, noisyData(spectrumIndex, :), 'LineWidth', 1.5);
title('Noisy Spectrum');
xlabel('Wavelength (arbitrary units)');
ylabel('Intensity');
subplot(3, 1, 3);
plot(x, filteredData(spectrumIndex, :), 'LineWidth', 1.5);
title('Filtered Spectrum');
xlabel('Wavelength (arbitrary units)');
ylabel('Intensity');
```

---

## Installation

### Prerequisites
- MATLAB R2016a or later

### Setup

1. Save `SavGol.m`, `computeSavGolCoeffs.m`, and `applySavGolRow.m` (subfunctions included in `SavGol.m`) into a folder on your MATLAB path.
2. Add the directory:

    ```matlab
    addpath('path/to/SavGol');
    ```

3. Confirm availability:

    ```matlab
    which SavGol
    ```

### Dependencies
- Built-in MATLAB functions: `pinv`, `conv`, `polyfit`, `polyval`, `factorial`, `spdiags`, and vector operations

---

## License

Released under the **MIT License**

---

## Authors

- **Adrián Gómez-Sánchez** 
- **Date of Creation**: December 14, 2024

---

## Changelog

- **v1.0 (2024-12-14)**:  
  Initial implementation of Savitzky–Golay filtering with edge handling options.

---

## Keywords

- Savitzky–Golay  
- smoothing  
- differentiation  
- SavGol  
- MATLAB  
- spectral preprocessing  
- moving window  
- least squares  
- edge handling  
- digital filter

