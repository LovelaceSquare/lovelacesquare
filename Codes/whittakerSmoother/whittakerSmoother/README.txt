# WhittakerSmoother: Penalized Least Squares Smoothing (MATLAB)

Applies Whittaker smoothing to each row of a 2D data matrix using penalized least squares with customizable smoothing strength and difference order.

**Reference**:  
Whittaker, E. T. (1922). *On a new method of graduation*.  
*Proceedings of the Edinburgh Mathematical Society*, 41, 63–75.  

---

## Overview

The `WhittakerSmoother` function implements the Whittaker smoother—a penalized least-squares technique that balances fidelity to the data with smoothness of the output signal.

For each input row vector `y`, the smoothed signal `z` is found by solving:

    (I + lambda * D' * D) * z = y

Where:

- `I`: Identity matrix of size equal to number of columns  
- `D`: Finite-difference operator of order `d`, penalizing sharp fluctuations  
- `lambda`: Smoothing parameter — higher values enforce smoother results (penalize deviation from a polynomial of degree `d - 1`)  

This method:

- Is applied row-wise to an `[nRows × nCols]` matrix  
- Suppresses noise while preserving general shape and trends  
- Solves a sparse linear system for each row  
- Is ideal for moderate-sized spectra or time-series datasets

---

## Inputs

- `noisyMatrix` (matrix): Input matrix `[nRows × nCols]` of noisy signals  
- `lambda` (scalar): Smoothing parameter  
- `d` (integer): Order of the finite-difference operator

## Outputs

- `smoothedMatrix`: Output matrix of same size, with smoothed signals

---

## Usage Example

Paste into MATLAB:

```matlab
% Demonstrate WhittakerSmoother on noisy spectral data
nRows = 5;             % Number of spectra (rows)
nCols = 500;           % Number of data points per spectrum (columns)
lambda = 50;           % Smoothing parameter
d = 1;                 % Order of the finite difference operator
noiseLevel = 0.1;      % Standard deviation of added Gaussian noise
x = linspace(0, 10, nCols);
inputMatrix = zeros(nRows, nCols);

for i = 1:nRows
    % Generate a spectrum with multiple Gaussian peaks
    inputMatrix(i, :) = ...
        exp(-((x - 3).^2) / (2 * 0.2^2)) + ...
        0.5 * exp(-((x - 6).^2) / (2 * 0.5^2)) + ...
        0.3 * exp(-((x - 8).^2) / (2 * 0.3^2));
    inputMatrix(i, :) = inputMatrix(i, :) + 0.1 * sin(2 * pi * x / 10);
end

noisyMatrix = inputMatrix + noiseLevel * randn(size(inputMatrix));

%% Apply Whittaker Smoother
smoothedMatrix = whittakerSmoother(noisyMatrix, lambda, d);

%% Visualization
spectrumIndex = 1;  % Index of the spectrum to visualize
figure;
subplot(3, 1, 1);
plot(x, inputMatrix(spectrumIndex, :), 'LineWidth', 1.5);
title('Original Spectrum');
xlabel('Wavelength (arbitrary units)');
ylabel('Intensity');
subplot(3, 1, 2);
plot(x, noisyMatrix(spectrumIndex, :), 'LineWidth', 1.5);
title('Noisy Spectrum');
xlabel('Wavelength (arbitrary units)');
ylabel('Intensity');
subplot(3, 1, 3);
plot(x, smoothedMatrix(spectrumIndex, :), 'LineWidth', 1.5);
title('Smoothed Spectrum');
xlabel('Wavelength (arbitrary units)');
ylabel('Intensity');
```

---

## Installation

### Prerequisites
- MATLAB R2016a or later

### Setup

1. Save `WhittakerSmoother.m` into a directory on your MATLAB path.
2. Add the directory:

    ```matlab
    addpath('path/to/WhittakerSmoother');
    ```

3. Confirm installation:

    ```matlab
    which WhittakerSmoother
    ```

### Dependencies
- Built-in MATLAB functions: `speye`, `diff`, and standard matrix operations

---

## License

Released under the **MIT License**

---

## Authors

- **Adrián Gómez-Sánchez** (GitHub: [@adriangomez](https://github.com/adriangomez))  
- **Date of Creation**: December 14, 2024

---

## Changelog

- **v1.0 (2024-12-14)**:  
  Initial implementation of Whittaker smoothing with parameterized `lambda` and difference order.

---

## Keywords

- Whittaker smoother  
- penalized least squares  
- smoothing  
- MATLAB  
- spectral smoothing  
- signal processing  
- finite differences  
- noise reduction  
- trend estimation  

