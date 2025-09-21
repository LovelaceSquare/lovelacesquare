# AsLS: Asymmetric Least Squares Baseline Correction (MATLAB)

Performs Asymmetric Least Squares (AsLS) baseline correction on spectral data using MATLAB.

**Reference**:  
Eilers, P. H. C., & Boelens, H. F. M. (2005). *Baseline correction with asymmetric least squares smoothing*. Leiden University Medical Centre Report, 1(1), 5.  

---

## Overview

The `AsLS` function implements the asymmetric least squares (AsLS) smoothing algorithm to estimate and subtract a smooth baseline from each row of a two-dimensional spectral dataset.

It minimizes the following objective function:

    ||w·(y – z)||² + λ·||D²z||²

Where:
- `y`: the original spectrum  
- `z`: the estimated baseline  
- `D²`: the second-order difference operator  
- `w`: a weight vector that penalizes positive residuals more than negative ones based on asymmetry parameter `p`  
- `λ`: smoothing parameter  

The algorithm iteratively updates the weights (`w`) over a fixed number of iterations (default: 10) until convergence.

---

## Inputs

- `X` (matrix): Spectral data (rows = spectra, columns = spectral channels)  
- `lambda` (scalar): Smoothing parameter (larger = smoother baseline)  
- `p` (scalar): Asymmetry parameter (`0 < p < 1`, smaller = more negative bias)  

## Outputs

- `correctedData`: The baseline-corrected spectra  
- `baseline`: The estimated baseline for each spectrum  

---

## Usage Example

(Insert into MATLAB script or command window)

    % Example: simulate noisy spectral data
    X = peaks(50, 200) + randn(50, 200) * 0.05;

    % Set AsLS parameters
    lambda = 1e6;    % smoothing strength
    p      = 0.001;  % asymmetry (smaller = more negative weight)

    % Apply baseline correction
    [correctedData, baseline] = AsLS(X, lambda, p);

    % Plot result for the first spectrum
    figure;
    plot(X(1,:), 'k--', 'DisplayName', 'Original');
    hold on;
    plot(baseline(1,:), 'r-', 'DisplayName', 'Estimated Baseline');
    plot(correctedData(1,:), 'b-', 'DisplayName', 'Corrected');
    legend;
    xlabel('Spectral Channel');
    ylabel('Intensity');
    title('AsLS Baseline Correction Example');
    hold off;

---

## Installation

### Prerequisites
- MATLAB R2016a or later

### Setup

1. Save `AsLS.m` (and subfunction `asls_baseline`) into a directory.
2. Add the folder to your MATLAB path:

        addpath('path/to/AsLS');

3. Verify:

        which AsLS

### Dependencies
- Only MATLAB built-ins: `diff`, `speye`, `spdiags`, standard matrix operations.

---

## License

Released under the **MIT License**.

---

## Authors

- **Adrián Gómez-Sánchez**  
- **Created**: December 16, 2024  
- **Reviewed by**: Lovelace’s Square  

---

## Changelog

- **v1.0 (2024-12-16)**: Initial release.

---

## Keywords

- Asymmetric Least Squares  
- Baseline correction  
- Spectral preprocessing  
- MATLAB  
- Smoothing  
- Weighted least squares  
- Signal processing  
- Iterative algorithm  
- Denoising

