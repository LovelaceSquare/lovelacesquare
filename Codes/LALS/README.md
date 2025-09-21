# LALS: Local Asymmetric Least Squares Baseline Correction (MATLAB)

The Local Asymmetric Least Squares (LALS) algorithm extends asymmetric least squares (AsLS) smoothing to one-dimensional data vectors by allowing different asymmetry parameters (`pVals`) and second-derivative smoothing penalties (`lambdasAsym`) within specified intervals. Outside these intervals, a uniform smoothing penalty (`lambdaWhit`) and a global first-derivative penalty (`mu`) are enforced.

Through an Iteratively Reweighted Least Squares (IRLS) procedure, the method updates a weight matrix based on residuals to adaptively estimate a smooth, locally tailored baseline.

**Reference of the original ALS**:  
Eilers, P. H. C., & Boelens, H. F. M. (2005). *Baseline correction with asymmetric least squares smoothing*. Leiden University Medical Centre Report, 1(1), 5.

---

## Overview

The `LALS` function implements a Local Asymmetric Least Squares (LALS) baseline correction algorithm tailored for one-dimensional data vectors, such as spectra or chromatograms.

Key features:
- **Interval-specific asymmetry (`pVals`)**: Different sensitivity to positive/negative residuals across intervals.  
- **Interval-specific smoothing (`lambdasAsym`)**: Finer curvature control in defined segments.  
- **Global penalties**:
  - `lambdaWhit`: uniform second-derivative penalty outside intervals  
  - `mu`: global first-derivative penalty discouraging abrupt slope changes  

The method solves the following weighted system in each IRLS iteration:

    (W + Dᵀ * diag(lambdaVec) * D + μ * Lᵀ * L) * b = W * y

Where:
- `y`: observed data vector  
- `W`: diagonal weight matrix updated per iteration  
- `D`, `L`: second- and first-order difference operators  
- `lambdaVec`: vector combining local and global smoothing terms  

Each iteration includes:
1. Weighted system assembly  
2. Solution of normal equations for new baseline estimate  
3. Convergence check using `tol` or `maxIter`  
4. Weight update based on sign of residuals and `pVals`  

This method is ideal for data with region-dependent baseline drift (e.g., multi-modal spectra).

---

## Inputs

- `y` (column vector): 1D signal to be baseline-corrected  
- `intervals` (n×2 matrix): Start–end indices for each local region  
- `pVals` (n×1 vector): Asymmetry parameters (0 < p < 1) per interval  
- `lambdasAsym` (n×1 vector): Second-derivative penalties per interval  
- `lambdaWhit` (scalar): Global second-derivative penalty outside intervals  
- `mu` (scalar): First-derivative global smoothness penalty  
- `maxIter` (integer): Maximum number of IRLS iterations  
- `tol` (scalar): Convergence tolerance  

## Outputs

- `estimatedBaseline`: Baseline vector of same length as `y`  
- `weights`: Final diagonal of the weight matrix used in IRLS  

---

## Usage Example

(Insert into MATLAB script or command window)

    %% Setup Parameters
    nRows = 5; nCols = 500;
    intervals = [140,160;280,320;375,425];
    pVals = [1e-4;1e-4;5e-4]; lambdasAsym = [1e5;1e5;1e5];
    lambdaWhit = 10; mu = 10; maxIter = 50; tol = 1e-6;
    baselineAmplitude = 0.4; noiseLevel = 0.05;
    rng(42);

    %% Generate Synthetic Data
    x = linspace(0,10,nCols);
    trueBaseline = zeros(nRows,nCols); data = zeros(nRows,nCols);
    for i=1:nRows
        spectrum = exp(-((x-3).^2)/0.1^2) + 0.8*exp(-((x-6).^2)/0.2^2) + 0.5*exp(-((x-8).^2)/0.2^2);
        baseline = baselineAmplitude*(0.5*sin(2*pi*x/20) + 0.1*x);
        data(i,:) = spectrum + baseline + noiseLevel*randn(1,nCols);
        trueBaseline(i,:) = baseline;
    end

    %% Apply LALS on One Spectrum
    spectrumIndex = 1; y = data(spectrumIndex,:)';
    [estimatedBaseline, weights] = LALS(y, intervals, pVals, lambdasAsym, lambdaWhit, mu, maxIter, tol);
    correctedSpectrum = y - estimatedBaseline;

    %% Visualization and Evaluation
    figure;
    subplot(3,1,1);
    plot(x, data(spectrumIndex,:), 'LineWidth',1.5); hold on;
    plot(x, trueBaseline(spectrumIndex,:), '--', 'LineWidth',1.5);
    title('Original Spectrum with True Baseline'); legend('Noisy','True Baseline');

    subplot(3,1,2);
    plot(x, data(spectrumIndex,:), 'LineWidth',1.5); hold on;
    plot(x, estimatedBaseline, '--', 'LineWidth',1.5);
    title('Estimated Baseline (LALS)'); legend('Noisy','Estimated Baseline');

    subplot(3,1,3);
    plot(x, correctedSpectrum, 'LineWidth',1.5);
    title('Baseline-Corrected Spectrum'); legend('Corrected');

    mseBaseline = mean((trueBaseline(spectrumIndex,:)' - estimatedBaseline).^2);
    fprintf('Mean Squared Error: %.6e', mseBaseline);

---

## Installation

### Prerequisites
- MATLAB R2016a or later (uses sparse matrices and solvers)

### Setup

1. Save `LALS.m` into a folder on your MATLAB path.
2. Add the folder:

        addpath('path/to/LALS');

3. Verify installation:

        which LALS

### Dependencies
- Built-in MATLAB functions: `spdiags`, matrix division (`\`), `warnings`, `norm`, `ones`, `length`

---

## License

Released under the **MIT License**.

---

## Authors

- **Adrián Gómez-Sánchez**  
- **Created**: December 16, 2024  

---

## Changelog

- **v1.0 (2024-12-16)**: Initial implementation of per-interval LALS baseline correction

---

## Keywords

- Baseline correction  
- Asymmetric least squares  
- AsLS  
- Local asymmetric least squares  
- LALS  
- Spectral preprocessing  
- MATLAB  
- Detrending


