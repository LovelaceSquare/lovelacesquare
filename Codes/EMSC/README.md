# EMSC: Extended Multiplicative Scatter Correction for Spectral Data (MATLAB)

Performs Extended Multiplicative Scatter Correction on 2D spectral data, removing scatter and modeling baseline drift with polynomial terms.

**Reference**:  
Martens, Harald, and E. Stark. “Extended multiplicative signal correction and spectral interference subtraction: new preprocessing methods for near infrared spectroscopy.”  
*Journal of Pharmaceutical and Biomedical Analysis*, 9(8), 625–635 (1991).

---

## Overview

The `EMSC` function addresses the challenge of multiplicative and additive scatter effects in spectral datasets—common artifacts in near-infrared (NIR) and other spectroscopic measurements—while optionally correcting baseline drift.

It extends classic MSC by including polynomial terms in the regression model, fitting each spectrum `x_i` as:

    x_i(λ) = a_i + b_i·r(λ) + ∑_{p=1}^P c_{i,p}·λ^p + ε_i(λ)

Where:
- `r(λ)`: Reference spectrum (mean, median, or user-supplied)  
- `a_i`, `b_i`: Additive and multiplicative scatter coefficients  
- `c_{i,p}`: Polynomial baseline coefficients of order `P`  
- `ε_i(λ)`: Residual term  

The corrected spectrum is computed as:

    x̂_i(λ) = (x_i(λ) – a_i – ∑ c_{i,p}·λ^p) / b_i

This process improves the quality of spectra for downstream tasks such as calibration, classification, or multivariate decomposition.

---

## Inputs

- `X` (matrix): Spectral data matrix `[samples × variables]`  
- `refMode` (string): `'Mean'`, `'Median'`, or `'External'`  
- `polyOrder` (integer): Polynomial order for baseline modeling  
- `externalRef` (optional vector): Reference spectrum (required if `refMode = 'External'`)  

## Outputs

- `Xcorr`: Corrected spectra matrix `[samples × variables]`  
- `ref`: Reference spectrum used in the correction `[1 × variables]`  

---

## Usage Example

(Insert into MATLAB script or command window)

    % X is a [40×500] matrix of spectral measurements
    X = rand(40,500) + linspace(0,2,500);  % simulate baseline drift

    % 1) Mean reference, linear baseline (polyOrder = 1)
    [XcorrMean, refMean] = EMSC(X, 'Mean', 1);

    % 2) Median reference, quadratic baseline (polyOrder = 2)
    [XcorrMed, refMed] = EMSC(X, 'Median', 2);

    % 3) External reference with cubic baseline
    externalRef = load('cert_ref.mat','spectrum').spectrum;  % 1×500 vector
    [XcorrExt]()

