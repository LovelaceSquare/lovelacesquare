# MSC: Multiplicative Scatter Correction for Spectral Data (MATLAB)

Performs classic Multiplicative Scatter Correction (MSC) on row-wise spectral data in MATLAB.

**Reference**:  
Isaksson, T., & Næs, T. (1988). *The effect of multiplicative scatter correction (MSC) and linearity improvement in NIR spectroscopy*.  
Applied Spectroscopy, 42(7), 1273–1284.

---

## Overview

The `MSC` function removes both additive and multiplicative scatter effects from two-dimensional spectral data (`samples × wavelengths`) by aligning each spectrum to a reference.

Supported reference modes:
1. **Mean Spectrum** (default): uses the average across all samples  
2. **Reference Index**: uses a specified row from the input matrix  
3. **Reference Spectrum**: uses a user-supplied reference vector  

Each spectrum `x_i` is modeled as:

    x_i = a_i + b_i·r + ε_i

Where:
- `r`: reference spectrum  
- `a_i`: additive offset  
- `b_i`: multiplicative scale  
- `ε_i`: residual  

The corrected spectrum is computed as:

    (x_i – a_i) / b_i

This produces spectra corrected for scattering distortions.

---

## Inputs

- `X` (matrix): Spectral data `[samples × variables]`  
- `mode` (optional, string): `'Mean'` (default), `'Reference Index'`, or `'Reference Spectrum'`  
- `param` (optional): Sample index (integer) or reference vector `[1 × variables]`, depending on `mode`  

## Outputs

- `Xcorr`: MSC-corrected spectra `[samples × variables]`  
- `ref`: Reference spectrum used `[1 × variables]`  

---

## Usage Example

(Insert into MATLAB script or command window)

    % X is a [50×400] matrix of spectral measurements
    X = rand(50,400);

    % 1) Mean spectrum reference (default)
    [XcorrMean, refMean] = MSC(X);

    % 2) Specific sample as reference (sample #10)
    [Xcorr10, ref10] = MSC(X, 'Reference Index', 10);

    % 3) External spectrum reference
    externalSpec = rand(1,400);
    [XcorrExt, refExt] = MSC(X, 'Reference Spectrum', externalSpec);

    % Inspect outputs
    size(XcorrMean)   % → [50 400]
    size(refMean)     % → [1 400]

---

## Installation

### Prerequisites
- MATLAB R2016a or later

### Setup

1. Save `MSC.m` and its helper subfunction `msc_correction` into a folder on your MATLAB path.
2. Add the folder:

        addpath('path/to/functions');

3. Verify installation:

        which MSC

### Dependencies
- Built-in MATLAB functions only

---

## License

Released under the **MIT License**.

---

## Authors

- **Adrián Gómez-Sánchez**  
- **Created**: December 14, 2024  
- **Reviewed & Updated**: Lovelace’s Square – Version 1.3 (July 14, 2025)

---

## Changelog

- **v1.3 (2025-07-14)**: Added explicit-spectrum reference option; retained full documentation  
- **v1.0 (2024-12-14)**: Initial release and review by Lovelace’s Square

---

## Keywords

- Multiplicative scatter correction  
- MSC  
- NIR spectroscopy  
- Spectral preprocessing  
- MATLAB  
- Chemometrics  
- Baseline correction  
- Scatter removal  
- Data normalization  
- Signal correction
