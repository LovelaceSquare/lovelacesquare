# birfi_ls: Tikhonov-Regularized IRF Estimation from Fluorescence Decays (MATLAB)

Estimates the instrument response function (IRF) from fluorescence decay data using Tikhonov-regularized Hankel least-squares.

**Reference**:  
Gómez-Sánchez, Adrián, et al. (2024). *Blind instrument response function identification from fluorescence decays*. Biophysical Reports, 4(2).

---

## Overview

The `birfi_ls` function performs blind IRF estimation by formulating the decay reconstruction as a regularized least-squares problem.

A Hankel matrix is constructed from the measured decay signal, and the IRF vector is solved from the system:

    (Hᵀ H + λ I) · h = Hᵀ d

Where:
- `H`: Hankel matrix built from the input decay vector `d`  
- `I`: Identity matrix of size `irf_size`  
- `λ`: Tikhonov regularization parameter (controls stability)  
- `h`: Estimated IRF of length `irf_size`  

By adding the identity term weighted by `λ`, the solution mitigates overfitting and noise amplification. The method is straightforward, computationally efficient, and suitable when the instrument response cannot be measured directly.

---

## Inputs

- `decay_signal` (row vector): Measured fluorescence decay  
- `irf_size` (integer): Desired number of time bins in the IRF  
- `lambda` (optional, scalar): Regularization parameter; defaults to `1e5`  

## Output

- `t_irf`: Estimated IRF as a row vector of length `irf_size`  

---

## Usage Example

(Insert into MATLAB script or command window)

    % Load or define a decay signal (row vector)
    decay_signal = load('decay_example.mat');  % user-provided fluorescence decay

    iSize = 50;          % desired IRF length (number of time bins)
    regParam = 1e5;      % Tikhonov regularization weight

    % Estimate the IRF
    t_irf = birfi_ls(decay_signal, iSize, regParam);

    % Plot the result
    timeBins = 0:(iSize-1);
    figure;
    plot(timeBins, t_irf, 'LineWidth', 1.5);
    xlabel('Time Bin Index');
    ylabel('IRF Amplitude');
    title('Estimated IRF via Tikhonov LS');
    grid on;

    % Using default regularization (lambda = 1e5)
    t_irf_default = birfi_ls(decay_signal, iSize);

---

## Installation

### Prerequisites
- MATLAB R2016b or newer  
- No additional toolboxes required

### Setup

1. Save `birfi_ls.m` into a directory on your MATLAB path.
2. Add the folder:

        addpath('path/to/functions');

3. Confirm availability:

        which birfi_ls

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

- **v1.0 (2024-12-16)**: Initial release with Tikhonov-regularized Hankel LS IRF estimation

---

## Keywords

- IRF estimation  
- Tikhonov regularization  
- Hankel matrix  
- Least squares  
- Fluorescence decay  
- Inverse problem  
- MATLAB  
- Regularized LS  
- Instrument response  
- Signal processing


