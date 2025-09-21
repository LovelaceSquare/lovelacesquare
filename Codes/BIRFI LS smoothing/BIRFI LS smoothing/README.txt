# birfi_ls_smoothing: Smoothed IRF Estimation via Tikhonov-Regularized Hankel LS (MATLAB)

Estimates the instrument response function (IRF) from fluorescence decay data by solving a Tikhonov-regularized Hankel system with smoothing.

**Reference**:  
Gómez-Sánchez, Adrián, et al. (2024). *Blind instrument response function identification from fluorescence decays*. Biophysical Reports, 4(2).

---

## Overview

The `birfi_ls_smoothing` function implements a blind IRF estimation algorithm for time-resolved fluorescence decay measurements.

It constructs a Hankel matrix from the measured decay signal and solves the regularized normal equations:

    (Hᵀ H + λ Dᵀ D) · h = Hᵀ d

Where:
- `H`: Hankel matrix derived from the decay vector `d`  
- `D`: Second-order finite-difference operator enforcing smoothness  
- `λ`: Tikhonov regularization parameter (penalty weight)  
- `h`: Estimated IRF vector of length `irf_size`  

By balancing data fidelity (via `Hᵀ H`) and smoothness (via `Dᵀ D`), the method recovers a stable, smooth IRF estimate even in the presence of noise.  
This approach is particularly suited for cases where the instrument response is unknown or direct measurement is impractical.

---

## Inputs

- `decay_signal` (row vector): Measured fluorescence decay  
- `irf_length` (integer): Desired length of the IRF vector  
- `penalty` (optional, scalar): Smoothing regularization parameter; defaults to `1e8`  

## Output

- `t_est`: Estimated IRF (row vector of length `irf_length`)  

---

## Usage Example

(Insert into MATLAB script or command window)

    % Given a fluorescence decay signal 'decay_signal' (1×N):
    decay_signal = load('decay_example.mat');  % user-provided data
    irf_length   = 100;          % desired IRF vector length
    penaltyParam = 1e8;          % smoothing regularization weight

    % Estimate the IRF
    t_est = birfi_ls_smoothing(decay_signal, irf_length, penaltyParam);

    % Plot the estimated IRF
    t = (0:irf_length-1);
    figure;
    plot(t, t_est, 'LineWidth', 1.5);
    xlabel('Time Bin Index');
    ylabel('IRF Amplitude');
    title('Estimated Instrument Response Function');
    grid on;

    % Using default regularization (penalty = 1e8)
    irf_default = birfi_ls_smoothing(decay_signal, irf_length);

---

## Installation

### Prerequisites
- MATLAB R2016b or later  
- No additional toolboxes required beyond base MATLAB

### Setup

1. Save `birfi_ls_smoothing.m` to a directory on your MATLAB path.
2. Add the folder:

        addpath('path/to/functions');

3. Verify availability:

        which birfi_ls_smoothing

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

- **v1.0 (2024-12-16)**: Initial implementation of Tikhonov-regularized IRF estimation with smoothing

---

## Keywords

- IRF estimation  
- Instrument response function  
- Tikhonov regularization  
- Hankel matrix  
- Smoothing penalty  
- Second-order difference  
- Fluorescence decay  
- MATLAB  
- Inverse problem  
- Signal processing

````
