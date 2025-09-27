# WienerFiltering: Frequency-Domain Spectral Denoising (MATLAB)

Applies Wiener filtering to 2D spectral datasets by estimating signal and noise power spectra and performing frequency-domain noise reduction using Welch’s method.

**Reference**:  
Wiener, N. (1949). *Extrapolation, Interpolation, and Smoothing of Stationary Time Series*. MIT Press.

---

## Overview

The `WienerFiltering` function implements a row-wise Wiener filter for spectral data matrices, leveraging power spectral density (PSD) estimation via Welch’s method (`pwelch`).

For each input spectrum:

1. **PSD Estimation**  
   Computes the two-sided PSD of the noisy signal using specified segment length, overlap, and FFT parameters.

2. **Noise Model**  
   Assumes flat noise PSD equal to the provided `noiseVar` estimate.

3. **Signal PSD Extraction**  
   Derives the signal PSD by subtracting the noise PSD from the observed PSD, clipping negative values to zero.

4. **Filter Design**  
   Constructs the Wiener filter transfer function:

       H(f) = Sxx(f) / (Sxx(f) + Snn(f))

   Where:
   - `Sxx(f)`: estimated signal power spectrum  
   - `Snn(f)`: noise power spectrum  

5. **Frequency-Domain Filtering**  
   Applies `H(f)` to the FFT of the original signal and performs inverse FFT to recover the filtered spectrum.

This approach adaptively attenuates noise-dominated frequencies while preserving signal features, improving signal-to-noise ratios in each spectrum.

---

## Inputs

- `data` (matrix): Input matrix `[nSamples × nPoints]` of noisy spectra  
- `noiseVar` (scalar): Estimated noise variance (assumed flat across frequency)  
- `segmentLength` (integer): Segment length for `pwelch`  
- `overlap` (integer): Overlap between segments  
- `nfft` (integer): Number of FFT points

## Outputs

- `filteredData`: Output matrix with filtered spectra (same size as input)

---

## Usage Example

Paste into MATLAB:

```matlab
% Demonstrate WienerFiltering on synthetic spectra

% Generate sample data: 4 noisy Gaussian spectra
nSamples = 4; nPoints = 512;
x = linspace(400, 700, nPoints);
data = zeros(nSamples, nPoints);
for i = 1:nSamples
    data(i,:) = exp(-((x-550).^2)/(20 + 5*i)^2) + 0.05*randn(1, nPoints);
end

% Wiener filter parameters
noiseVar = 0.05;         % estimated noise variance
segmentLength = 128;     % pwelch segment length
overlap = 64;            % overlap between segments
nfft = 512;              % FFT length

% Apply Wiener filtering
filteredData = WienerFiltering(data, noiseVar, segmentLength, overlap, nfft);

% Plot original vs. filtered for first spectrum
figure;
plot(x, data(1,:), 'k--', 'DisplayName', 'Original'); hold on;
plot(x, filteredData(1,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Filtered');
legend; xlabel('Wavelength (nm)'); ylabel('Intensity');
title('Wiener Filtering Example');
hold off;
```

---

## Installation

### Prerequisites
- MATLAB R2016a or later  
- Signal Processing Toolbox (required for `pwelch`)

### Setup

1. Save `WienerFiltering.m` into a directory on your MATLAB path.
2. Add the directory:

    ```matlab
    addpath('path/to/WienerFiltering');
    ```

3. Confirm availability:

    ```matlab
    which WienerFiltering
    ```

### Dependencies
- MATLAB built-in functions: `pwelch`, `fft`, `ifft`, standard array operations

---

## License

Released under the **MIT License**

---

## Authors

- **Adrián Gómez-Sánchez** (Email: gomez.sanchez.adr@gmail.com)  
- **Date of Creation**: December 14, 2024

---

## Changelog

- **v1.0 (2024-12-14)**:  
  Initial implementation of Wiener filtering with Welch PSD estimation and configurable parameters.

---

## Keywords

- Wiener filter  
- spectral denoising  
- power spectral density  
- Welch’s method  
- MATLAB  
- noise reduction  
- frequency-domain filtering  
- signal processing  
- spectral preprocessing  
- stationary time series

