# fourierFilter: Frequency-Domain Filtering for Time-Series Signals (MATLAB)

Applies frequency-domain filtering to time-series signals using FFT with interval masking.

**Reference**:  
Gómez-Sánchez, Adrián. (2024). *fourierFilter function*. Lovelace’s Square.  
https://lovelacesquare.org/

---

## Overview

The `fourierFilter` function implements a flexible and efficient frequency-domain filtering method for 2D time-series data, where each row represents an individual signal.

Processing steps:
1. **FFT computation**: Each signal is zero-padded to the next power of two and transformed to the frequency domain using FFT.  
2. **Frequency masking**: A mask is created over the FFT result, retaining only the frequencies specified in the `freqIntervals` parameter.  
3. **Symmetric filtering**: The mask is applied to both positive and corresponding negative frequencies to ensure real-valued outputs.  
4. **Signal reconstruction**: The filtered signal is converted back to the time domain using IFFT and truncated to the original signal length.  

This approach is particularly suitable for band-pass or band-stop filtering without relying on time-domain filter design. It includes automatic merging of overlapping intervals for robust frequency selection.

---

## Inputs

- `X` (matrix): Time-series data matrix `[signals × time points]`  
- `freqIntervals` (n×2 matrix): Frequency intervals to retain (in Hz), one row per range  
- `dt` (scalar): Sampling interval in seconds (e.g., 0.01 for fs = 100 Hz)  

## Outputs

- `Xfiltered`: Filtered signals (same size as `X`)  
- `fVec`: Frequency vector (in Hz) corresponding to FFT bins  

---

## Usage Example

(Insert into MATLAB script or command window)

    % Example data: 10 time-series signals with 1000 points each
    X = randn(10, 1000);

    % Define frequency ranges to retain: 0–15 Hz and 40–60 Hz
    intervals = [0 15; 40 60];

    % Sampling interval (seconds), e.g., 0.01 implies fs = 100 Hz
    dt = 0.01;

    % Apply frequency filtering
    [Xfiltered, fVec] = fourierFilter(X, intervals, dt);

    % Plot one example
    plot(fVec, abs(fft(X(1,:), 2^nextpow2(size(X,2))))(1:length(fVec))); hold on;
    plot(fVec, abs(fft(Xfiltered(1,:), 2^nextpow2(size(Xfiltered,2))))(1:length(fVec)));
    legend('Original','Filtered');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    title('Frequency Filtering Example');

---

## Installation

### Prerequisites
- MATLAB R2016b or newer  
- No external libraries required

### Setup

1. Save the `fourierFilter.m` file and the included helper `mergeOverlappingRanges` in your MATLAB path.
2. Add the folder if needed:

        addpath('path/to/your/functions');

3. Verify installation:

        which fourierFilter

---

## License

Released under the **MIT License**.

---

## Authors

- **Adrián Gómez-Sánchez**  
- **Created**: December 14, 2024  
- **Reviewed by**: Lovelace’s Square  

---

## Changelog

- **v1.0 (2024-12-14)**: Initial release with full FFT-based filtering and interval merging

---

## Keywords

- Frequency filtering  
- Fourier  
- Smoothing  
- FFT  
- MATLAB  
- Band-pass  
- Time-series  
- Signal processing  
- Inverse FFT

