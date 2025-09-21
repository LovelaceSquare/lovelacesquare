# CosmicPeakCorrection: Spike Removal from Spectral Data via Derivative Thresholding (MATLAB)

Removes cosmic ray–induced spike artifacts from spectral data via derivative thresholding and interpolation in MATLAB.

**Reference**:  
Gómez-Sánchez, Adrián & Rocha de Oliveira, Rodrigo. (2024). *CosmicPeakCorrection function*. Lovelace’s Square.  
https://lovelacesquare.org/

---

## Overview

The `CosmicPeakCorrection` function addresses the common problem of narrow, high-amplitude spikes in spectral measurements, typically caused by cosmic ray impacts on detectors.

Key steps:
- Computes an Nth-order finite difference (derivative) along each spectrum (row-wise).
- Flags data points where the absolute derivative exceeds a user-defined threshold.
- Removes a specified number of neighboring channels around each detected spike.
- Reconstructs the signal using linear interpolation (with nearest-value boundary handling).

This method preserves true spectral features while effectively suppressing sharp spike artifacts that can interfere with downstream analyses like peak fitting or chemometric modeling.

---

## Inputs

- `X` (matrix): Spectral data matrix (rows = samples, columns = wavelengths)  
- `derivativeOrder` (integer): Order of the finite difference used for detection  
- `channelsToRemove` (integer): Number of neighboring channels (on each side) to remove around each detected spike  
- `threshold` (scalar): Absolute derivative magnitude threshold for spike detection  

## Outputs

- `cleanedData` (matrix): Corrected spectral data after spike removal  
- `peakMask` (logical matrix): Binary mask indicating locations of corrected (interpolated) points  

---

## Usage Example

(Insert into MATLAB script or command window)

    % Load or generate your spectral data: rows = samples, columns = wavelengths
    X = randn(50, 300) + sin(linspace(0,10,300));  

    % Set parameters
    derivativeOrder   = 1;    % 1st-order derivative sharpens single-channel spikes
    channelsToRemove  = 2;    % remove two channels on each side of detected peaks
    threshold         = 5.0;  % derivative magnitude threshold

    % Run cosmic spike correction
    [cleanedData, peakMask] = CosmicPeakCorrection(X, derivativeOrder, channelsToRemove, threshold);

    % Visualize results for the first sample
    figure;
    subplot(2,1,1);
    plot(X(1,:), 'k--', 'DisplayName','Original');
    hold on;
    plot(cleanedData(1,:), 'b-', 'DisplayName','Corrected');
    legend;
    xlabel('Wavelength Channel');
    ylabel('Intensity');
    title('Spectrum Before and After Cosmic Peak Correction');

    subplot(2,1,2);
    imagesc(peakMask);
    colorbar;
    xlabel('Channel');
    ylabel('Sample Index');
    title('Mask of Corrected Points');

---

## Installation

### Prerequisites
- MATLAB R2016a or later (uses built-in `diff`, `find`, `fillmissing`)

### Setup

1. Save `CosmicPeakCorrection.m` into a folder on your MATLAB path.
2. Add the folder:

        addpath('path/to/your/functions');

3. Verify availability:

        which CosmicPeakCorrection

### Dependencies
- No additional toolboxes required.

---

## License

Released under the **MIT License**.

---

## Authors

- **Adrián Gómez-Sánchez**  
- **Rodrigo Rocha de Oliveira**  
- **Created**: December 14, 2024  
- **Reviewed by**: Lovelace’s Square  

---

## Changelog

- **v1.0 (2024-12-14)**: Initial release, reviewed by Lovelace’s Square

---

## Keywords

- Cosmic ray correction  
- Spectral data cleaning  
- Spike removal  
- MATLAB  
- Derivative thresholding  
- Interpolation  
- Data preprocessing  
- Spectroscopy  
- Signal processing  
