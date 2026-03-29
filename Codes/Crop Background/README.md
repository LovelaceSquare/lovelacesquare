# cropBackground: Background Pixel Removal via Intensity Thresholding (MATLAB)

Crops background pixels from a 3D image cube by intensity thresholding and returns retained and discarded pixel indices.

**Reference**:  
Gómez-Sánchez, Adrián. (2024). *cropBackground function*. Lovelace’s Square.  
https://lovelacesquare.org/

---

## Overview

The `cropBackground` function processes hyperspectral or multi-channel image data stored as a 3D array (`rows × cols × channels`) by summing intensities across all channels to form a global intensity map.

Processing steps:
- Computes a per-pixel total intensity by summing along the third dimension.
- Thresholds the resulting intensity map using user-defined `minThresh` and `maxThresh`.
- Discards pixels outside this range (typically background or saturated regions).
- Converts spatial indices to linear indices for tracking.
- Outputs a cropped 2D matrix of size `(N_retained × channels)` containing only retained pixels.

This facilitates downstream analysis by focusing on meaningful spectral data, removing artifacts from low- or high-intensity background regions.

---

## Inputs

- `imageCube` (3D array): Hyperspectral or multi-channel image (`rows × cols × channels`)  
- `minThresh` (scalar): Minimum total intensity to retain a pixel  
- `maxThresh` (scalar): Maximum total intensity to retain a pixel  

## Outputs

- `croppedMatrix`: 2D matrix of retained pixels (`N_retained × channels`)  
- `keptIdx`: Linear indices of retained pixels  
- `removedIdx`: Linear indices of discarded pixels  

---

## Usage Example

(Insert into MATLAB script or command window)

    % Define thresholds to cut out dark background
    minT = 0.2 * max(grayImg(:));   % discard pixels darker than 20% of max
    maxT = max(grayImg(:));         % keep up to the brightest

    % Apply the cropBackground function
    [croppedMatrix, keptIdx, removedIdx] = cropBackground(imageCube, minT, maxT);

---

## Installation

### Prerequisites
- MATLAB R2016b or newer  
- No additional toolboxes required beyond core MATLAB

### Setup

1. Save `cropBackground.m` in a directory on your MATLAB path.
2. Add the folder if needed:

        addpath('path/to/your/functions');

3. Verify installation:

        which cropBackground

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

- **v1.0 (2024-12-14)**: Initial release; implements global intensity threshold cropping and returns pixel indices

---

## Keywords

- Image processing  
- Hyperspectral  
- Background removal  
- Thresholding  
- MATLAB  
- 3D image cube  
- Pixel cropping  
- Data preprocessing

