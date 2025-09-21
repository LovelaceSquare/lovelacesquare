## Title & Short Description & Reference

**imageRegistration & WarpImageSimilarity (with demo_ImageRegistration)** – A small set of MATLAB files that line up (“register”) a moving image with a reference image, even if the moving one has been shifted, rotated, or resized. Includes a ready-to-run demo that creates a fake 100-band hyperspectral scene, applies a known transformation, and then recovers it automatically.

**Reference**: Piqueras, Sara, et al. "Handling different spatial resolutions in image fusion by multivariate curve resolution-alternating least squares for incomplete image multisets." *Analytical Chemistry*, 90.11 (2018): 6757–6765. https://doi.org/10.1021/acs.analchem.8b00630

## Overview

### Why you might need this
- Images taken at different times, angles, or resolutions often do not align perfectly—they may be shifted (translation), rotated, or zoomed (scaled).
- To compare or fuse such images, they must be aligned pixel by pixel—a process known as **registration**.

### What’s inside
- **`imageRegistration.m`** – Searches for the best four transformation parameters: horizontal shift (`tx`), vertical shift (`ty`), rotation (`theta` in degrees), and zoom (`scale`).
  - Starts from either a default or user-defined initial guess, with optional interactive clicking for rough alignment.
  - Keeps parameter values within user-defined limits during the search.
  - Measures alignment quality using the squared difference only in areas where both images contain valid (non-NaN) values.

- **`WarpImageSimilarity.m`** – Applies the final transformation to move, rotate, and resize the image. Works on both grayscale and hyperspectral (multi-band) data.

- **`demo_ImageRegistration.m`** – Demonstrates the full workflow:
  1. Builds a simulated hyperspectral cube (256×256×100) with three geometric shapes and known spectral profiles.
  2. Applies a known transformation (shift, rotation, scaling).
  3. Registers summed grayscale images to estimate that transformation.
  4. Applies the recovered transformation to the entire hyperspectral cube.
  5. Displays numerical accuracy and visual overlay of the result.

### How it works (in plain language)
1. **Initialize the guess** – Either use a default `[0 0 0 1]` (no transformation) or manually select corresponding points for a better start.
2. **Optimize** – MATLAB’s built-in `fminsearch` adjusts the parameters to minimize the mismatch between the moving and fixed images.
3. **Warp the image** – Once the best parameters are found, the image is transformed. Each pixel in the output image is mapped back to its source position in the original image using inverse transformation. Bilinear interpolation estimates sub-pixel values, and pixels outside the image bounds are set to `NaN`.
4. **Repeat for every band** – For hyperspectral images, the same transformation is applied to all bands to maintain spatial consistency.

## Usage

```matlab
%% 1. Quick manual run with mouse picking
fixed  = imread('reference.png');    % your reference
moving = imread('to_align.png');     % image to align
[registered, pBest] = imageRegistration(moving, fixed);
imshowpair(fixed, registered);       % red/green overlay

%% 2. Fully automatic with your own limits (no GUI)
p0 = [0 0 0 1];
lb = [-150 -150 -30 0.8];    % lower bounds
ub = [ 150  150  30 1.2];    % upper bounds
[registered, pBest] = imageRegistration(moving, fixed, p0, lb, ub, false);

%% 3. Align a hyperspectral cube using the recovered parameters
load hyperCube.mat            % variable: hyperCube (HxWxN)
aligned = zeros(size(hyperCube));
for k = 1:size(hyperCube,3)
    aligned(:,:,k) = WarpImageSimilarity(hyperCube(:,:,k), pBest);
end

%% 4. Full demonstration
demo_ImageRegistration;
```

### What you’ll see in the demo
- Console output with the applied distortion, the recovered parameters, and the percent error.
- Three figures:
  1. Original summed image (fixed)
  2. Distorted summed image (moving)
  3. Overlay (red = fixed, green = aligned); yellow = good match

Typical output:
```
Ground truth:       [60.0  -40.0  20.0°  1.200]
Recovered transform [-57.4  19.1  -19.9° 0.8333]
Errors:             [0.18% 0.33% 0.48% 0.00%]
✓ Registration successful!
```


## License

MIT License – free to use, share, and modify with attribution.

## Authors & Contact

Adrián Gómez-Sánchez
Created: 2025-08-02  
Reviewed by: Lovelace’s Square

## Changelog

- **v1.0 (2025-08-02)** – Initial public release

## Keywords

image, registration, alignment, MATLAB, rotation, translation, scaling, hyperspectral imaging