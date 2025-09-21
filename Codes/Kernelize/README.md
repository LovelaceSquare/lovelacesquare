# kernelize: Signal Kernelization via Exponential and Impulse Convolution (MATLAB)

Applies kernelization to multivariate signals by convolving each sample with a bank of normalized exponential and impulse kernels to enhance feature extraction in trilinear decomposition.

**Reference**:  
Gómez-Sánchez, Adrián, et al. (2023). *Kernelizing: A way to increase accuracy in trilinear decomposition analysis of multiexponential signals*.  
Analytica Chimica Acta, vol. 1273, article 341545.

---

## Overview

The `kernelize` function transforms a two-dimensional data matrix `D` (samples × timepoints) into a three-way array by applying a set of user-specified kernels via convolution.

It generates:
- **Impulse kernels** at the first and last positions to capture local spikes  
- **Paired decaying exponentials** mirrored about the center to detect transient behaviors  
- **A symmetric exponential** (if an odd number of kernels is requested)  

All kernels are normalized to unit maximum to ensure comparability.  
Convolution is performed using the `'valid'` option, producing an output length of:

    num_timepoints – kernel_width + 1

The result is a 3D array of size `[samples × kernels × timepoints′]`, ready for trilinear or multicomponent analysis.

---

## Inputs

- `D` (matrix): Input data matrix `[samples × timepoints]`  
- `num_kernels` (integer): Number of kernels to apply (≥2)  
- `kernel_width` (integer): Width of each kernel (≤ number of columns in `D`)  

## Output

- `D_kernelized`: 3D array `[samples × kernels × timepoints′]`  

---

## Usage Example

(Insert into MATLAB script or command window)

    % Example data: 10 samples, 100 timepoints
    D = randn(10,100);

    % Specify number of kernels and kernel width
    num_kernels  = 5;    % must be integer ≥2
    kernel_width = 20;   % ≤ number of columns in D

    % Apply kernelization
    D_kernelized = kernelize(D, num_kernels, kernel_width);

    % Inspect size: [10 × 5 × (100-20+1)]
    disp(size(D_kernelized));  % → [10  5  81]

---

## Installation

### Prerequisites
- MATLAB R2016a or later  
- No additional toolboxes required

### Setup

1. Save `kernelize.m` into a directory on your MATLAB path.
2. Add the folder:

        addpath('path/to/your/functions');

3. Verify availability:

        which kernelize

### Dependencies
- Built-in MATLAB functions: `conv`, `exp`, `linspace`, `max`, `warning`, basic matrix operations

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

- **v1.0 (2024-12-14)**: Initial release reviewed by Lovelace’s Square

---

## Keywords

- Kernelization  
- Convolution  
- Signal processing  
- MATLAB  
- Multivariate analysis  
- Three-way data  
-Trilinear
- Exponential kernels  
- Feature extraction  
- Trilinear decomposition  
- Data transformation

