# pcaFilter: PCA-Based Dimensionality Reduction and Reconstruction (MATLAB)

Performs PCA‐based dimensionality reduction and reconstruction of a matrix in MATLAB.

**Reference**:  
Gómez-Sánchez, Adrián. (2024). pcaFilter function. Lovelace’s Square. https://lovelacesquare.org/

---

## Overview

The `pcaFilter` function applies Principal Component Analysis (PCA) to an input data matrix to filter out lower‐variance components (often corresponding to noise) and reconstruct the data using only the top principal components.

Internally, it performs a singular value decomposition (SVD) of the input matrix:

    X = U·S·Vᵗ

It then retains the first *k* columns of `U` and `V` and the leading *k×k* block of `S`, reconstructing the matrix as:

    X̂ = Uₖ·Sₖ·Vₖᵗ

where *k* = `numComponents`.  
This produces a low‐rank approximation that preserves the dominant structure of the data while suppressing noise and irrelevant variation.

---

## Inputs

- `inputMatrix`: Numeric matrix `[samples × variables]`  
- `numComponents`: Number of principal components to retain

## Outputs

- `filteredMatrix`: Reconstructed matrix of same size as input

---

## Usage Example

Paste into MATLAB:

    % Example: Filter a 100×200 data matrix keeping 5 principal components
    inputMatrix   = randn(100,200) + peaks(100,200)*0.1;  % simulated noisy data
    numComponents = 5;

    % Run PCA filter
    filteredMatrix = pcaFilter(inputMatrix, numComponents);

    % Verify dimensions
    disp(size(filteredMatrix));  % should display [100 200]

    % Compare original vs. filtered (e.g., plot first row)
    figure;
    plot(inputMatrix(1,:), 'k--', 'DisplayName','Original');
    hold on;
    plot(filteredMatrix(1,:), 'b-', 'DisplayName','Filtered');
    legend;
    xlabel('Variable Index');
    ylabel('Value');
    title('Original vs. PCA-Filtered Signal');

---

## Installation / Setup

### Prerequisites
- MATLAB R2016a or later (no additional toolboxes required)

### Setup

1. Save `pcaFilter.m` to a folder on your MATLAB path
2. Add path and verify:

        addpath('path/to/your/functions');
        which pcaFilter

### Dependencies

- Built-in MATLAB functions only

---

## License

Released under the **MIT License**

---

## Authors

- **Adrián Gómez-Sánchez**  
- **Date of Creation**: December 14, 2024  
- **Reviewed by**: Lovelace’s Square

---

## Changelog

- **v1.0 (2024-12-14)**: Initial release reviewed by Lovelace’s Square

---

## Keywords

- PCA  
- Principal component analysis  
- Dimensionality reduction  
- MATLAB  
- Signal filtering  
- Denoising  
- SVD  
- Low-rank approximation  
- Eigenvectors  
- Eigenvalues


