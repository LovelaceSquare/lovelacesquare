# I_SVD: Iterative SVD-Based PCA Imputation (MATLAB)

Performs iterative SVD-based PCA imputation to handle missing values in data matrices by alternating between low-rank approximation and data update until convergence.

**Reference**:  
Roweis, S. (1997). *EM algorithms for PCA and SPCA*. In Advances in Neural Information Processing Systems, 10.

---

## Overview

The `I_SVD` function implements an iterative Singular Value Decomposition (SVD) algorithm for Principal Component Analysis (PCA) imputation of missing data. Given a matrix with NaN entries, it:

1. **Initializes** missing values with column means.  
2. **Performs** an economy-sized SVD to extract the top `nComp` principal components.  
3. **Reconstructs** the matrix from these components, yielding a low-rank approximation.  
4. **Updates** the missing entries in the original matrix with the reconstructed values.  
5. **Repeats** steps 2–4 until the Lack of Fit (LOF) change falls below tolerance `tol` or until `maxIter` iterations.

During each iteration, explained variance (`r2`) and lack of fit (`lofc`) metrics are computed to monitor convergence.

---

## Inputs

- `D` (matrix): Data matrix with missing values as `NaN`  
- `nComp` (integer): Number of principal components to use  
- `maxIter` (integer): Maximum number of iterations  
- `tol` (scalar): Tolerance for LOF convergence  

## Outputs

- `Dimp`: Imputed data matrix (same size as `D`)  
- `T`: Scores matrix from PCA  
- `P`: Loadings matrix from PCA  
- `r2`: Vector of explained variance per iteration  
- `lofc`: Vector of LOF values per iteration  

---

## Usage Example

(Insert into MATLAB script or command window)

    % test_I_SVD.m
    % Demonstrate I_SVD imputation on a matrix with missing entries

    % Generate synthetic data with missing values
    rng(0);
    D_true = randn(100, 50);
    D = D_true;
    missingRate = 0.1;
    mask = rand(size(D)) < missingRate;
    D(mask) = NaN;

    % Set I_SVD parameters
    nComp = 5;         % number of principal components
    maxIter = 100;     % maximum iterations
    tol = 1e-8;        % convergence tolerance

    % Run imputation
    [Dimp, T, P, r2, lofc] = I_SVD(D, nComp, maxIter, tol);

    % Plot convergence metrics
    figure;
    subplot(2,1,1);
    plot(r2, '-o'); xlabel('Iteration'); ylabel('Explained Variance (R^2)');
    title('Convergence of Explained Variance');
    subplot(2,1,2);
    plot(lofc, '-o'); xlabel('Iteration'); ylabel('Lack of Fit (LOF)');
    title('Convergence of Lack of Fit');

    % Compare original and imputed heatmaps for missing entries
    figure;
    imagesc(D); colorbar; title('Original Data with NaNs');
    figure;
    imagesc(Dimp); colorbar; title('Imputed Data');

---

## Installation

### Prerequisites
- MATLAB R2016a or later

### Setup

1. Save `I_SVD.m` and ensure `lofNaN.m` (or include the provided subfunction) is on your MATLAB path.
2. Add the folder:

        addpath('path/to/I_SVD');

3. Verify installation:

        which I_SVD

### Dependencies
- Built-in MATLAB functions: `svd`, `mean`, `isnan`, `nan`, array indexing  
- Optional: `lofNaN` script (if used separately)

---

## License

Released under the **MIT License**.

---

## Authors

- **Adrián Gómez-Sánchez**  
- GitHub: [@adriangomez](https://github.com/adriangomez)  
- **Created**: January 2, 2025

---

## Changelog

- **v1.0 (2025-01-02)**: Initial release

---

## Keywords

- PCA imputation  
- SVD  
- Missing data  
- Iterative algorithm  
- MATLAB  
- Expectation  
- Maximization  
- EM  
- Low-rank approximation  
- Data preprocessing

