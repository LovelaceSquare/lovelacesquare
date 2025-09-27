# OALS: Orthogonalized Alternating Least Squares for PCA with Missing Data (MATLAB)

MATLAB function implementing Orthogonalized Alternating Least Squares (O-ALS) for PCA missing value imputation and data decomposition.

**Reference**:  
Gómez-Sánchez, Adrián, et al. "Solving the missing value problem in PCA by Orthogonalized-Alternating Least Squares (O-ALS)."  
*Chemometrics and Intelligent Laboratory Systems* (2024): 105153.

---

## Overview

This function addresses the common issue of missing values in principal component analysis (PCA) by iteratively estimating scores and loadings while enforcing orthogonality and normalization constraints.

The Orthogonalized Alternating Least Squares (OALS) algorithm proceeds as follows:

1. **Initialization**  
   - An initial loadings matrix is provided by the user or generated randomly.  
   - The first component is seeded with the column-wise means of the data to capture the dominant trend.

2. **Orthogonalization**  
   - Both loadings and scores are orthogonalized via Gram–Schmidt projection to ensure each component captures unique variance.

3. **Alternating Updates**  
   - **Scores Update**: Given current loadings, compute scores by least squares (`ScoresLS.m`)  
   - **Loadings Update**: Given updated scores, compute loadings by least squares (`LoadingsLS.m`)

4. **Normalization & Rescaling**  
   - Loadings vectors are normalized to unit norm  
   - Scores are rescaled accordingly to preserve explained variance

5. **Convergence Assessment**  
   - At each iteration, explained variance (r²) and lack of fit (LOF) are computed (`lofNaN.m`)  
   - Algorithm halts when relative LOF improvement < 1e-12 or max iterations reached

This method is well-suited for chemometrics, spectroscopy, or any domain as long as the dataset is bilinear and the number of component is well approximated.

---

## Inputs

- `D`: Data matrix with NaNs `[samples × variables]`  
- `maxIter`: Maximum number of iterations  
- `nComp`: Number of principal components  
- `P_init` *(optional)*: Initial loadings matrix `[nComp × variables]`

## Outputs

- `Dr`: Reconstructed/imputed data matrix  
- `T`: Scores matrix `[samples × nComp]`  
- `P`: Loadings matrix `[nComp × variables]`  
- `r2`: Explained variance per iteration  
- `lofc`: Lack of fit per iteration

---

## Usage Example

Paste into MATLAB:

    % Example 1: Random initialization with 3 components, 100 iterations
    [Dr, T, P, r2, lofc] = OALS(D, 100, 3);

    % Example 2: User-provided initial loadings for 2 components
    P_init = rand(2, size(D,2));
    [Dr, T, P, r2, lofc] = OALS(D, 100, 2, P_init);

---

## Installation

### Prerequisites
- MATLAB R2021a or later

### Setup

1. Save `OALS.m`, `ScoresLS.m`, `LoadingsLS.m`, and `lofNaN.m` in the same directory.
2. Add the folder to your MATLAB path:

        addpath(genpath('path/to/OALS_folder'))

### Dependencies

- `ScoresLS.m`: Computes scores from data and loadings  
- `LoadingsLS.m`: Computes loadings from data and scores  
- `lofNaN.m`: Calculates explained variance and lack of fit, handling NaNs

---

## License

Released under the **MIT License**

---

## Authors

- **Adrián Gómez-Sánchez**  
- **Date of Creation**: January 2, 2025  
- **Reviewed by**: Lovelace’s Square

---

## Changelog

- **v2.0 (2025-01-02)**:  
  Added support for user-provided initial loadings, multiple random initializations with best-LOF selection, and orthogonalization of scores.

---

## Keywords

- Orthogonalized Alternating Least Squares  
- OALS  
- Principal Component Analysis  
- PCA  
- Missing value imputation  
- Data decomposition  
- MATLAB  
- Chemometrics

