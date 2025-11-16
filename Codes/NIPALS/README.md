# NIPALS: Nonlinear Iterative Partial Least Squares for PCA with Missing Data (MATLAB)

MATLAB function implementing the NIPALS (Nonlinear Iterative Partial Least Squares) algorithm for PCA with missing value handling and data decomposition.

**Reference**:
Wold, H. "Nonlinear estimation by iterative least square procedures."
*Research Papers in Statistics: Festschrift for J. Neyman* (1966): 411-444.

---

## Overview

This function performs Principal Component Analysis (PCA) on data matrices containing missing values (NaNs) using the classic NIPALS algorithm. The method iteratively extracts principal components while automatically handling missing data.

The NIPALS algorithm proceeds as follows:

1. **Initialization**
   - For each component, the loading profile is initialized with the column means of the deflated data matrix.
   - The initialization automatically skips NaN values when computing means.

2. **Iterative Refinement**
   - **Scores Update**: Given current loadings, compute scores by least squares (`ScoresLS.m`)
   - **Loadings Update**: Given updated scores, compute loadings by least squares (`LoadingsLS.m`)

3. **Normalization & Rescaling**
   - Loadings vectors are normalized to unit norm
   - Scores are rescaled accordingly to preserve explained variance

4. **Convergence Assessment**
   - At each iteration, the change in loadings is monitored
   - Algorithm halts when change < 1e-12 or max iterations reached

5. **Deflation**
   - After each component is extracted, the data matrix is deflated by subtracting the contribution of that component (T_j * P_j)
   - Deflation is performed only on non-missing values

6. **Reconstruction**
   - The final reconstructed/imputed matrix is computed as Dr = T * P

This method is well-suited for chemometrics, spectroscopy, or any domain where missing data is common and a low-rank bilinear approximation is appropriate.

---

## Inputs

- `D`: Data matrix with NaNs `[samples × variables]`
- `it`: Maximum number of iterations per component
- `n`: Number of principal components to extract

## Outputs

- `Dr`: Reconstructed/imputed data matrix `[samples × variables]`
- `T`: Scores matrix `[samples × n]`
- `P`: Loadings matrix `[n × variables]`

---

## Usage Example

Paste into MATLAB:

```matlab
% Example 1: Extract 3 components with max 100 iterations per component
[Dr, T, P] = NIPALS(D, 100, 3);

% Example 2: Extract 5 components with max 200 iterations per component
[Dr, T, P] = NIPALS(D, 200, 5);

% Example 3: Reconstruct data using extracted components
D_reconstructed = T * P;
```

---

## Installation

### Prerequisites
- MATLAB R2016b or later (uses 'omitnan' option)

### Setup

1. Save `NIPALS.m`, `ScoresLS.m`, `LoadingsLS.m`, and `lofNaN.m` in the same directory.
2. Add the folder to your MATLAB path:

```matlab
addpath(genpath('path/to/NIPALS_folder'))
```

### Dependencies

- `ScoresLS.m`: Computes scores from data and loadings
- `LoadingsLS.m`: Computes loadings from data and scores
- `lofNaN.m`: Calculates explained variance and lack of fit, handling NaNs (optional, for quality assessment)

---

## Algorithm Details

### Missing Data Handling

The NIPALS algorithm naturally handles missing data (NaNs) by:
- Skipping NaN values when computing least squares estimates in both `ScoresLS.m` and `LoadingsLS.m`
- Using only available data points for each row/column during updates
- Performing deflation only on non-missing values

### Convergence

The convergence threshold is set to `1e-12` for the change in loadings between iterations. This can be adjusted in the main function if needed (line 76 of NIPALS.m).


## Limitations

- Performance degrades with high percentages of missing data (>50%)
- Assumes the data follows a low-rank bilinear structure
- In the presence of missing data, orthogonality of the extracted components is not guaranteed.
---

## License

Released under the **MIT License**

---

## Authors

- **Adrián Gómez-Sánchez**
- **Date of Creation**: January 2, 2025
- **Reviewed by**: Lovelace's Square

---

## Keywords

- NIPALS
- Nonlinear Iterative Partial Least Squares
- Principal Component Analysis
- PCA
- Missing value imputation
- Data decomposition
- MATLAB
- Chemometrics
