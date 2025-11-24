# PCA-ALSQR: Principal Component Analysis using Alternating Least Squares with QR Orthogonalization (MATLAB)

MATLAB function implementing PCA-ALSQR (Alternating Least Squares with QR orthogonalization) for PCA with missing value handling and guaranteed orthogonality.

**Reference**:
Gómez-Sánchez, Adrián, et al. "Solving the missing value problem in PCA by Orthogonalized-Alternating Least Squares (O-ALS)."
*Chemometrics and Intelligent Laboratory Systems* (2024): 105153.

---

## Overview

This function performs Principal Component Analysis (PCA) on data matrices containing missing values (NaNs) using an efficient variant of the Orthogonalized-Alternating Least Squares (O-ALS) algorithm. The method achieves the same PCA subspace as O-ALS with reduced computational cost by enforcing orthogonality only once after convergence.

The PCA-ALSQR algorithm proceeds as follows:

1. **Initialization**
   - Loadings matrix P is initialized with the first component set to column means
   - Additional components are randomly initialized
   - Custom initialization can be provided via P_init parameter

2. **Alternating Least Squares Iterations**
   - **Scores Update**: Given current loadings, compute scores by least squares (`ScoresLS.m`)
   - **Loadings Update**: Given updated scores, compute loadings by least squares (`LoadingsLS.m`)
   - **Quality Assessment**: Compute explained variance (r²) and lack of fit using `lofNaN.m`

3. **Convergence Assessment**
   - Monitor relative change in lack of fit between iterations
   - Algorithm halts when change < 1e-12 or max iterations reached

4. **QR Orthogonalization**
   - After ALS convergence, perform QR decomposition on scores: T = QR
   - Transform loadings: M = R * P
   - This ensures orthogonality of the final components

5. **Optional SVD Refinement**
   - If useSVD=true (default), perform SVD on M = U * S * V'
   - Reorder components by descending singular values
   - Build final orthogonal scores: T = Q * U * sqrt(S)
   - Build final orthogonal loadings: P = sqrt(S) * V'

6. **Normalization**
   - Loadings vectors are normalized to unit norm
   - Scores are rescaled accordingly to preserve explained variance

7. **Reconstruction**
   - The final reconstructed/imputed matrix is computed as Dr = T * P

This method is well-suited for chemometrics, spectroscopy, or any domain where missing data is common and orthogonal components are required for interpretability and further analysis.

---

## Inputs

- `D`: Data matrix with NaNs `[samples × variables]`
- `iter`: Maximum number of ALS iterations
- `nPC`: Number of principal components to extract
- `P_init`: (Optional) Initial loadings matrix `[nPC × variables]`
- `useSVD`: (Optional) Logical flag (default=true). If true, performs SVD refinement after QR orthogonalization

## Outputs

- `Dr`: Reconstructed/imputed data matrix `[samples × variables]`
- `T`: Orthogonal scores matrix `[samples × nPC]`
- `P`: Orthogonal loadings matrix `[nPC × variables]`
- `r2`: Explained variance per iteration `[1 × #iterations]`
- `lofc`: Lack of fit per iteration `[1 × #iterations]`

---

## Usage Example

Paste into MATLAB:

```matlab
% Example 1: Extract 3 components with max 100 iterations (with SVD refinement)
[Dr, T, P, r2, lofc] = PCA_ALSQR(D, 100, 3);

% Example 2: Extract 5 components without SVD refinement
[Dr, T, P, r2, lofc] = PCA_ALSQR(D, 200, 5, [], false);

% Example 3: Use custom initialization for loadings
P_init = randn(3, size(D,2));  % Custom initialization for 3 components
[Dr, T, P, r2, lofc] = PCA_ALSQR(D, 100, 3, P_init);

% Example 4: Monitor convergence
[Dr, T, P, r2, lofc] = PCA_ALSQR(D, 100, 3);
figure;
subplot(2,1,1); plot(r2); title('Explained Variance vs Iteration');
subplot(2,1,2); plot(lofc); title('Lack of Fit vs Iteration');
```

---

## Installation

### Prerequisites
- MATLAB R2016b or later (uses 'omitnan' option)

### Setup

1. Save `PCA_ALSQR.m`, `ScoresLS.m`, `LoadingsLS.m`, and `lofNaN.m` in the appropriate directories:
   - `PCA_ALSQR.m` in the main folder
   - `ScoresLS.m`, `LoadingsLS.m`, `lofNaN.m` in the `dependencies` subfolder

2. Add the folder to your MATLAB path:

```matlab
addpath(genpath('path/to/PCA_ALSQR_folder'))
```

### Dependencies

- `ScoresLS.m`: Computes scores from data and loadings using least squares
- `LoadingsLS.m`: Computes loadings from data and scores using least squares
- `lofNaN.m`: Calculates explained variance (r²) and lack of fit, handling NaNs

---

## Algorithm Details

### Missing Data Handling

The PCA-ALSQR algorithm naturally handles missing data (NaNs) by:
- Skipping NaN values when computing least squares estimates in both `ScoresLS.m` and `LoadingsLS.m`
- Using only available data points for each row/column during updates
- Computing quality metrics (r², LOF) only on observed values

### Orthogonality

Unlike NIPALS, which does not guarantee orthogonality with missing data:
- PCA-ALSQR enforces strict orthogonality through QR decomposition
- Orthogonal components are essential for interpretation and downstream analysis
- The method achieves the same subspace as O-ALS but with lower computational cost

### Convergence

The convergence criterion is based on relative change in lack of fit:
```
conv = (lofc(it-1) - lofc(it)) / lofc(it)
```
Iterations stop when `abs(conv) < 1e-12` or when `iter` is reached.

### SVD Refinement

When `useSVD=true` (default):
- Components are reordered by descending variance (singular values)
- Provides additional numerical stability
- Ensures components are ordered by importance

When `useSVD=false`:
- Only QR orthogonalization is performed
- Slightly faster but components may not be variance-ordered

---

## Comparison with NIPALS

| Feature | NIPALS | PCA-ALSQR |
|---------|--------|-----------|
| Missing data | ✓ | ✓ |
| Orthogonality | ✗ (with NaNs) | ✓ |
| Computational cost | Lower | Moderate |
| Component ordering | Variance-ordered | Variance-ordered (with SVD) |
| Convergence monitoring | Per-component | Global |
| Best for | Quick exploration | Publication-ready analysis |

---

## Limitations

- Performance degrades with high percentages of missing data (>50%)
- Assumes the data follows a low-rank bilinear structure
- Slightly more computationally expensive than NIPALS due to orthogonalization step
- Convergence may be slower than NIPALS for some datasets

---

## License

Released under the **MIT License**

---

## Authors

- **Adrián Gómez-Sánchez**
- **Date of Creation**: March 1, 2025
- **Version**: v 1.0
- **Reviewed by**: Lovelace's Square

---

## Keywords

- PCA-ALSQR
- Alternating Least Squares
- QR Orthogonalization
- O-ALS
- Principal Component Analysis
- PCA
- Missing value imputation
- Data decomposition
- Orthogonal components
- MATLAB
- Chemometrics
