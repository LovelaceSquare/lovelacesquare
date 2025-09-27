# OALS: Recovery of Saturated Spectral Peaks (MATLAB)

Recovers saturated spectral signals in MATLAB by treating clipped values as missing (NaNs) and reconstructing them with Orthogonalized Alternating Least Squares (O-ALS).

**Reference**:  
Gómez-Sánchez, Adrián, et al. *Solving the missing value problem in PCA by Orthogonalized-Alternating Least Squares (O-ALS).* Chemometrics and Intelligent Laboratory Systems (2024): 105153.  

---

## Overview

The `OALS` function applies Orthogonalized Alternating Least Squares (O-ALS) to estimate a bilinear model of spectral data:

```
X ≈ T · P
```

where:  
- `T` = scores (sample contributions)  
- `P` = loadings (spectral profiles)  

In the case of detector saturation:  
- Clipped channels are replaced with `NaN`  
- OALS iteratively estimates `T` and `P` under orthogonality constraints  
- Full spectra are reconstructed as `X̂ = T·P`, providing estimates for the saturated regions  

This recovery is possible because spectra lie in a **low-dimensional subspace**: even if some values are missing (saturated), the incomplete point still belongs to that subspace, and missing values can be inferred from it.

---

## Inputs

- `D`: Data matrix `[samples × variables]` with NaNs at saturated entries  
- `iter`: Maximum number of O-ALS iterations  
- `nPC`: Number of principal components (rank of the bilinear model)  
- `P` *(optional)*: Initial guess for loadings `[nPC × variables]`  

---

## Outputs

- `Dr`: Reconstructed data matrix (saturated values recovered)  
- `T`: Scores matrix `[samples × nPC]`  
- `P`: Loadings matrix `[nPC × variables]`  
- `r2`: Explained variance per iteration  
- `lofc`: Lack of fit per iteration  

---

## Usage Example

```matlab
% Example: Recover saturated peaks in simulated spectra

% Step 1: Replace saturated (clipped) values by NaN
D_withNaNs = D_clipped;
D_withNaNs(D_clipped >= ADC_max) = NaN;

% Step 2: Apply OALS with 3 components and 200 iterations
[Dr, T, P, r2, lofc] = OALS(D_withNaNs, 200, 3);

% Step 3: Substitute only missing entries in the clipped data
D_recovered = D_clipped;
missingMask = isnan(D_withNaNs);
D_recovered(missingMask) = Dr(missingMask);
```

---

## Installation / Setup

### Prerequisites
- MATLAB R2016a or later  

### Setup

1. Save `OALS.m` and dependencies (`ScoresLS.m`, `LoadingsLS.m`, `lofNaN.m`) in a folder  
2. Add path and verify:

```matlab
addpath('path/to/your/functions');
```

---

## Practical Considerations

- **Bilinear model and multiple samples**: OALS requires several spectra; single-spectrum recovery is ill-posed.  
- **Coverage of saturated regions**: a variable must be observed unsaturated in at least one sample; if all spectra saturate, it cannot be recovered.  
- **Rank selection (nPC)**: must be chosen or estimated beforehand; too small or too large leads to bias.  
- **Selective peaks**: for a dominant, localized peak, restrict to that region and use `nPC=1` for stability.  
- **Initialization and local minima**: OALS is nonconvex; run multiple initializations and compare final LOF.  

---

## License

Released under the **MIT License**

---

## Authors

- **Adrián Gómez-Sánchez**  
- **Date of Creation**: September 27, 2025  
- **Reviewed by**: Lovelace’s Square  

---

## Changelog

- **v1.0 (2025-09-27)**: Creation

---

## Keywords

- O-ALS  
- PCA with missing data  
- Saturation recovery  
- Chemometrics  
- Low-rank reconstruction  

