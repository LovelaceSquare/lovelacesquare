## Title & Short Description & Reference

**PURE**: Legacy weighted-SIMPLISMA algorithm for selecting the purest variables in multivariate data matrices.  
**Reference**: Windig, W., & Guilment, J. (1991). Interactive Self-Modeling Mixture Analysis. *Analytical Chemistry*, 63(14), 1425–1432.

## Overview

The `pure` function implements a legacy variant of the SIMPLISMA (SIMPle-to-use Interactive Self-modeling Mixture Analysis) algorithm for selecting "pure" variables—those that are likely dominated by single chemical species in multivariate data. Originally developed for spectroscopic applications, this method evaluates each variable's purity based on its standard deviation and mean, accounting for noise and collinearity.

This version reproduces the behavior of the classic `pure.m` script popularized in the 1990s, with the following elements:

- **Purity Score**:  
  p_j = sigma_j / (mu_j + n)  
  where sigma_j and mu_j are the standard deviation and mean of variable j, and n is a fixed noise offset.

- **Extra Legacy Weight**:  
  w_j = (sigma_j^2 + mu_j^2) / (sigma_j^2 + (mu_j + n)^2)

- **Weighted Purity**:  
  P1(j) = w_j * p_j

- **Collinearity Penalty**:  
  From step 2 onward, each candidate variable is evaluated by the determinant of the correlation matrix built from already-selected variables plus the candidate. Before correlation, each column is scaled by:  
  l_j = sqrt(sigma_j^2 + (mu_j + n)^2)

- **Selection Criterion**:  
  At each iteration k, the variable j that maximizes:  
  Pk(j) = P1(j) * det(Rk(j))  
  is selected, where Rk(j) is the correlation submatrix including candidate j and the k–1 previously selected variables.

The function returns the raw data profiles and indices of the selected variables. Profiles are normalised row-wise to unit maximum.

## Usage

```matlab
% Example data: rows = spectra, cols = variables
X = randn(100, 50) + linspace(0,10,50);  % simulated example

% Extract 3 pure components with 1% noise threshold
[pureProfiles, pureIndices] = pure(X, 3, 1);

% Visualise
figure, plot(pureProfiles'); title('PURE-selected profiles');
disp('Pure variable indices:'), disp(pureIndices);
```

## Installation / Setup

1. **Prerequisites**  
   - MATLAB R2016a or later.

2. **Setup**  
   - Place `pure.m` in a directory on your MATLAB path:  
     ```matlab
     addpath('path/to/pure');
     ```

3. **Dependencies**  
   - No external toolboxes required.

## License

This project is released under the **MIT License**.

## Author

Legacy SIMPLISMA pure.m implementations circulated from multiple sources; the most widely used open version appears in the Barcelona MCR-ALS toolboxes (Jaumot, Tauler, de Juan). Proprietary implementations exist in PLS_Toolbox (Eigenvector).

Adrián Gómez-Sánchez — polished and documented this legacy-style implementation.

Date of Creation: August 4, 2025

Reviewed by: Lovelace’s Square


## Changelog

- **v1.0 (2025-08-04)**: First public release, reproduces 1990s `pure.m` logic with additional documentation and plotting support.

## Keywords

SIMPLISMA, self-modeling mixture analysis, pure variable selection, spectral unmixing, correlation weighting
