# whittakerSmootherImpute: Iterative Smoothing with NaN Imputation (MATLAB)

Applies iterative Whittaker smoothing with NaN imputation to 2D data matrices in MATLAB, based on the method described by Edmund T. Whittaker in *“On a new method of graduation”*,  
*Proceedings of the Edinburgh Mathematical Society*, **41** (1922): 63–75.

---

## Overview

The `whittakerSmootherImpute` function addresses the common problem of smoothing noisy one-dimensional signals (rows of a matrix) while simultaneously handling missing data (`NaN`s).  
Traditional smoothing methods require complete data or pre-imputed values, which can introduce bias or artifacts.

This function combines:

1. **Initial Imputation**:  
   Linearly interpolates missing values in each row to provide a starting estimate.

2. **Whittaker Smoothing**:  
   Implements the penalized least squares approach of Whittaker (1922), minimizing:

       ||y − z||^2 + lambda * ||D^d * z||^2

   where:
   - `y`: observed (or imputed) signal  
   - `z`: smooth estimate  
   - `lambda`: smoothing parameter controlling fidelity vs. smoothness  
   - `D^d`: d-th order finite-difference matrix  

3. **Iterative Re-Imputation**:  
   After smoothing, replaces only the originally missing entries with their smoothed values and repeats until convergence (change < tolerance) or maximum iterations reached.

By alternating between smoothing and selective imputation, the algorithm converges to a smooth and self-consistent signal that respects both observed and missing-value structure.

---

## Inputs

- `data` (matrix): Input matrix containing NaNs  
- `lambda` (scalar): Smoothing parameter  
- `diffOrder` (integer): Finite difference order (`d`)  
- `maxIterations` (integer): Maximum number of imputation iterations  
- `tolerance` (scalar): Convergence threshold

## Outputs

- `smoothedMatrix`: Final smoothed matrix  
- `imputedMatrix`: Final version of `data` with NaNs replaced

---

## Usage Example

Paste into MATLAB:

```matlab
% Define your data matrix (rows = separate signals)
data = [
    1, 2, 3, NaN, 5;
    2, NaN, 4, 5, 6
];

% Set parameters
lambdaValue   = 10;    % Larger => smoother
diffOrder     = 2;     % Second‐order differences (curvature penalty)
maxIterations = 5;     % Maximum iterations for imputation
tolerance     = 1e-4;  % Convergence threshold

% Run smoother + imputer
[smoothedMatrix, imputedMatrix] = whittakerSmootherImpute( ...
    data, lambdaValue, diffOrder, maxIterations, tolerance);

% Inspect results
disp('Smoothed Output:');
disp(smoothedMatrix);

disp('Final Imputed Matrix:');
disp(imputedMatrix);
```

---

## Installation

### Prerequisites
- MATLAB R2016a or later (no additional toolboxes required)

### Setup

1. Save `whittakerSmootherImpute.m` into a folder on your MATLAB path.
2. Add the directory:

    ```matlab
    addpath('path/to/your/functions');
    ```

3. Verify availability:

    ```matlab
    which whittakerSmootherImpute
    ```

### Dependencies
- Only built-in MATLAB functions

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

- **v1.0 (2024-12-14)**:  
  Initial release and review by Lovelace’s Square.

---

## Keywords

- Whittaker  
- smoothing  
- iterative imputation  
- missing data  
- signal proces

