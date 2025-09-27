# autoscale: Z-Score Normalization of 2D Data (MATLAB)

Centers and scales a 2D data matrix to zero mean and unit variance by row or column in MATLAB.

**Reference**:  
Gómez-Sánchez, Adrián. (2021). *autoscale function*. Lovelace’s Square.  
https://lovelacesquare.org/

---

## Overview

The `autoscale` function standardizes each feature (column) or observation (row) of a 2D numeric matrix to have zero mean and unit variance. This z-score normalization step is essential in multivariate analyses (e.g., PCA, clustering, regression) to avoid bias from varying feature scales.

The function:
- Computes the mean and standard deviation across the selected dimension.
- Replaces any zero standard deviations with one (to avoid division by zero) and issues a warning.
- Applies element-wise centering and scaling.

The scaling direction can be specified as:
- `'column'` (default): each feature is scaled independently.
- `'row'`: each sample is scaled independently.

---

## Inputs

- `X` (matrix): 2D numeric data to be standardized  
- `dim` (optional, string): `'column'` (default) or `'row'`  

## Outputs

- `X_scaled`: Standardized version of `X`  
- `params`: Struct with fields:
  - `mean`: Mean values used for centering  
  - `std`: Standard deviations used for scaling  

---

## Usage Example

(Insert into MATLAB script or command window)

    % Example 1: Column-wise autoscaling (default)
    X = [1 2 3; 4 5 6; 7 8 9];
    [scaledCols, paramsCols] = autoscale(X);
    % Result:
    % scaledCols =
    %   [-1.2247  -1.2247  -1.2247
    %     0         0         0
    %     1.2247   1.2247   1.2247]
    % paramsCols.mean = [4 5 6]
    % paramsCols.std  = [2.4495 2.4495 2.4495]

    % Example 2: Row-wise autoscaling
    Y = [10 20 30; 5 5 5];
    [scaledRows, paramsRows] = autoscale(Y, 'row');
    % Result:
    % scaledRows(1,:) has zero mean and unit variance
    % scaledRows(2,:) is all zeros (zero variance row)
    % paramsRows.mean = [20; 5]
    % paramsRows.std  = [8.1649; 0]

---

## Installation

### Prerequisites
- MATLAB R2016b or later (due to implicit expansion)

### Setup

1. Save `autoscale.m` into a directory on your MATLAB path.
2. Add the folder:

        addpath('path/to/your/functions');

3. Confirm the function is available:

        which autoscale

### Dependencies
- None. Uses only core MATLAB functions.

---

## License

Released under the **MIT License**.

---

## Authors

- **Adrián Gómez-Sánchez**  
- **Created**: March 18, 2021  
- **Reviewed by**: Lovelace’s Square  

---

## Changelog

- **v2.0 (2021-03-18)**: Added default `'column'` behavior, row-wise scaling, and zero-variance handling with warning  
- **v1.0**: Initial release

---

## Keywords

- autoscaling  
- z-score normalization  
- standardization  
- MATLAB  
- preprocessing  
- zero mean  
- unit variance  
- data scaling  
- column scaling  
- row scaling

