# normMatrix: Matrix Normalization Utility (MATLAB)

Normalizes an input matrix using specified norms and dimensions.

**Reference**:  
Gómez-Sánchez, Adrián. (2024). *normMatrix* function.  
Lovelace’s Square. https://lovelacesquare.org/

---

## Overview

The `normMatrix` function provides flexible normalization of numerical matrices in MATLAB by supporting a variety of common vector and matrix norms.

It accepts an input matrix and applies one of six normalization types:
- `'max'`: maximum value  
- `'euclidean'`: Euclidean norm  
- `'l1'`: L1 norm  
- `'l2'`: L2 norm  
- `'linf'`: L-infinity norm  
- `'frobenius'`: Frobenius norm  

Normalization can be applied along:
- Each row  
- Each column  
- The entire matrix (`'all'`)  

Internally:
- Norms are computed depending on type and dimension (e.g., `max(data(:))`, `norm(data,2)`, `sqrt(sum(data.^2,2))`, etc.)
- A warning is issued if any denominator is zero (to catch potential NaNs)
- Normalization is done using `bsxfun` for dimension-wise scaling

This utility is suitable for:
- Preprocessing for machine learning or chemometric models  
- Scaling features for comparability  
- Standardizing signals or measurements  

---

## Inputs

- `data` (matrix): Input numeric array `[M × N]`  
- `normType` (string): One of `'max'`, `'euclidean'`, `'l1'`, `'l2'`, `'linf'`, `'frobenius'`  
- `dim` (string): `'row'`, `'column'`, or `'all'`  

## Outputs

- `normalized`: Output matrix after normalization (same size as input)

---

## Usage Example

Paste into MATLAB:

    % Example: Normalize a random matrix row-wise using Euclidean norm
    data = randn(5,4);                    % 5×4 matrix of random values
    normalizedRow = normMatrix(data, 'euclidean', 'row');

    % Example: Normalize the entire matrix by its Frobenius norm
    normalizedAll = normMatrix(data, 'frobenius', 'all');

    % Example: Normalize each column by its maximum absolute value
    normalizedCol = normMatrix(data, 'linf', 'column');

    % Display original and normalized first row
    disp('Original first row:');
    disp(data(1,:));
    disp('Row-normalized (Euclidean) first row:');
    disp(normalizedRow(1,:));

---

## Installation

### Prerequisites
- MATLAB R2016a or later

### Setup

1. Save `normMatrix.m` in a directory on your MATLAB path.
2. Add the directory:

        addpath('path/to/normMatrix');

3. Verify availability:

        which normMatrix

### Dependencies
- Built-in MATLAB functions only: `max`, `sum`, `norm`, `bsxfun`, etc.

---

## License

Released under the **MIT License**

---

## Authors

- **Adrián Gómez-Sánchez**  
- **Date of Creation**: December 14, 2024

---

## Changelog

- **v1.0 (2024-12-14)**:  
  Initial release with support for `'max'`, `'euclidean'`, `'l1'`, `'l2'`, `'linf'`, and `'frobenius'` norms across `'all'`, `'row'`, and `'column'` dimensions. Includes warnings for zero denominators.

---

## Keywords

- matrix normalization  
- max normalization  
- Euclidean norm  
- L1 norm  
- L2 norm  
- L-infinity norm  
- Frobenius norm  
- MATLAB  
- data preprocessing  
- feature scaling

