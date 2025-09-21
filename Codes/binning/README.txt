# binning: N-Dimensional Binning by Sum or Mean (MATLAB)

Combines elements of an N-dimensional array into larger “bins” by summing or averaging.

**Reference**:  
Gómez-Sánchez, Adrián. (2025). *binning function*. Lovelace’s Square.  
https://lovelacesquare.org/

---

## Overview

The `binning` function reduces the resolution of an N-dimensional array by grouping adjacent elements along each dimension, as specified by the user’s `binVector`.

For each dimension:
- The function permutes that dimension to the front.
- It reshapes the data into blocks of size `binVector(dim)`.
- Then it either sums or averages those blocks according to the chosen `mode` (`'sum'` or `'mean'`).

The implementation is optimized to avoid full-array permutations. Only the current dimension is permuted into focus and then restored, making it efficient even for large N-D arrays.

If the array size is not an exact multiple of the bin sizes, the function truncates the data to the largest compatible size before binning.

---

## Inputs

- `X` (array): Input N-dimensional array to be binned  
- `binVector` (vector): Bin size for each dimension (must divide or truncate dimensions of `X`)  
- `mode` (string): Either `'sum'` or `'mean'`  

## Output

- `Xbinned`: The binned version of the input array, with reduced resolution along specified dimensions  

---

## Usage Example

(Insert into MATLAB script or command window)

    % Example 1: 2D sum binning
    X = rand(100, 200);            % 100×200 data matrix
    binVec = [4, 5];               % group every 4 rows and 5 columns
    mode   = 'sum';                % sum within each bin
    Xbinned = binning(X, binVec, mode);
    size(Xbinned)                  % returns [25 40]

    % Example 2: 3D mean binning
    Y = rand(30, 40, 50);          % 3D data
    binVec = [3, 4, 5];            % bin sizes along dimensions 1, 2, and 3
    mode   = 'mean';               % average within each bin
    Ybinned = binning(Y, binVec, mode);
    size(Ybinned)                  % returns [10 10 10]

---

## Installation

### Prerequisites
- MATLAB R2016a or later (uses built-in `permute`, `reshape`, `sum`, `mean`)

### Setup

1. Save `binning.m` into a folder on your MATLAB path.
2. Add the folder:

        addpath('path/to/your/functions');

3. Verify the function is found:

        which binning

### Dependencies
- No additional toolboxes required.

---

## License

Released under the **MIT License**.

---

## Authors

- **Adrián Gómez-Sánchez**  
- **Created**: April 4, 2025  

---

## Changelog

- **v2.0 (2025-04-04)**: Optimized selective-dimension permutation for faster binning; preserves truncation logic  
- **v1.0**: Initial implementation using full-dimension rotations

---

## Keywords

- binning  
- downsampling  
- MATLAB  
- reshape  
- permute  
- sum  
- mean  
- data reduction

