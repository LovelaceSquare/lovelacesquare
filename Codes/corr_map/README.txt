# corr_map

Display a heatmap of correlation coefficients between spectra. The function
computes the pairwise correlation matrix across rows and visualizes it using a
color-coded image with optional numeric labels. This quick view helps identify
groups of similar spectra or outliers in a dataset.
## Inputs
- `C` (matrix): data matrix with samples in rows.

## Outputs
- `figureHandle` (handle): handle to the correlation map figure.

## Example
```matlab
corr_map(C);
```
## What the code does
Display a heatmap of correlation coefficients between spectra. The function

## How to use it
Run `corr_map.m` in MATLAB. A basic demonstration is provided in `test_corr_map.m`.

## Installation/setup instructions
MATLAB R2018b or later. Add this folder to your MATLAB path.

## Usage examples
See `test_corr_map.m` for a usage example.

## Contact information
contact@lovelacesquare.org

## Authors
Adrián Gómez-Sánchez

## License
MIT
