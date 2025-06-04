# binning

Group adjacent elements of an N‑D array and sum or average them according to the
specified bin sizes. The routine partitions the array into blocks and reduces
each block to a single value, effectively decreasing resolution while preserving
overall trends. This is commonly performed on hyperspectral images prior to
analysis to improve signal‑to‑noise ratios.

## Inputs
- `data` (array): array to bin.
- `binVector` (vector): bin size for each dimension.
- `mode` ("sum" or "mean"): how to combine elements.

## Outputs
- `dataBinned` (array): binned array with reduced size.

## Example
```matlab
B = binning(rand(100,200),[4 5],'sum');
```

## What the code does
Group adjacent elements of an N‑D array and sum or average them according to the

## How to use it
Run `binning.m` in MATLAB. A basic demonstration is provided in `test.m`.

## Installation/setup instructions
MATLAB R2018b or later. Add this folder to your MATLAB path.

## Usage examples
See `test.m` for a usage example.

## Contact information
contact@lovelacesquare.org

## Authors
Adrián Gómez-Sánchez

## License
MIT
