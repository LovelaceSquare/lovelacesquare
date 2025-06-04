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

## Output
- `dataBinned` (array): binned array with reduced size.

## Example
```matlab
B = binning(rand(100,200),[4 5],'sum');
```

## Author
Adrián Gómez-Sánchez

MIT License.
