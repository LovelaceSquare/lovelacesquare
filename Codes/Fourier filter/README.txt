# fourierFilter

Applies frequency-domain filtering to each row of a time-series matrix using FFT.

## Inputs
- `data` (matrix): signals arranged as rows.
- `freqIntervals` (Nx2 array): frequency bands [low high] in Hz to keep.
- `dt` (double, optional): sampling interval in seconds.

## Outputs
- `filteredData` (matrix): signals after filtering.
- `freqVector` (vector): frequencies corresponding to the FFT bins.

## Example
```matlab
intervals = [0 15; 40 60];
[Xfilt, f] = fourierFilter(X, intervals, 0.01);
```

## Author
Adrián Gómez-Sánchez

MIT License.
