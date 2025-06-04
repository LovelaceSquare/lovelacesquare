# CosmicPeakCorrection

Detects abrupt spikes in spectral data (e.g., cosmic rays) using derivatives and replaces affected channels via interpolation.

## Inputs
- `data` (matrix): spectral data.
- `derivativeOrder` (int): derivative order used for spike detection.
- `channelsToRemove` (int): number of neighboring channels to remove.
- `threshold` (double): derivative threshold for peak detection.

## Outputs
- `correctedData` (matrix): spectra with spikes removed.
- `peakMask` (matrix): logical mask of corrected points.

## Example
```matlab
[cleaned, mask] = CosmicPeakCorrection(X,1,2,5);
```

## Authors
Adrián Gómez-Sánchez and Rodrigo Rocha de Oliveira

MIT License.
