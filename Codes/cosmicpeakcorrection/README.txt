# CosmicPeakCorrection

Detect abrupt spikes in spectral data (e.g., cosmic rays) using derivatives and
replace the affected channels by interpolation. Derivatives amplify sudden
intensity jumps so spikes are easily found with a threshold. After detection, a
small window around each peak is marked as invalid and filled in by linear
interpolation, restoring a smooth spectrum while leaving genuine peaks intact.

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
=======
## What the code does
COSMICPEAKCORRECTION  Removes cosmic spikes from spectral data by

## How to use it
Run `CosmicPeakCorrection.m` in MATLAB. A basic demonstration is provided in `test.m`.

## Installation/setup instructions
MATLAB R2018b or later. Add this folder to your MATLAB path.

## Usage examples
See `test.m` for a usage example.

## Contact information
contact@lovelacesquare.org

## Authors
% Authors: Adrián Gómez-Sánchez and Rodrigo Rocha de Oliveira

## License
MIT
