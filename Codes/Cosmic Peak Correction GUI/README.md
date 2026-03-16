# Cosmic Peak Correction

Interactive MATLAB GUI for detecting and removing **cosmic ray spikes** from spectral datasets using derivative-based detection and linear interpolation.

Uses the **AppBase + uihtml** architecture: HTML/CSS/JS frontend with MATLAB backend.

## Features

- **Derivative-Based Detection**: Configurable derivative order (1st, 2nd, ...) with absolute value thresholding
- **Automatic Threshold**: Data-driven starting threshold from `10^floor(log10(max(|derivative|)))`
- **Manual Threshold Tuning**: Logarithmic slider + click-to-set on the derivative chart
- **All Spectra Overlay**: View all individual spectra overlaid (with subsampling for large datasets)
- **Signal-by-Signal Mode**: Browse individual spectra with navigation badges
- **Corrected Toggle**: Show corrected only or both original + corrected
- **Derivative Chart**: Pointwise max |derivative| across all samples with sort mode (rank vs channel)
- **Channel Expansion**: Configurable number of channels to remove around each detected spike
- **Batch Correction**: Apply correction to all spectra simultaneously
- **Progress Feedback**: Step-by-step progress bar for loading, correction, and export
- **Zoom/Pan**: Mouse wheel zoom + drag pan on both charts
- **Export**: Corrected data exported as MATLAB workspace variable
- **Dark Mode**: Toggle via Ctrl+D or the + menu
- **Resizable Panels**: Drag handle between controls and charts
- **Demo Data**: Built-in Raman-like spectra with synthetic cosmic spikes

## Installation

Add the `CosmicPeakCorrection` folder to your MATLAB path:

```matlab
addpath('path/to/CosmicPeakCorrection');
```

## Quick Start

```matlab
% Empty (load data from GUI)
CosmicPeakCorrection();

% With a numeric matrix from workspace
CosmicPeakCorrection(spectra);

% With struct (data + optional wavelength)
s.data = spectra;
s.wavelength = wavelength;
CosmicPeakCorrection(s);
```

## Input Data

### Spectral Matrix
- 2D numeric matrix: samples (rows) x channels (columns)
- Must be real, finite (no NaN or Inf), and have at least 3 columns
- Supported types: `double`, `single`, integer types (converted to `double` internally)

### X-Axis Vector (optional)
- 1D numeric vector matching the number of columns
- Selectable from the Load Data dialog (e.g., wavelength, wavenumber, Raman shift)
- If not provided, channel indices (1, 2, ..., N) are used

## Algorithm

The correction follows five steps per sample:

1. Compute `abs(diff(spectrum, d))` where *d* is the derivative order.
2. Compare each value against the threshold.
3. Flag channels where |derivative| > threshold as spike positions.
4. Expand each flagged position by *k* channels on each side.
5. Set flagged channels to NaN, then `fillmissing('linear')` to interpolate.

| Parameter | Effect | Default | Typical Range |
|-----------|--------|---------|---------------|
| Derivative order | Higher orders amplify narrow spikes | 1 | 1 -- 3 |
| Channels to remove | Expansion radius around each spike | 1 | 0 -- 5 |
| Threshold | Lower = more aggressive detection | Auto | 10^-2 -- 10^6 |

## Architecture

```
CosmicPeakCorrection/
├── CosmicPeakCorrection.m                 # AppBase class (MATLAB backend)
├── CosmicPeak_test.m                      # Demo data generator
├── ui/
│   └── cosmic_peak_correction_ui.html     # HTML/CSS/JS frontend
├── business_logic/
│   ├── @CosmicPeakCorrector/
│   │   └── CosmicPeakCorrector.m          # Detection & interpolation logic
│   └── @DataValidator/
│       └── DataValidator.m                # Input validation
└── README.md
```

### Three-Layer Separation

| Layer | Files | Responsibility |
|-------|-------|----------------|
| **UI** | `ui/cosmic_peak_correction_ui.html` | Layout, sliders, charts (Canvas 2D), modals, dark mode |
| **Backend** | `CosmicPeakCorrection.m` | Action dispatcher, workspace I/O, orchestration |
| **Business Logic** | `@CosmicPeakCorrector`, `@DataValidator` | Algorithm, validation (no GUI code) |

## Programmatic Usage (No GUI)

```matlab
% Create test data
CosmicPeak_test;  % generates 'spectra' and 'wavelength' in workspace

% Use the corrector directly
corrector = CosmicPeakCorrector();

% Correct all spectra
[correctedData, correctionMask] = corrector.correct( ...
    spectra, 1, 1, 100);

% Compute automatic threshold
threshold = corrector.autoThreshold(spectra, 1);
```

## Requirements

- MATLAB R2020b or newer (for `uihtml`)
- No toolbox dependencies

## Author

Adrian Gomez-Sanchez

## License

MIT
