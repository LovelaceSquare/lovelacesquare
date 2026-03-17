# AsLS GUI

Interactive MATLAB GUI for baseline correction using **Asymmetric Least Squares Smoothing** (AsLS). Estimates and subtracts baselines from spectral datasets with real-time preview.

Uses the **AppBase + uihtml** architecture: HTML/CSS/JS frontend with MATLAB backend.

## Features

- **Two-Parameter Model**: Lambda (smoothness) and p (asymmetry) with logarithmic sliders
- **Real-Time Preview**: Baseline preview updates live as you drag sliders
- **Signal-by-Signal Mode**: Preview on individual spectra or the mean spectrum
- **Show All Spectra**: Overlay all spectra with optional subsampling for large datasets
- **Batch Correction**: Apply baseline correction to every spectrum in the dataset
- **Progress Feedback**: Step-by-step progress bar for loading, correction, and export
- **Big Dataset Guard**: Warning and subsampling options for datasets > 1000 samples
- **Zoom/Pan**: Mouse wheel zoom + drag pan on charts
- **Export**: Corrected data and baselines exported as MATLAB workspace variables
- **Dark Mode**: Toggle via Ctrl+D or the + menu
- **Resizable Panels**: Drag handle between controls and charts
- **Demo Data**: Built-in synthetic spectra for quick exploration

## Installation

Add the `AsLS` folder to your MATLAB path:

```matlab
addpath('path/to/AsLS_GUI');
```

## Quick Start

```matlab
% Empty (load data from GUI)
AsLS_GUI();

% With a numeric matrix from workspace
AsLS_GUI(spectra);

% With struct (data + optional wavelength)
s.data = spectra;
s.wavelength = wavelength;
AsLS_GUI(s);
```

## Input Data

### Spectral Matrix
- 2D numeric matrix: samples (rows) x channels (columns)
- Must be real, finite (no NaN or Inf), and have at least 4 columns
- Supported types: `double`, `single`, integer types (converted to `double` internally)

### X-Axis Vector (optional)
- 1D numeric vector matching the number of columns
- Selectable from the Load Data dialog (e.g., wavelength, wavenumber, m/z)
- If not provided, channel indices (1, 2, ..., N) are used

## Algorithm

AsLS (Eilers, 2005) fits a smooth baseline **below** the signal using penalized least squares with asymmetric weights:

| Parameter | Symbol | Effect | Typical Range |
|-----------|--------|--------|---------------|
| Smoothness | &lambda; | Larger = smoother baseline | 10^3 &ndash; 10^8 |
| Asymmetry | p | Smaller = baseline pushed below peaks | 10^-4 &ndash; 10^-1 |
| Iterations | &mdash; | More = better convergence | 5 &ndash; 50 |

The preview shows the baseline on the **mean spectrum** (or a selected individual spectrum in signal-by-signal mode). Apply runs the correction on every spectrum independently.

## Architecture

```
AsLS/
├── AsLS_GUI.m                              # AppBase class (MATLAB backend)
├── AsLS_test.m                             # Demo data generator
├── ui/
│   └── asls_baseline_correction_ui.html    # HTML/CSS/JS frontend
├── business_logic/
│   ├── @AsLSCorrector/
│   │   └── AsLSCorrector.m                # AsLS solver (preview + batch correct)
│   └── @DataValidator/
│       └── DataValidator.m                 # Input validation
└── README.md
```

### Three-Layer Separation

| Layer | Files | Responsibility |
|-------|-------|----------------|
| **UI** | `ui/asls_baseline_correction_ui.html` | Layout, sliders, charts (Canvas 2D), modals, dark mode |
| **Backend** | `AsLS_GUI.m` | Action dispatcher, workspace I/O, orchestration |
| **Business Logic** | `@AsLSCorrector`, `@DataValidator` | Algorithm, validation (no GUI code) |

## Programmatic Usage (No GUI)

```matlab
% Create test data
AsLS_test;  % generates 'spectra' and 'wavelength' in workspace

% Use the corrector directly
corrector = AsLSCorrector();

% Preview baseline on mean spectrum
[meanSpec, baseline] = corrector.previewBaseline( ...
    spectra, 1e6, 0.001, struct.empty, false, 1e5, 10);

% Correct all spectra
[correctedData, baselineData] = corrector.correct( ...
    spectra, 1e6, 0.001, struct.empty, false, 1e5, 10);
```

## Requirements

- MATLAB R2020b or newer (for `uihtml`)
- No toolbox dependencies

## References

- Eilers, P. H. C. (2003). A Perfect Smoother. *Analytical Chemistry*, 75(14), 3631-3636.
- Eilers, P. H. C., & Boelens, H. F. M. (2005). Baseline Correction with Asymmetric Least Squares Smoothing.

## Author

Adrian Gomez-Sanchez

## License

MIT
