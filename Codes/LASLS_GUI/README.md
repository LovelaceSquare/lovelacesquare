# LAsLS GUI: Interactive Baseline Correction Interface

Interactive graphical interface for Local Asymmetric Least Squares (LAsLS) baseline correction. Built with MATLAB AppBase plus `uihtml` (HTML/JS frontend, MATLAB backend).

LAsLS extends Asymmetric Least Squares (AsLS) by allowing interval-specific asymmetry and smoothing parameters for baseline correction. The GUI provides real-time preview, interactive interval drawing, and batch processing for multi-signal datasets.

The GUI can also be used as a standard AsLS corrector, a Whittaker smoother, or a local Whittaker smoother, depending on the parameter configuration.

Reference: Eilers, P. H. C., and H. F. M. Boelens. Baseline correction with asymmetric least squares smoothing. Leiden University Medical Centre Report 1(1), 2005.

---

## Quick Start

Launch from MATLAB:

```matlab
LASLS
```

1. Load data. Click `Load` to import a matrix from the workspace (rows = samples, columns = channels). Alternatively, use `Generate demo data` from the `+` menu.
2. Draw intervals. Click `Draw` and drag on the preview chart to define peak regions. Each interval gets its own local smoothing and asymmetry parameters.
3. Tune parameters. Use the `Global` tab for baseline smoothness, asymmetry, and tension. Use the `Interval` tab to fine-tune individual intervals.
4. Preview. The baseline estimate updates in real time as parameters change. The corrected signal is shown in the bottom chart.
5. Apply and export. Click `Apply` to compute baselines for all samples, then `Export` to save corrected data and baselines to the workspace. Export names are validated, duplicate names are rejected, and existing workspace variables require overwrite confirmation.

---

## Features

- Real-time baseline preview
- Interactive interval drawing
- Automatic peak detection
- Signal-by-signal navigation
- Per-signal parameter mode
- Batch processing
- Import and export of interval definitions
- Safer workspace export with live overwrite warnings and confirmation
- Custom x-axis support
- Editable axis labels
- Interactive zoom and pan
- Dark mode
- Built-in tutorial

---

## Contents

| File/Folder | Description |
|-------------|-------------|
| `LASLS.m` | Main app class (MATLAB AppBase) |
| `LASLS_test.m` | Synthetic data generator for testing |
| `business_logic/@LASLSCorrector/` | Core IRLS baseline correction algorithm |
| `business_logic/@DataValidator/` | Input validation |
| `ui/lasls_baseline_correction_ui.html` | HTML/JS frontend |

---

## Installation

Add with subfolders to the MATLAB path:

```matlab
addpath(genpath('path/to/LASLS_GUI'))
```

### Dependencies

Built-in MATLAB functions only. No external toolboxes are required.

---

## License

Released under the MIT License.

---

## Authors

- Adrian Gomez-Sanchez -- Universitat de Barcelona and Lovelace Square, Barcelona, Spain
- Berta Torres-Cobos -- University of Copenhagen, Denmark and Lovelace Square, Barcelona, Spain
- Rodrigo Rocha de Oliveira -- Universitat de Barcelona and Lovelace Square, Barcelona, Spain
