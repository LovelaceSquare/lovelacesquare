# LAsLS GUI: Interactive Baseline Correction Interface

Interactive graphical interface for Local Asymmetric Least Squares (LAsLS) baseline correction. Built with MATLAB AppBase + uihtml (HTML/JS frontend, MATLAB backend).

LAsLS extends Asymmetric Least Squares (AsLS) by allowing interval-specific asymmetry and smoothing parameters for baseline correction. The GUI provides real-time preview, interactive interval drawing, and batch processing for multi-signal datasets.

The GUI can also be used as a standard AsLS corrector, a Whittaker smoother, or a local Whittaker smoother, depending on the parameter configuration.

**Reference**: Eilers, P. H. C., & Boelens, H. F. M. (2005). *Baseline correction with asymmetric least squares smoothing*. Leiden University Medical Centre Report, 1(1), 5.

---

## Quick Start

Launch from MATLAB:

```matlab
LASLS
```

1. **Load data** -- Click *Load* to import a matrix from the workspace (rows = samples, columns = channels). Alternatively, use *Generate demo data* from the *+* menu.
2. **Draw intervals** -- Click *Draw* and drag on the preview chart to define peak regions. Each interval gets its own local smoothing and asymmetry parameters.
3. **Tune parameters** -- Use the *Global* tab for baseline smoothness, asymmetry, and tension. Use the *Interval* tab to fine-tune individual intervals.
4. **Preview** -- The baseline estimate updates in real time as parameters change. The corrected signal is shown in the bottom chart.
5. **Apply & Export** -- Click *Apply* to compute baselines for all samples, then *Export* to save corrected data and baselines to the workspace.

---

## Features

- **Real-time baseline preview** -- Baseline and corrected signal update instantly as parameters change.
- **Interactive interval drawing** -- Click and drag to define peak regions directly on the chart.
- **Automatic peak detection** -- Built-in peak finder to assist with interval definition.
- **Signal-by-signal navigation** -- Browse individual samples with prev/next controls and editable sample index.
- **Per-signal parameter mode** -- Assign different global parameters to each sample for heterogeneous datasets.
- **Batch processing** -- Apply correction to all samples at once.
- **Import/Export intervals** -- Save and reload interval definitions across sessions.
- **Custom x-axis** -- Load a wavelength or retention time vector for meaningful axis labels.
- **Editable axis labels** -- Double-click axis labels to customize them.
- **Interactive zoom & pan** -- Scroll to zoom, drag to pan on both charts.
- **Dark mode** -- Toggle via the *+* menu.
- **Built-in tutorial** -- Step-by-step guided tour of all interface features.

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

Add with subfolders to your MATLAB path:

```matlab
addpath(genpath('path/to/LASLS_GUI'));
```

### Dependencies

Built-in MATLAB functions only: `spdiags`, matrix division (`\`), `norm`. No external toolboxes required.

---

## License

Released under the **MIT License**.

---

## Authors

- **Adrian Gomez-Sanchez** -- Universitat de Barcelona & Lovelace's Square, Barcelona, Spain
- **Berta Torres-Cobos** -- University of Copenhagen, Denmark & Lovelace's Square, Barcelona, Spain
- **Rodrigo Rocha de Oliveira** -- Universitat de Barcelona & Lovelace's Square, Barcelona, Spain
