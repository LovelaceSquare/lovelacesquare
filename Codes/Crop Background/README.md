# CropBackground

Crop background pixels from a 3D image cube using multiple selection methods.

Uses the **AppBase + uihtml** architecture.

## Features

- **Threshold mode** — min/max intensity sliders with real-time local preview; drag threshold lines directly on the distribution chart
- **Selection mode** — draw rectangle, polygon, or freehand regions on the heatmap to select pixels
- **Auto-Detect mode** — automatic particle/cell detection using Otsu or adaptive thresholding with morphological cleanup and optional watershed separation
- **Spectral band selection** — view the mean spectrum, draw ranges to select which wavelengths to use; right-click to remove a range
- **Preprocessing** — apply intensity transforms (square, log, square root) before any method
- Full-resolution heatmap with zoom (mouse wheel) and pan (click+drag)
- Sorted intensity distribution chart with zoom/pan and draggable threshold markers
- Jet colormap with colorbar
- Dark mode toggle (Ctrl+D)
- Export cropped matrix and pixel indices to workspace
- Toast notifications for all operations

## Installation

```matlab
addpath('path/to/CropBackground');
addpath(genpath('path/to/CropBackground/business_logic'));
```

## Quick Start

```matlab
% Generate test data (20 particles, 4 spectral types)
CropBackground_demo

% Launch the GUI
app = CropBackground(cube);

% Or launch empty and load from workspace
app = CropBackground();
```

## Input Data

| Field  | Size                       | Description                      |
|--------|----------------------------|----------------------------------|
| `cube` | `[rows x cols x channels]` | 3D numeric image cube (required) |

## Cropping Methods

### Threshold

1. Compute global intensity: `sum(imageCube, 3)`.
2. Optionally preprocess: square (x²), log(1+x), or sqrt(x).
3. Set min/max thresholds via sliders, text input, or by dragging the red lines on the distribution chart.
4. Pixels with intensity in `[min, max]` are retained; the rest are discarded.

### Selection (Drawing Tools)

| Tool      | Usage                                           |
|-----------|-------------------------------------------------|
| Rectangle | Click + drag to draw a bounding box              |
| Polygon   | Click to add vertices, double-click to close     |
| Freehand  | Click + drag to draw a free-form region          |

Pixels inside the shape are retained. Click **Invert** to keep the outside instead.

### Auto-Detect

Automatic object detection using the Image Processing Toolbox:

| Parameter      | Description                              | Default |
|---------------|------------------------------------------|---------|
| Method         | `otsu` (global) or `adaptive` (local)   | otsu    |
| Sensitivity    | Adaptive threshold sensitivity (0–1)     | 0.5     |
| Smoothing      | Gaussian sigma before thresholding       | 1.0     |
| Min area       | Minimum object size in pixels            | 50      |
| Watershed      | Separate touching particles              | off     |

Pipeline: normalize → smooth → threshold → morphological open/close → fill holes → remove small objects → clear border → (optional) watershed.

## Spectral Band Selection

Before applying any method, you can select which spectral bands to use for the intensity image:

1. Click **Select Bands** — the bottom chart switches to the mean spectrum
2. **Click + drag** on the spectrum to select a wavelength range (shaded blue)
3. Draw more ranges — they accumulate (union)
4. **Right-click** on a range to remove it
5. Click **Clear Ranges** to remove all selections
6. Click **All Bands** to reset and use the full spectrum

The intensity image is recomputed using only the selected bands: `sum(cube(:,:,selectedBands), 3)`. This is applied before preprocessing.

## Preprocessing Transforms

| Transform   | Formula       | Use case                              |
|-------------|---------------|---------------------------------------|
| None        | raw sum       | Default                               |
| Square      | x²            | Enhance contrast between bright/dim   |
| Log         | log(1+x)      | Compress large dynamic range          |
| Square root | sqrt(x)       | Mild compression                      |

## Architecture

```
CropBackground/
  CropBackground.m              Main GUI class (controller)
  CropBackground_demo.m         Test script (generates realistic demo data)
  ui/
    crop_background_ui.html     HTML/CSS/JS frontend (view)
  business_logic/
    @BackgroundCropper/
      BackgroundCropper.m       Intensity, thresholding, extraction, preprocessing
    @ParticleDetector/
      ParticleDetector.m        Automatic detection (Otsu/adaptive + morphology)
    @DataValidator/
      DataValidator.m           Input validation
```

**Three-Layer Separation**

| Layer      | File                        | Responsibility                                      |
|------------|-----------------------------|-----------------------------------------------------|
| View       | `crop_background_ui.html`   | HTML/CSS/JS frontend, heatmap, charts, drawing tools |
| Controller | `CropBackground.m`          | AppBase GUI, routes events between view and logic    |
| Model      | `BackgroundCropper.m`, `ParticleDetector.m`, `DataValidator.m` | Pure computation, no UI dependency |

## Programmatic Usage

```matlab
app = CropBackground();
s.cube = myImageCube;
app.setInputData(s);

% User interacts with GUI, then:
result = app.getData();
croppedMatrix = result.croppedMatrix;   % [nRetainedPixels x nChannels]
retainedIdx   = result.retainedIdx;     % Linear indices of kept pixels
discardedIdx  = result.discardedIdx;    % Linear indices of removed pixels
```

## Keyboard Shortcuts

| Shortcut     | Action                      |
|--------------|-----------------------------|
| Ctrl+D       | Toggle dark mode            |
| Escape       | Cancel current drawing      |
| Mouse wheel  | Zoom heatmap / distribution |
| Click+drag   | Pan heatmap / distribution  |

## Requirements

- MATLAB R2022a or later
- **Image Processing Toolbox** required only for Auto-Detect mode

## Author

Adrian Gomez-Sanchez

## License

MIT

## Changelog

### v1.1 (2026-03-29)
- Added Selection mode (rectangle, polygon, freehand drawing tools, multi-shape)
- Added Auto-Detect mode (Otsu/adaptive thresholding + watershed)
- Added spectral band selection (draw ranges on mean spectrum, right-click to remove)
- Added preprocessing transforms (square, log, sqrt)
- Heatmap fills container with zoom/pan support
- Draggable threshold lines on distribution chart
- Live local JS preview during threshold adjustment
- Toast notifications for all operations
- Demo with real plastic sample image + synthetic NIR spectra
- Removed orchestrator/Next Module pattern (standalone GUI)

### v1.0 (2026-03-16)
- Initial release with threshold-based cropping
