%% LASLS_test - Generate synthetic chromatogram matrix for testing
% Creates a 20x4000 matrix with multiple synthetic chromatograms.
% Each sample has slightly different peak intensities and noise.
% Use as demo data for the LASLS GUI.
%
% Authors: Adrián Gómez-Sánchez, Berta Torres-Cobos, Rodrigo Rocha de Oliveira
% License: MIT
% Repository: https://github.com/LovelaceSquare/lovelacesquare

rng(42);

% Grid
n = 4000;
nSamples = 20;
x = linspace(50, 130, n)';

% Pre-allocate output
chromatogram = zeros(nSamples, n);

for s = 1:nSamples
    % Baseline with slight variations per sample
    t = x - x(1);
    T = x(end) - x(1);
    tn = t / T;

    % Add sample-specific variation to baseline
    baselineVariation = 1 + 0.1 * randn();

    solvent = (950 + 50*randn()) * exp(-t / (7 + 0.5*randn()));
    drift = (40 + 5*randn()) * tn + (120 + 15*randn()) * tn.^2;
    sig = 1 ./ (1 + exp(-(x - 95) / 1.7));
    late_drift = (900 + 80*randn()) * tn .* sig;
    wander = smoothdata(randn(size(x)), 'gaussian', 801);
    wander = (12 + 3*randn()) * wander / std(wander);
    baseline = (90 + 10*randn()) + solvent + drift + late_drift + wander;
    step_sig = 1 ./ (1 + exp(-(x - 88) / 0.08));
    baseline = baseline + (220 + 30*randn()) * step_sig;
    baseline = smoothdata(baseline, 'gaussian', 301);
    baseline = baseline * baselineVariation;

    % Peaks with sample-specific intensity variation
    gauss = @(xv, pos, h, w) h .* exp(-((xv - pos).^2) ./ (2*w.^2));
    peaks = zeros(size(x));

    % Major peaks with variable intensities
    major = [56.5 150 0.08; 65.2 200 0.09; 74.3 130 0.08; 83.8 170 0.10;
             95.2 145 0.11; 106.5 185 0.11; 115.8 135 0.10; 124.2 160 0.12];
    for i = 1:size(major, 1)
        % Vary height by 10-30% per sample
        heightVar = major(i,2) * (0.8 + 0.4*rand());
        % Slight position shift (up to 0.2 units)
        posVar = major(i,1) + 0.2*(rand()-0.5);
        peaks = peaks + gauss(x, posVar, heightVar, major(i,3));
    end

    % Minor peaks
    rng(456 + s*100);  % Different seed per sample for minor peaks
    positions = sort(51 + 78 * rand(41, 1));
    heights = 3 + 5 * lognrnd(2.5, 1.5, 41, 1);
    heights(heights < 200) = 0;
    heights = heights .* (0.7 + 0.6*rand(41, 1));  % Variable heights
    widths = max(0.02, 0.04 + 0.18 * (positions - 50) / 80 + 0.1 * rand(41, 1));

    for i = 1:41
        p = gauss(x, positions(i), heights(i), widths(i));
        tail = zeros(size(x));
        idx = x >= positions(i);
        tail(idx) = (heights(i) * 0.35) * exp(-(x(idx) - positions(i)) / 0.55);
        peaks = peaks + p + tail;
    end

    % Leading edge artifact
    peaks = peaks + gauss(x, 0, 260 + 40*randn(), 7.5);

    % Combine and add Poisson noise
    clean = baseline + peaks;
    offset = abs(min(clean)) + 100;
    chromatogram(s, :) = poissrnd(clean + offset) - offset;
end

% Clear all variables except chromatogram
clearvars -except chromatogram

fprintf('Created: chromatogram (%dx%d) - %d samples, %d channels\n', size(chromatogram,1), size(chromatogram,2), size(chromatogram,1), size(chromatogram,2));
fprintf('Run LASLS to load and process this data.\n');
