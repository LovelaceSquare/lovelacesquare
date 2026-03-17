% AsLS_test. Generate synthetic spectra with baselines for testing.
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-16
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 2.0
%
% Creates spectra with polynomial baselines and Gaussian peaks. Use as
% demo data for the AsLS_GUI baseline correction tool.
%
% OUTPUT (left in workspace):
%   spectra    - [3000 x 500] double (samples x channels)
%   wavelength - [1 x 500] double
%
% EXAMPLE:
%   AsLS_test;
%   s.data = spectra; s.wavelength = wavelength;
%   app = AsLS_GUI(s);
%
% Disclaimer:
%   Author and Lovelace's Square are not responsible for any issues,
%   inaccuracies, or data loss arising from the use of this script.

rng(42);

nSamples = 3000;
nChannels = 500;
wavelength = linspace(400, 900, nChannels);
t = linspace(0, 1, nChannels);

spectra = zeros(nSamples, nChannels);

gauss = @(x, mu, h, w) h .* exp(-((x - mu).^2) ./ (2*w.^2));

for s = 1:nSamples
    % Polynomial baseline with per-sample variation
    a0 = 500 + 80*randn();
    a1 = 200 + 40*randn();
    a2 = 300 + 60*randn();
    a3 = -150 + 30*randn();
    baseline = a0 + a1*t + a2*t.^2 + a3*t.^3;

    % Exponential decay component
    baseline = baseline + (400 + 50*randn()) .* exp(-3*t);

    % 5 Gaussian peaks with variable intensities
    peaks = zeros(1, nChannels);
    peakPositions = [480, 550, 630, 720, 810];
    peakHeights   = [120, 250, 180, 300, 150];
    peakWidths    = [1, 2, 4, 1, 0.5];

    for p = 1:length(peakPositions)
        h = peakHeights(p) * (0.7 + 0.6*rand());
        pos = peakPositions(p) + 3*(rand()-0.5);
        peaks = peaks + gauss(wavelength, pos, h, peakWidths(p));
    end

    % Combine with Gaussian noise
    spectra(s, :) = baseline + peaks + 3*randn(1, nChannels);
end

clearvars -except spectra wavelength

fprintf('Created: spectra (%dx%d), wavelength (1x%d)\n', ...
    size(spectra,1), size(spectra,2), length(wavelength));
fprintf('Run AsLS_GUI to load and process this data.\n');
