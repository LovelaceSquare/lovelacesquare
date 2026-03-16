% CosmicPeak_test. Generate synthetic Raman spectra with cosmic ray spikes.
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-16
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 2.0
%
% This script generates 40 synthetic Raman spectra (9 Gaussian peaks over
% a fluorescence background) and injects sharp cosmic ray spikes into ~30%
% of the samples. The resulting data can be used as demo input for the
% CosmicPeakCorrection_GUI GUI.
%
% OUTPUT (left in workspace):
%   spectra    - [40 x 800] double (samples x channels)
%   wavelength - [1 x 800] double (Raman shift, cm-1)
%
% EXAMPLE:
%   CosmicPeak_test;
%   s.data = spectra; s.wavelength = wavelength;
%   app = CosmicPeakCorrection_GUI(s);
%
% Disclaimer:
%   Authors and Lovelace's Square are not responsible for any issues,
%   inaccuracies, or data loss arising from the use of this script.

rng(42);

nSamples = 40;
nChannels = 800;
wavelength = linspace(200, 3500, nChannels);

gauss = @(x, mu, h, w) h .* exp(-((x - mu).^2) ./ (2*w.^2));

spectra = zeros(nSamples, nChannels);

for s = 1:nSamples
    % Fluorescence background (broad)
    bg = 500 + 200*randn() + gauss(wavelength, 1800, 400+50*randn(), 800);

    % Raman peaks (typical organic compound)
    peaks = zeros(1, nChannels);
    peakPos = [520, 780, 1000, 1200, 1450, 1600, 2200, 2900, 3100];
    peakH   = [80, 60, 120, 90, 150, 200, 50, 180, 100];
    peakW   = [15, 12, 20, 18, 25, 22, 30, 40, 35];

    for p = 1:length(peakPos)
        h = peakH(p) * (0.7 + 0.6*rand());
        peaks = peaks + gauss(wavelength, peakPos(p)+3*(rand()-0.5), h, peakW(p));
    end

    spectra(s, :) = bg + peaks + 5*randn(1, nChannels);
end

% Insert cosmic ray spikes: 2-5 spikes per ~30% of samples
nAffected = round(nSamples * 0.3);
affectedIdx = randperm(nSamples, nAffected);

for i = 1:nAffected
    s = affectedIdx(i);
    nSpikes = randi([2 5]);
    spikePos = randi([5, nChannels-4], 1, nSpikes);

    for j = 1:nSpikes
        spikeHeight = 2000 + 3000*rand();  % Very intense, narrow spike
        col = spikePos(j);
        spectra(s, col) = spectra(s, col) + spikeHeight;
        % Some spikes affect 1-2 adjacent channels
        if rand() > 0.5 && col < nChannels
            spectra(s, col+1) = spectra(s, col+1) + spikeHeight*0.3;
        end
    end
end

fprintf('Created: spectra (%dx%d), wavelength (1x%d)\n', ...
    size(spectra,1), size(spectra,2), length(wavelength));
fprintf('%d samples contain cosmic ray spikes.\n', nAffected);
fprintf('Run CosmicPeakCorrection_GUI to detect and remove spikes.\n');

clearvars -except spectra wavelength
