%% Test Script for EMSC Function
% Demonstrates the three reference-selection modes of EMSC
%   • 'Mean'    – dataset mean spectrum
%   • 'Median'  – dataset median spectrum
%   • 'External'– user-supplied reference spectrum
% A quadratic baseline is injected so polyOrder = 2 is required.
%
% Author:  Adrián Gómez-Sánchez
% Date Created: 2024-12-14
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v1.4  (2025-07-14)

%% 1) Parameters for Simulated Spectral Data
rng default                                    % reproducible
nSamples      = 20;
nWavelengths  = 1000;
wlMin = 800;   wlMax = 2500;
wl     = linspace(wlMin, wlMax, nWavelengths);
wlN    = linspace(0, 1,   nWavelengths);       % 0–1 axis for baseline

%% 2) “True” Spectra (Gaussian peaks, no baseline)
peakPos = [950 1200 1400 1900 2100 2300];
peakAmp = [1   0.8  1.2  0.9  1.1  0.7];
peakWid = [50  60   55   70   65   60];

trueSpectra = zeros(nSamples, nWavelengths);
for i = 1:nSamples
    for p = 1:numel(peakPos)
        A = peakAmp(p)*(1 + 0.05*randn);       % ±5 % amplitude jitter
        trueSpectra(i,:) = trueSpectra(i,:) + ...
            A*exp(-0.5*((wl-peakPos(p))/peakWid(p)).^2);
    end
end

%% 3) Inject Quadratic Baseline  (curvature only)
quadCoeff = 0.6*randn(nSamples,1);             % random curvature
baseline  = quadCoeff .* wlN.^2;               % N×λ via implicit expansion

%% 4) Multiplicative & Additive Scatter
mult = 0.85 + 2*rand(nSamples,1);              % 0.85–2.85
add  = 3*rand(nSamples,1);                     % 0–3

rawSpectra = (trueSpectra + baseline).*mult + add;

%% 5) Apply EMSC  (quadratic baseline → polyOrder = 2)
polyOrder = 2;

[correctedMean,   refMean]   = EMSC(rawSpectra,'Mean',    polyOrder);
[correctedMedian, refMedian] = EMSC(rawSpectra,'Median',  polyOrder);
customRef = trueSpectra(7,:);                              % external “gold”
[correctedExt,   refExt]    = EMSC(rawSpectra,'External', polyOrder, customRef);

%% 6) Mean-Squared Error versus Ground Truth
mseMean   = mean((correctedMean   - trueSpectra).^2,'all');
mseMedian = mean((correctedMedian - trueSpectra).^2,'all');
mseExt    = mean((correctedExt    - trueSpectra).^2,'all');
fprintf('MSE  Mean   reference : %.4e\n', mseMean);
fprintf('MSE  Median reference : %.4e\n', mseMedian);
fprintf('MSE  External reference: %.4e\n', mseExt);

%% 7) Visualisation
% 7-A  Global view
figure;
subplot(4,1,1);
plot(wl, rawSpectra','LineWidth',1);
title('Raw Spectra (peaks + quadratic baseline)'); grid on; xlim([wlMin wlMax]);

subplot(4,1,2);
plot(wl, correctedMean','LineWidth',1);
title('EMSC – Mean Reference (polyOrder = 2)'); grid on; xlim([wlMin wlMax]);

subplot(4,1,3);
plot(wl, correctedMedian','LineWidth',1);
title('EMSC – Median Reference (polyOrder = 2)'); grid on; xlim([wlMin wlMax]);

subplot(4,1,4);
plot(wl, correctedExt','LineWidth',1);
title('EMSC – External Reference (polyOrder = 2)'); grid on;
xlabel('Wavelength (nm)'); xlim([wlMin wlMax]);

% 7-B  Overlay one representative spectrum
s = randi(nSamples);
figure;
plot(wl, trueSpectra(s,:), 'k-', 'LineWidth',2); hold on;
plot(wl, rawSpectra(s,:),  'c--','LineWidth',1.5);
plot(wl, correctedMean(s,:),'r-.','LineWidth',1.4);
plot(wl, correctedMedian(s,:),'m:','LineWidth',1.4);
plot(wl, correctedExt(s,:),'g-','LineWidth',1.4); hold off;
legend({'True','Raw','Mean ref','Median ref','External ref'},'Location','Best');
title(sprintf('Sample #%d – EMSC Variants (quadratic baseline removed)', s));
xlabel('Wavelength (nm)'); ylabel('Intensity (a.u.)'); grid on; xlim([wlMin wlMax]);
