classdef CosmicPeakCorrector
% CosmicPeakCorrector. Detect and correct cosmic ray spikes in spectral data.
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-16
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 2.0
%
% Detects cosmic ray spikes using derivative-based thresholding and removes
% them via linear interpolation. Processing is applied per sample (row-wise).
%
% The algorithm computes abs(diff(data, order, 2)), flags channels exceeding
% the threshold, expands the mask by k channels on each side, sets flagged
% values to NaN, and fills gaps with fillmissing('linear').
%
% METHODS:
%   correct           - Full correction pipeline: detect + interpolate
%   computeDerivative - Compute absolute derivative of data
%   autoThreshold     - Compute automatic threshold from data statistics
%
% EXAMPLE:
%   corrector = CosmicPeakCorrector();
%   [corrected, mask] = corrector.correct(spectra, 1, 1, 100);
%   threshold = corrector.autoThreshold(spectra, 1);
%
% Disclaimer:
%   Authors and Lovelace's Square are not responsible for any issues,
%   inaccuracies, or data loss arising from the use of this software.

    methods
        function [corrected, correctionMask] = correct(~, data, derivativeOrder, channelsToRemove, threshold)
            % CORRECT Detect cosmic peaks and replace with interpolated values.
            %
            %   [corrected, correctionMask] = correct(obj, data, derivativeOrder, channelsToRemove, threshold)
            %
            %   Inputs:
            %       data             - [nSamples x nChannels] spectral matrix
            %       derivativeOrder  - Order of derivative (integer >= 1)
            %       channelsToRemove - Number of channels to expand on each side of peak
            %       threshold        - Detection threshold for absolute derivative
            %
            %   Outputs:
            %       corrected      - [nSamples x nChannels] corrected data matrix
            %       correctionMask - [nSamples x nChannels] logical mask of corrected positions

            % Compute derivative along columns (dimension 2)
            derivativeData = diff(data, derivativeOrder, 2);
            absDerivativeData = abs(derivativeData);

            % Detect peaks above threshold
            peakMask = absDerivativeData > threshold;
            [sampleIdx, peakLocs] = find(peakMask);

            % Map peak locations back to original data indices
            originalPeakChannels = peakLocs + derivativeOrder;

            % Filter out-of-bounds indices
            validIdx = originalPeakChannels <= size(data, 2);
            sampleIdx = sampleIdx(validIdx);
            originalPeakChannels = originalPeakChannels(validIdx);

            % Build correction mask with channel expansion
            correctionMask = false(size(data));
            for k = 1:length(originalPeakChannels)
                sample = sampleIdx(k);
                peakChannel = originalPeakChannels(k);
                startIdx = max(peakChannel - channelsToRemove, 1);
                endIdx = min(peakChannel + channelsToRemove, size(data, 2));
                correctionMask(sample, startIdx:endIdx) = true;
            end

            % Set masked values to NaN and interpolate
            corrected = data;
            corrected(correctionMask) = NaN;
            corrected = fillmissing(corrected, 'linear', 2, 'EndValues', 'nearest');
        end

        function [absDerivative, wavelengthDeriv] = computeDerivative(~, data, derivativeOrder, wavelength)
            % COMPUTEDERIVATIVE Compute the absolute derivative of spectral data.
            %
            %   [absDerivative, wavelengthDeriv] = computeDerivative(obj, data, derivativeOrder, wavelength)
            %
            %   Inputs:
            %       data            - [nSamples x nChannels] spectral matrix
            %       derivativeOrder - Order of derivative (integer >= 1)
            %       wavelength      - [1 x nChannels] wavelength vector
            %
            %   Outputs:
            %       absDerivative   - [nSamples x (nChannels - derivativeOrder)] absolute derivative
            %       wavelengthDeriv - [1 x (nChannels - derivativeOrder)] corresponding wavelengths

            derivativeData = diff(data, derivativeOrder, 2);
            absDerivative = abs(derivativeData);
            wavelengthDeriv = wavelength(derivativeOrder+1:end);
        end

        function threshold = autoThreshold(~, data, derivativeOrder)
            % AUTOTHRESHOLD Compute an automatic threshold based on data.
            %
            %   threshold = autoThreshold(obj, data, derivativeOrder)
            %
            %   Computes 10^floor(log10(max(abs(diff(data, order, 2))))) as a
            %   reasonable starting threshold for peak detection.
            %
            %   Inputs:
            %       data            - [nSamples x nChannels] spectral matrix
            %       derivativeOrder - Order of derivative (integer >= 1)
            %
            %   Outputs:
            %       threshold       - Suggested threshold value

            derivativeData = diff(data, derivativeOrder, 2);
            absMax = max(abs(derivativeData(:)));
            if absMax > 0
                threshold = 10^floor(log10(absMax));
            else
                threshold = 1;
            end
        end
    end
end
