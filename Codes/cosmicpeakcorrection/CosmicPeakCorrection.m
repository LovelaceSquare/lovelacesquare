function [correctedData, peakMask] = CosmicPeakCorrection(data, derivativeOrder, channelsToRemove, threshold)
% COSMICPEAKCORRECTION  Removes cosmic spikes from spectral data by 
% detecting sharp peaks in the derivative and interpolating over them.
%
% Authors: Adrián Gómez-Sánchez and Rodrigo Rocha de Oliveira
% Date Created: 2024-12-14
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0
%
% This function cleans a spectral data matrix by detecting and correcting
% abrupt spikes typically caused by cosmic rays. It operates by computing
% a row-wise derivative of the data, flagging points where the change
% exceeds a specified threshold, and then removing a defined number of
% neighboring channels around each detected spike. The gaps are filled by
% linear interpolation across each spectrum.
%
% INPUTS:
%     data            - Spectral data matrix [nSamples x nWavelengths].
%     derivativeOrder - Positive integer order of derivative (1 or 2 typical).
%     channelsToRemove- Integer: number of channels to remove on each side of a detected peak.
%     threshold       - Positive threshold. Points with derivative beyond this are flagged as cosmic peaks.
%
% OUTPUTS:
%   correctedData - The cleaned data matrix, same size as 'data'.
%   peakMask      - A logical matrix indicating which data points were replaced by interpolation.
%
% EXAMPLE:
%   % Suppose your spectral data is in a matrix 'X' (rows = samples, cols = wavelengths)
%   order = 1; 
%   chanRem = 2; 
%   thresh = 5;
%   [cleaned, mask] = CosmicPeakCorrection(X, order, chanRem, thresh);
%
% DISCLAIMER:
%   Authors and Lovelace's Square are not responsible for any 
%   issues, inaccuracies, or data loss arising from the use of this function.

    %% Validate Input Arguments
    if nargin < 4
        error('Usage: CosmicPeakCorrection(data, derivativeOrder, channelsToRemove, threshold)');
    end
    if ~isnumeric(data) || ndims(data) ~= 2 || isempty(data)
        error('Input "data" must be a non-empty 2D numeric matrix.');
    end
    if derivativeOrder < 1 || floor(derivativeOrder) ~= derivativeOrder
        error('Derivative order must be a positive integer (1, 2, 3, ...).');
    end
    [nSamples, nWavelengths] = size(data);
    if derivativeOrder >= nWavelengths
        error('Derivative order cannot exceed or equal the number of columns in "data".');
    end
    if channelsToRemove < 1 || floor(channelsToRemove) ~= channelsToRemove
        error('channelsToRemove must be a positive integer.');
    end
    if threshold <= 0
        error('threshold must be a positive number.');
    end

    % If data is complex, convert to real (as in the GUI)
    if ~isreal(data)
        warning('Complex data detected; using only the real part.');
        data = real(data);
    end

    %% 1) Compute the specified derivative
    % diff(..., derivativeOrder, 2) computes row-wise derivative along columns
    derivativeData = diff(data, derivativeOrder, 2); 
    % derivativeData is size [nSamples x (nWavelengths-derivativeOrder)]

    %% 2) Identify spikes above the absolute threshold
    absDerivativeData = abs(derivativeData);
    peakMaskShort = absDerivativeData > threshold; 
    % size(peakMaskShort) = [nSamples x (nWavelengths-derivativeOrder)]

    % The actual spike channels in the original data offset by 'derivativeOrder'
    % (e.g., derivative(1st order) at col k corresponds to difference between col k and col k+1 in original)
    [rowIdx, colIdxShort] = find(peakMaskShort);
    peakChannels = colIdxShort + derivativeOrder; 
    % Ensure channels do not exceed original matrix size
    peakChannels(peakChannels > nWavelengths) = [];

    %% 3) Remove channels around each detected peak
    % Create a logical mask over the original data size
    peakMask = false(nSamples, nWavelengths);

    for iPeak = 1:length(peakChannels)
        r = rowIdx(iPeak);
        c = peakChannels(iPeak);
        startIdx = max(c - channelsToRemove, 1);
        endIdx   = min(c + channelsToRemove, nWavelengths);
        peakMask(r, startIdx:endIdx) = true;
    end

    %% 4) Fill removed channels by linear interpolation along each row
    correctedData = data;
    % Mark peak positions with NaN
    correctedData(peakMask) = NaN;
    % Use fillmissing to interpolate. 'linear' along dimension 2 (columns),
    % 'EndValues','nearest' for boundary handling
    correctedData = fillmissing(correctedData, 'linear', 2, 'EndValues','nearest');

end
