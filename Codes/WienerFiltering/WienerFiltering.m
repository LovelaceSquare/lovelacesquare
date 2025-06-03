function [filteredData] = WienerFiltering(...
    data, noiseVar, segmentLength, overlap, nfft)
% WienerFiltering. Apply Wiener filtering to spectral data.
%
% REFERENCES:
%   Wiener, Norbert. "Extrapolation, Interpolation, and Smoothing of
%   Stationary Time Series." MIT Press, 1949.
%
% Author:  Adrián Gómez-Sánchez.
% Contact: gomez.sanchez.adr@gmail.com
% Date:    2024-12-14
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: 1.0
%
% This function applies a Wiener filter to each row of spectral data by
% estimating the power spectral density (PSD) using Welch’s method and
% then filtering the signal in the frequency domain, effectively reducing
% noise based on a given noise variance estimate.
%
% INPUTS:
%   data          - 2D numeric matrix of size [nSamples x nWavelengths].
%                   Rows represent different samples/spectra, columns are the
%                   wavelength or spectral points.
%   noiseVar      - Estimated noise variance (scalar > 0). Larger values
%                   lead to more aggressive filtering.
%   segmentLength - Positive integer indicating the segment length used in
%                   Welch's PSD estimation (pwelch). Typically 128, 256, etc.
%   overlap       - Non-negative integer specifying the overlap between
%                   successive segments in pwelch. Must be < segmentLength.
%   nfft          - Positive integer for the FFT length in pwelch and
%                   Wiener filtering. Must be >= number of columns in data.
%
% OUTPUT:
%   filteredData  - The Wiener-filtered signals, size [nSamples x nWavelengths].
%
% EXAMPLE:
%   % Suppose 'X' is a [30 x 500] matrix with 30 spectra, each 500 points long.
%   noiseVar      = 0.01;
%   segmentLength = 256;
%   overlap       = 128;
%   nfftVal       = 512;
%   Xfilt = WienerFiltering(X, noiseVar, segmentLength, overlap, nfftVal);
%
% DEPENDENCIES:
%   pwelch  (MATLAB).
%
% DISCLAIMER:
%   Authors and Lovelace's Square are not responsible for any issues,
%   inaccuracies, or data loss arising from the use of this function.


    %% 1) Validate Inputs
    if nargin < 5
        error(['Usage: WienerFiltering(data, noiseVar, segmentLength, overlap, nfft)\n',...
               'All 5 inputs are required.']);
    end
    if ~ismatrix(data) || isempty(data) || ~isnumeric(data)
        error('Input "data" must be a non-empty 2D numeric matrix [samples x wavelengths].');
    end
    [nSamples, nWavelengths] = size(data);

    if ~isscalar(noiseVar) || noiseVar < 0
        error('noiseVar must be a non-negative scalar.');
    end
    if ~isscalar(segmentLength) || segmentLength < 1
        error('segmentLength must be a positive integer.');
    end
    if ~isscalar(overlap) || overlap < 0
        error('overlap must be a non-negative integer.');
    end
    if overlap >= segmentLength
        error('overlap must be strictly less than segmentLength.');
    end
    if ~isscalar(nfft) || nfft < nWavelengths
        error('nfft must be an integer >= number of columns in data.');
    end

    % If data is complex, convert to real (mimics the GUI's behavior).
    if ~isreal(data)
        warning('Complex data detected; converting to real part only.');
        data = real(data);
    end

    %% 2) Apply Wiener filter row-by-row
    filteredData = zeros(size(data));
    for iRow = 1:nSamples
        y = data(iRow, :).';       % column vector
        xhat = localWienerFilter(y, noiseVar, segmentLength, overlap, nfft);
        filteredData(iRow, :) = xhat';
    end
end

%--------------------------------------------------------------------------
function x_hat = localWienerFilter(y, noiseVar, segmentLength, overlap, nfft)
% LOCALWIENERFILTER  Applies the Wiener filter logic to a single spectral row.
%
%   y          : Column vector with the spectral data.
%   noiseVar   : Estimated noise variance.
%   segmentLength, overlap, nfft
%              : Parameters for PSD estimation and Wiener filtering.
%
%   x_hat      : The Wiener-filtered signal, same length as y.

    % Ensure NFFT is at least nextpow2(length(y)) if the user chooses so
    if nfft < length(y)
        nfft = 2^nextpow2(length(y));
    end

    % 1) Estimate power spectral density (PSD) of the noisy signal using Welch's method
    [Pyy, ~] = pwelch(y, segmentLength, overlap, nfft, 'twosided');

    % 2) Noise PSD (assume flat = noiseVar)
    Snn = noiseVar * ones(size(Pyy));

    % 3) Estimate signal PSD: Sxx = Pyy - Snn, but clip at eps to avoid negatives
    Sxx = max(Pyy - Snn, eps);

    % 4) Compute Wiener filter frequency response H(f) = Sxx / (Sxx + Snn)
    H = Sxx ./ (Sxx + Snn);

    % 5) Apply filter in frequency domain: X_hat = Y * H
    Y = fft(y, nfft);
    X_hat = Y .* H;

    % 6) Inverse FFT
    x_hat_full = ifft(X_hat, 'symmetric');
    x_hat = x_hat_full(1:length(y));  % keep original length
end
