function [correctedData, referenceSpec] = MSC(data, refType, refIndex)
% MSC.  Perform Multiplicative Scatter Correction on a 2D data matrix.
%
% REFERENCES:
%   Isaksson, T., & Næs, T. (1988). The effect of multiplicative scatter
%   correction (MSC) and linearity improvement in NIR spectroscopy.
%   Applied Spectroscopy, 42(7), 1273-1284.
%
% Author:  Adrián Gómez-Sánchez
% Date Created: 2024-12-14
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0
%
% This function removes additive and multiplicative scatter effects
% from row-wise spectral data by aligning each row to a chosen reference
% spectrum (mean or a specific row). Internally, it uses the same *structure*
% as the EMSC function, but it forces the polynomial order to 0. This makes
% it equivalent to classic MSC.
%
% INPUTS:
%   data      - 2D numeric matrix [nSamples x nWavelengths]. Each row is
%               a spectrum or data series.
%   refType   - String specifying the reference spectrum to use:
%               * 'Mean Spectrum'   => Use the mean of all rows.
%               * 'Reference Index' => Use data(refIndex,:) as reference.
%   refIndex  - Positive integer specifying which row to use when
%               refType = 'Reference Index'. If out-of-range, defaults
%               to 1. If refType = 'Mean Spectrum', this value is still
%               required but ignored.
%
% OUTPUTS:
%   correctedData - MSC-corrected data, same size as input 'data'.
%   referenceSpec - The actual reference spectrum used [1 x nWavelengths].
%
% EXAMPLE:
%   % Suppose X is [50 x 400].
%   % 1) Use mean spectrum as reference:
%   [XcorrMean, refMean] = MSC(X, 'Mean Spectrum', 1);
%
%   % 2) Use sample #10 as reference:
%   [Xcorr10, ref10] = MSC(X, 'Reference Index', 10);
%
% DISCLAIMER:
%   Authors and Lovelace's Square are not responsible for any issues,
%   inaccuracies, or data loss arising from the use of this function.

    %% 1) Validate Inputs
    if nargin < 3
        error('Usage: MSC(data, refType, refIndex)');
    end
    if ~isnumeric(data) || ndims(data) ~= 2 || isempty(data)
        error('Input "data" must be a non-empty 2D numeric matrix [samples x wavelengths].');
    end
    [nSamples, nWavelengths] = size(data);

    if ~ischar(refType) && ~isstring(refType)
        error('refType must be a character array or string: "Mean Spectrum" or "Reference Index".');
    end

    % Convert complex data to real if necessary
    if ~isreal(data)
        warning('Complex data detected; using only the real part.');
        data = real(data);
    end

    %% 2) Determine the reference spectrum
    switch lower(refType)
        case 'mean spectrum'
            referenceSpec = mean(data, 1);
        case 'reference index'
            refIndex = round(refIndex);
            if refIndex < 1 || refIndex > nSamples
                warning('Reference index %d out of range. Defaulting to 1.', refIndex);
                refIndex = 1;
            end
            referenceSpec = data(refIndex, :);
        otherwise
            error('Unknown refType "%s". Use "Mean Spectrum" or "Reference Index".', refType);
    end

    %% 3) Force standard MSC (i.e., EMSC with polyOrder = 0)
    polyOrder = 0;  % <--- This is the key step to make it "classic" MSC

    correctedData = msc_correction(data, referenceSpec, polyOrder);

end


%--------------------------------------------------------------------------
function Xcorr = msc_correction(X, refSpectrum, polyOrder)
% MSC_CORRECTION. Internal helper that mimics the EMSC approach
% but forces polyOrder=0 for classic MSC.
%
%   X           : [nSamples x nWavelengths] input data.
%   refSpectrum : [1 x nWavelengths], reference.
%   polyOrder   : forced to 0 here for standard MSC.
%
% Returns:
%   Xcorr : MSC-corrected data, same size as X.

    [nSamples, nVars] = size(X);

    % Build design matrix M:
    % Column 1: offset (ones)
    % Column 2: reference spectrum
    % (No polynomial terms beyond offset, because polyOrder=0)
    M = [ones(nVars,1), refSpectrum(:)];

    % Precompute pseudo-inverse of M
    A = (M' * M) \ (M');

    % Correct each row
    Xcorr = zeros(size(X));
    for i = 1:nSamples
        % Solve for offset (a) and scale factor (b)
        a = A * (X(i,:))';  % a(1) = offset, a(2) = scale factor

        offset   = a(1);
        scaleFac = a(2);

        % No extra polynomial baseline because polyOrder=0
        % so the correction is: (row - offset) / scaleFac
        Xcorr(i,:) = (X(i,:) - offset) / scaleFac;
    end
end
