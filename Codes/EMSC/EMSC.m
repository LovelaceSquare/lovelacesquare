function [correctedData, referenceSpec] = EMSC(data, refType, polyOrder)
% EMSC.  Perform Extended Multiplicative Scatter Correction on a 2D data matrix.
%
% REFERENCES:
%   Martens, Harald, and E. Stark. "Extended multiplicative signal correction
%   and spectral interference subtraction: new preprocessing methods for near
%   infrared spectroscopy." Journal of Pharmaceutical and Biomedical Analysis
%   9.8 (1991): 625-635.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 2024-12-14
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v1.0
%
% This function removes multiplicative and additive scatter effects from
% row-wise spectral data (or analogous measurements) by aligning each row
% to a chosen reference spectrum (mean or median) and optionally modeling
% baseline shifts with polynomial terms.
%
%
% INPUTS:
%   data      - 2D numeric matrix, size [nSamples x nWavelengths].
%               Each row is a separate spectrum or data series.
%   refType   - String specifying the reference spectrum: 'Mean' or 'Median'.
%   polyOrder - Non-negative integer specifying the polynomial order for
%               baseline modeling (0 for no polynomial, 1 for linear, etc.).
%
% OUTPUTS:
%   correctedData - Matrix of the same size as 'data' after EMSC correction
%                   (row-wise).
%   referenceSpec - The actual reference spectrum used (1 x nWavelengths).
%
% EXAMPLE:
%   % Suppose X is [40 x 500], i.e., 40 spectra each with 500 points:
%   refT = 'Mean';      % or 'Median'
%   pOrd = 1;           % linear baseline
%   [Xcorr, refSpec] = emsc(X, refT, pOrd);
%
%   % Xcorr is the EMSC-corrected data. refSpec is the computed reference spectrum.
%
%
% DISCLAIMER:
%   Authors and Lovelace's Square are not responsible for any issues,
%   inaccuracies, or data loss arising from the use of this function.

    %% 1) Validate Inputs
    if nargin < 3
        error('Usage: emsc(data, refType, polyOrder). All inputs are required.');
    end
    if ~isnumeric(data) || ndims(data) ~= 2 || isempty(data)
        error('Input "data" must be a non-empty 2D numeric matrix [samples x wavelengths].');
    end
    [nSamples, nWavelengths] = size(data);

    if ~ischar(refType) && ~isstring(refType)
        error('refType must be a character array or string: "Mean" or "Median".');
    end
    if ~any(strcmpi(refType, {'Mean','Median'}))
        error('refType must be "Mean" or "Median".');
    end
    if ~isnumeric(polyOrder) || polyOrder < 0 || floor(polyOrder) ~= polyOrder
        error('polyOrder must be a non-negative integer (0, 1, 2, ...).');
    end

    % If data is complex, warn and convert to real
    if ~isreal(data)
        warning('Complex data detected; using only the real part.');
        data = real(data);
    end

    %% 2) Determine the reference spectrum
    switch lower(refType)
        case 'mean'
            referenceSpec = mean(data, 1);
        case 'median'
            referenceSpec = median(data, 1);
    end

    %% 3) Perform EMSC correction
    correctedData = emsc_correction(data, referenceSpec, polyOrder);

end

%--------------------------------------------------------------------------
function Xcorr = emsc_correction(X, refSpectrum, polyOrder)
% EMSC_CORRECTION.  Internal helper to apply EMSC row-wise.
%
%   X          : [nSamples x nWavelengths] input data.
%   refSpectrum: [1 x nWavelengths], the chosen reference (Mean or Median).
%   polyOrder  : Non-negative integer for polynomial baseline order.
%
% Returns:
%   Xcorr : EMSC-corrected data, same size as X.

    [nSamples, nVars] = size(X);

    % Build design matrix M
    M = ones(nVars, 1);          % Column for constant offset
    M = [M, refSpectrum(:)];     % Column for reference spectrum
    xAxis = linspace(0, 1, nVars)';  % Normalized axis for polynomial baseline
    if polyOrder > 0
        for p = 1:polyOrder
            M = [M, xAxis.^p];
        end
    end

    % Precompute pseudo-inverse of M
    A = (M' * M) \ (M');

    % Correct each spectrum
    Xcorr = zeros(size(X));
    for i = 1:nSamples
        a = A * (X(i,:))';

        % a(1) = offset, a(2) = scale factor
        offset   = a(1);
        scaleFac = a(2);

        % Subtract polynomial baseline (if any)
        polyBaseline = 0;
        if polyOrder > 0
            polyCoeffs = a(3 : 2+polyOrder);
            for p = 1:polyOrder
                polyBaseline = polyBaseline + polyCoeffs(p)*(xAxis.^p)';
            end
        end

        % EMSC-corrected row
        Xcorr(i,:) = (X(i,:) - offset - polyBaseline) / scaleFac;
    end
end
