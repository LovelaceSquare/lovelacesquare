function [correctedData, referenceSpec] = MSC(data, refType, refIndex)
% MSC.  Perform classic Multiplicative Scatter Correction on a 2-D data matrix.
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
% Version: v1.3   % (2025-07-14) – add explicit-spectrum option, keep full docs
%
% This function removes additive and multiplicative scatter effects
% from row-wise spectral data by aligning each row to a chosen reference
% spectrum (mean, a specific row, or a user-supplied spectrum).
%
% INPUTS:
%   data      – 2-D numeric matrix [nSamples × nWavelengths]. Each row is
%               a spectrum or data series.
%   refType   – (optional) string/char specifying the reference mode:
%                 * 'Mean Spectrum'        (default)
%                 * 'Reference Index'
%                 * 'Reference Spectrum'   (NEW: pass a 1 × nWavelengths vector)
%   refIndex  – (optional)
%                 • when refType = 'Reference Index': positive integer row #
%                 • when refType = 'Reference Spectrum': numeric row vector
%                 • ignored for 'Mean Spectrum'
%
% OUTPUTS:
%   correctedData – MSC-corrected data (same size as *data*).
%   referenceSpec – The spectrum actually used as reference [1 × nWavelengths].
%
% EXAMPLE:
%   % X is [50 × 400]
%   [XcorrMean, refMean] = MSC(X);                       % mean reference
%   [Xcorr10,  ref10]   = MSC(X, 'Reference Index', 10); % sample #10
%   extSpec = rand(1,400);                               % external spectrum
%   [XcorrExt, refExt]  = MSC(X, 'Reference Spectrum', extSpec);
%
% DISCLAIMER:
%   Authors and Lovelace's Square are not responsible for any issues,
%   inaccuracies, or data loss arising from the use of this function.
% -------------------------------------------------------------------------

%% 0) Handle optional arguments & validate primary input
if nargin < 1
    error('Usage: MSC(data, [refType], [refIndex|refSpectrum])');
end
if nargin < 2 || isempty(refType),  refType = 'Mean Spectrum';  end
if nargin < 3,                     refIndex = [];              end

validateattributes(data, {'numeric'}, {'2d','nonempty'}, mfilename,'data',1);
[nSamples, nWavelengths] = size(data);

% Convert complex data to real part only (rare in NIR work, but safe)
if ~isreal(data)
    warning('Complex data detected; using only the real part.');
    data = real(data);
end

%% 1) Determine the reference spectrum
refType = lower(string(refType));

switch refType
    case 'mean spectrum'
        referenceSpec = mean(data, 1);

    case 'reference index'
        if isempty(refIndex), refIndex = 1; end
        validateattributes(refIndex,{'numeric'}, ...
            {'scalar','integer','positive','<=',nSamples}, ...
            mfilename,'refIndex',3);
        referenceSpec = data(refIndex, :);

    case {'reference spectrum','custom spectrum'}
        validateattributes(refIndex,{'numeric'}, ...
            {'row','numel',nWavelengths}, ...
            mfilename,'refSpectrum',3);
        referenceSpec = refIndex;

    otherwise
        error('Unknown refType "%s". Use ''Mean Spectrum'', ''Reference Index'', or ''Reference Spectrum''.', refType);
end

%% 2) Classic MSC correction
correctedData = msc_correction(data, referenceSpec);
end %-- MSC main function
% =====================================================================
function Xcorr = msc_correction(X, refSpectrum)
% Vectorised implementation of classic MSC.
%   X           : [nSamples × nWavelengths] input data.
%   refSpectrum : [1 × nWavelengths] reference spectrum.
[nSamples, nVars] = size(X);

% Design matrix: [offset  reference]
M = [ones(nVars,1), refSpectrum(:)];
A = (M' * M) \ M';        % 2 × nVars pseudo-inverse

ab = A * X.';             % 2 × nSamples  (offsets; scale factors)
offsets   = ab(1,:);
scaleFacs = ab(2,:);

% Guard against zero scale factors
zeroMask = (scaleFacs == 0);
if any(zeroMask)
    warning('%d scale factor(s) evaluated to zero – affected rows left unscaled.', nnz(zeroMask));
    scaleFacs(zeroMask) = 1;
end

Xcorr = (X - offsets.')./scaleFacs.';
end %-- msc_correction
