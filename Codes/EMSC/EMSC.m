function [correctedData, referenceSpec] = EMSC(data, refType, polyOrder, refSpectrum)
% EMSC.  Perform Extended Multiplicative Scatter Correction on a 2-D data matrix.
%
% REFERENCES:
%   Martens, Harald, and E. Stark. "Extended multiplicative signal correction 
%   and spectral interference subtraction: new preprocessing methods for near 
%   infrared spectroscopy." Journal of Pharmaceutical and Biomedical Analysis 
%   9 (8) (1991): 625-635.
%
% Author:  Adrián Gómez-Sánchez
% Date Created: 2024-12-14
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v1.2   % (2025-07-14) – add External reference + fix Xcorr assignment bug
%
% This function removes multiplicative and additive scatter effects from 
% row-wise spectral data by aligning each row to a chosen reference spectrum 
% (dataset mean, dataset median, or user-supplied) and, optionally, modelling 
% baseline drift with a polynomial of arbitrary order.
%
% INPUTS:
%   data         – 2-D numeric matrix [nSamples × nWavelengths]. Each row is
%                  a separate spectrum or data series.
%   refType      – String: 'Mean', 'Median', or 'External'.
%   polyOrder    – Non-negative integer specifying polynomial order for
%                  baseline modelling (0 = none, 1 = linear, 2 = quadratic, …).
%   refSpectrum  – Row vector [1 × nWavelengths] used only when
%                  refType = 'External'.
%
% OUTPUTS:
%   correctedData – Matrix of same size as *data* after EMSC correction.
%   referenceSpec – The actual reference spectrum used [1 × nWavelengths].
%
% EXAMPLE:
%   % X is [40 × 500] (40 spectra, 500 wavelengths):
%   [Xcorr1, ref1] = EMSC(X, 'Mean',    1);            % linear baseline
%   [Xcorr2, ref2] = EMSC(X, 'Median',  2);            % quadratic baseline
%   extRef = load('certified_ref.mat','spectrum').spectrum; % 1×500 reference
%   [Xcorr3, ref3] = EMSC(X, 'External', 2, extRef);
%
% DISCLAIMER:
%   Authors and Lovelace's Square are not responsible for any issues, 
%   inaccuracies, or data loss arising from the use of this function.
% -------------------------------------------------------------------------

%% 1) Validate Inputs
if nargin < 3
    error('Usage: EMSC(data, refType, polyOrder [, refSpectrum])');
end
validateattributes(data, {'numeric'}, {'2d','nonempty'}, mfilename,'data',1);
[nSamples, nWavelengths] = size(data);

if ~(ischar(refType) || isstring(refType))
    error('refType must be ’’Mean’’, ’’Median’’, or ’’External’’.');
end
refType = lower(string(refType));
if ~any(strcmp(refType, ["mean","median","external"]))
    error('refType must be ''Mean'', ''Median'', or ''External''.');
end

validateattributes(polyOrder, {'numeric'}, ...
    {'scalar','integer','>=',0}, mfilename,'polyOrder',3);

if refType == "external"
    if nargin < 4
        error('refSpectrum must be provided when refType = ''External''.');
    end
    validateattributes(refSpectrum, {'numeric'}, ...
        {'row','numel',nWavelengths}, mfilename,'refSpectrum',4);
end

% Convert complex data to real part only
if ~isreal(data)
    warning('Complex data detected; using only the real part.');
    data = real(data);
end

%% 2) Determine the reference spectrum
switch refType
    case "mean"
        referenceSpec = mean(data, 1);
    case "median"
        referenceSpec = median(data, 1);
    case "external"
        referenceSpec = refSpectrum;
end

%% 3) Perform EMSC correction
correctedData = emsc_correction(data, referenceSpec, polyOrder);

end %--- EMSC (main) ---

% =========================================================================
function Xcorr = emsc_correction(X, refSpectrum, polyOrder)
% EMSC_CORRECTION  Internal helper – apply EMSC row-wise to matrix X.

[nSamples, nVars] = size(X);

% ---------- build design matrix -----------------------------------------
M = [ones(nVars,1), refSpectrum(:)];      % columns: offset, scale
xAxis = linspace(0,1,nVars).';            % normalised wavelength axis
for p = 1:polyOrder
    M = [M, xAxis.^p];                    %#ok<AGROW>
end

% Moore-Penrose pseudo-inverse (small matrix, inexpensive)
A = (M' * M) \ M';                        % (nParams × nVars)

% ---------- solve LS for all spectra at once ----------------------------
coeff = A * X.';                          % (nParams × nSamples)
offsets   = coeff(1 ,:);
scaleFacs = coeff(2 ,:);
polyCoeffs = coeff(3:end ,:);             % empty if polyOrder == 0

% ---------- compute polynomial baseline for each sample -----------------
if polyOrder > 0
    polyTerms = zeros(polyOrder, nVars);
    for p = 1:polyOrder
        polyTerms(p,:) = (xAxis.^p).';
    end
    baseline = polyCoeffs.' * polyTerms;  % [nSamples × nVars]
else
    baseline = zeros(size(X));
end

% ---------- guard against division by zero ------------------------------
zeroScale = (scaleFacs == 0);
if any(zeroScale)
    warning('%d scale factor(s) were zero – affected rows set to scale = 1.', nnz(zeroScale));
    scaleFacs(zeroScale) = 1;
end

% ---------- final correction --------------------------------------------
Xcorr = (X - offsets.' - baseline) ./ scaleFacs.';

end %--- emsc_correction
