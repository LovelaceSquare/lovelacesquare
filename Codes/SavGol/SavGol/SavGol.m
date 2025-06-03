function [filteredData] = SavGol(data, windowSize, polyOrder, derivOrder, edgeMethod)
% SavGol.  Apply Savitzky-Golay filtering to a 2D data matrix.
%
% REFERENCES:
%   Savitzky, Abraham, and Marcel J.E. Golay. "Smoothing and differentiation
%   of data by simplified least squares procedures." Analytical Chemistry
%   36.8 (1964): 1627-1639.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 2024-12-14
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0
%
% This function smooths each row of the input data (or computes its derivatives)
% using a Savitzky-Golay filter with customizable window size, polynomial order,
% derivative order, and edge handling strategy.
%
%
% INPUTS:
%   data       - 2D numeric matrix, size [nSamples x nPoints].
%                Each row is a separate spectrum or data series.
%   windowSize - Odd positive integer, length of the Sav-Gol window.
%                Must be > polyOrder.
%   polyOrder  - Integer polynomial order (< windowSize).
%   derivOrder - Integer derivative order (0 = smoothing only).
%                Must satisfy derivOrder <= polyOrder.
%   edgeMethod - (Optional) String specifying edge handling. One of:
%                  'None', 'Reflection', 'Replication', 'Extrapolation'
%                Default = 'None'.
%
% OUTPUT:
%   filteredData - Matrix of the same size as 'data', after Sav-Gol filtering
%                  is applied row-wise.
%
% EXAMPLE:
%   % Suppose X is [40 x 500], 40 spectra each with 500 data points:
%   winSize  = 11;   % Must be odd
%   polyOrd  = 3;    % e.g., cubic
%   derivOrd = 0;    % smoothing only
%   edgeMeth = 'Reflection';
%   Xfilt = SavGol(X, winSize, polyOrd, derivOrd, edgeMeth);
%
%   % Xfilt is now the Sav-Gol smoothed version of X.
%
%
% DISCLAIMER:
%   Authors and Lovelace's Square are not responsible for any issues,
%   inaccuracies, or data loss arising from the use of this function.

    %% 1) Validate inputs
    if nargin < 4
        error('Usage: SavGol(data, windowSize, polyOrder, derivOrder, [edgeMethod])');
    end
    if ~ismatrix(data) || isempty(data) || ~isnumeric(data)
        error('Input "data" must be a non-empty 2D numeric matrix [samples x points].');
    end
    [nSamples, nPoints] = size(data);

    if ~(isscalar(windowSize) && mod(windowSize,2) == 1 && windowSize > 1)
        error('windowSize must be an odd positive integer > 1.');
    end
    if ~(isscalar(polyOrder) && polyOrder >= 0 && polyOrder < windowSize)
        error('polyOrder must be an integer >= 0 and < windowSize.');
    end
    if ~(isscalar(derivOrder) && derivOrder >= 0 && derivOrder <= polyOrder)
        error('derivOrder must be an integer >= 0 and <= polyOrder.');
    end
    if nargin < 5 || isempty(edgeMethod)
        edgeMethod = 'None';
    else
        validMethods = {'None','Reflection','Replication','Extrapolation'};
        if ~ismember(edgeMethod, validMethods)
            error('edgeMethod must be one of: %s', strjoin(validMethods, ', '));
        end
    end

    % If data is complex, convert to real (mimicking the GUI's behavior)
    if ~isreal(data)
        warning('Complex data detected. Using only the real part of input data.');
        data = real(data);
    end

    %% 2) Precompute Sav-Gol filter coefficients
    coeffs = computeSavGolCoeffs(windowSize, polyOrder, derivOrder);

    %% 3) Filter each row according to edge handling
    filteredData = zeros(size(data));
    for rowIdx = 1:nSamples
        y = data(rowIdx, :);
        filteredData(rowIdx, :) = applySavGolRow(y, coeffs, windowSize, edgeMethod, polyOrder);
    end
end

%--------------------------------------------------------------------------
function coeffs = computeSavGolCoeffs(windowSize, polyOrder, derivOrder)
% COMPUTESAVGOLCOEFFS  Compute the convolution coefficients for a
% Savitzky-Golay filter of given windowSize, polyOrder, and derivOrder.

    halfWin = (windowSize - 1)/2;  % e.g., 5 if windowSize=11
    xVals   = (-halfWin : halfWin)';
    A       = zeros(windowSize, polyOrder+1);

    % Build Vandermonde matrix for polynomial fitting in local window
    % Columns: [x^0, x^1, x^2, ..., x^polyOrder]
    for n = 0:polyOrder
        A(:, n+1) = xVals.^n;
    end

    % Pseudoinverse of A
    A_pinv = pinv(A);

    % We next account for derivative scaling factors
    % derivative order = derivOrder
    derivScale = zeros(polyOrder+1, 1);
    for k = derivOrder:polyOrder
        % factorial(k)/factorial(k-derivOrder)
        % derivative of x^k => k*(k-1)*... = factorial(k)/factorial(k-derivOrder)
        derivScale(k+1) = factorial(k) / factorial(k - derivOrder);
    end

    % Combine: each column in A_pinv is a polynomial fit,
    % so multiply row-by-row by derivScale, then sum across
    % We'll do dot(derivScale, A_pinv) row-wise
    % Actually we want the final conv. coefficients as row vector:
    coeffs = (derivScale' * A_pinv)';
end

%--------------------------------------------------------------------------
function yFilt = applySavGolRow(y, coeffs, windowSize, edgeMethod, polyOrder)
% APPLYSAVGOLROW  Filter a single row vector y using the Savitzky-Golay
% convolution 'coeffs' and the specified edge handling method.

    halfWin = (windowSize - 1)/2;
    nPoints = length(y);

    switch edgeMethod
        case 'None'
            % Zero-padding approach
            yPadded = [zeros(1, halfWin), y, zeros(1, halfWin)];
            yConv   = conv(yPadded, coeffs, 'same');
            yFilt   = yConv(halfWin+1 : end-halfWin);

        case 'Reflection'
            % Reflect edges
            leftExt  = fliplr(y(1 : halfWin));
            rightExt = fliplr(y(end-halfWin+1 : end));
            yPadded  = [leftExt, y, rightExt];
            yConv    = conv(yPadded, coeffs, 'same');
            yFilt    = yConv(halfWin+1 : end-halfWin);

        case 'Replication'
            % Replicate edges
            leftExt  = repmat(y(1), 1, halfWin);
            rightExt = repmat(y(end), 1, halfWin);
            yPadded  = [leftExt, y, rightExt];
            yConv    = conv(yPadded, coeffs, 'same');
            yFilt    = yConv(halfWin+1 : end-halfWin);

        case 'Extrapolation'
            % Extrapolate edges by polynomial fit
            xData = 1:nPoints;

            % Left side fit on first 'windowSize' points
            xLeft  = xData(1 : windowSize);
            yLeft  = y(1  : windowSize);
            pLeft  = polyfit(xLeft, yLeft, polyOrder);
            xExtraLeft = (0 : (halfWin-1)) + 1 - halfWin;
            % e.g., if halfWin=5, we want xExtraLeft = -3, -2, -1, 0 => shift for indexing
            yExtraLeft = polyval(pLeft, xExtraLeft);

            % Right side fit on last 'windowSize' points
            xRight  = xData(end-windowSize+1 : end);
            yRight  = y(end-windowSize+1 : end);
            pRight  = polyfit(xRight, yRight, polyOrder);
            xExtraRight = (nPoints+1 : nPoints+halfWin);
            yExtraRight = polyval(pRight, xExtraRight);

            yPadded = [yExtraLeft, y, yExtraRight];
            yConv   = conv(yPadded, coeffs, 'same');
            yFilt   = yConv(halfWin+1 : end-halfWin);

        otherwise
            error('Unknown edgeMethod: %s', edgeMethod);
    end
end
