function normalizedMatrix = normMatrix(data, normType, dimension)
% normMatrix. Normalizes the input matrix using specified norms.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 2024-12-14
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: 1.0
%
% The normMatrix function normalizes the input matrix using various norms,
% including maximum value ('max'), Euclidean ('euclidean'), L1 ('l1'), L2 ('l2'),
% L-infinity ('linf'), and Frobenius ('frobenius') norms. Note that the "euclidean"
% option applies row- or column-wise normalization using the Euclidean norm (i.e.,
% sqrt(sum(x_i^2))), while "l2" is defined only for normalization over the entire
% matrix. This distinction is provided so you can choose the normalization method
% that best suits your data structure.
%
% INPUTS:
%   data (array)      : The input matrix to be normalized.
%   normType (string) : The type of normalization:
%                        'max'        - Maximum value normalization.
%                        'euclidean'  - Euclidean norm normalization (row/column).
%                        'l1'         - L1 norm normalization.
%                        'l2'         - L2 norm normalization (entire matrix only).
%                        'linf'       - L-infinity norm normalization.
%                        'frobenius'  - Frobenius norm normalization (entire matrix only).
%   dimension (string): The dimension along which to perform normalization:
%                        'all'    - Normalization over the entire matrix.
%                        'row'    - Row-wise normalization.
%                        'column' - Column-wise normalization.
%
% OUTPUT:
%   normalizedMatrix (array): The normalized matrix based on the specified norm.
%
% USE EXAMPLE:
%   normalizedMatrix = normMatrix(data, 'euclidean', 'row');
%
% DISCLAIMER:
% Authors and Lovelace's Square are not responsible for any issues, inaccuracies,
% or data loss arising from the use of this function.

switch lower(normType)
    case 'max'
        switch lower(dimension)
            case 'all'
                denom = max(data(:));
                if denom == 0
                    warndlg('All elements are zero. Normalization will produce NaNs.','Normalization Warning');
                end
                normalizedMatrix = data ./ denom;
            case 'row'
                denom = max(data, [], 2);
                if any(denom == 0)
                    warndlg('At least one row has a zero maximum. Normalization will produce NaNs.','Normalization Warning');
                end
                normalizedMatrix = bsxfun(@rdivide, data, denom);
            case 'column'
                denom = max(data, [], 1);
                if any(denom == 0)
                    warndlg('At least one column has a zero maximum. Normalization will produce NaNs.','Normalization Warning');
                end
                normalizedMatrix = bsxfun(@rdivide, data, denom);
            otherwise
                error('Invalid dimension. Please choose from: ''all'', ''row'', ''column''.');
        end

    case 'euclidean'
        switch lower(dimension)
            case 'all'
                denom = norm(data, 'fro');
                if denom == 0
                    warndlg('The entire matrix norm is zero. Normalization will produce NaNs.','Normalization Warning');
                end
                normalizedMatrix = data ./ denom;
            case 'row'
                denom = sqrt(sum(data.^2, 2));
                if any(denom == 0)
                    warndlg('At least one row has zero Euclidean norm. Normalization will produce NaNs.','Normalization Warning');
                end
                normalizedMatrix = bsxfun(@rdivide, data, denom);
            case 'column'
                denom = sqrt(sum(data.^2, 1));
                if any(denom == 0)
                    warndlg('At least one column has zero Euclidean norm. Normalization will produce NaNs.','Normalization Warning');
                end
                normalizedMatrix = bsxfun(@rdivide, data, denom);
            otherwise
                error('Invalid dimension. Please choose from: ''all'', ''row'', ''column''.');
        end

    case 'l1'
        switch lower(dimension)
            case 'all'
                denom = norm(data, 1);
                if denom == 0
                    warndlg('The entire matrix L1 norm is zero. Normalization will produce NaNs.','Normalization Warning');
                end
                normalizedMatrix = data ./ denom;
            case 'row'
                denom = sum(abs(data), 2);
                if any(denom == 0)
                    warndlg('At least one row has zero L1 norm. Normalization will produce NaNs.','Normalization Warning');
                end
                normalizedMatrix = bsxfun(@rdivide, data, denom);
            case 'column'
                denom = sum(abs(data), 1);
                if any(denom == 0)
                    warndlg('At least one column has zero L1 norm. Normalization will produce NaNs.','Normalization Warning');
                end
                normalizedMatrix = bsxfun(@rdivide, data, denom);
            otherwise
                error('Invalid dimension. Please choose from: ''all'', ''row'', ''column''.');
        end

    case 'l2'
        switch lower(dimension)
            case 'all'
                % Note: 'l2' is defined for the entire matrix only. This is why there is both
                % 'euclidean' (for row/column normalization) and 'l2' (for entire matrix normalization).
                denom = norm(data, 2);
                if denom == 0
                    warndlg('The entire matrix L2 norm is zero. Normalization will produce NaNs.','Normalization Warning');
                end
                normalizedMatrix = data ./ denom;
            otherwise
                error('L2 norm is not defined for row-wise or column-wise normalization.');
        end

    case 'linf'
        switch lower(dimension)
            case 'all'
                denom = norm(data, inf);
                if denom == 0
                    warndlg('The entire matrix L-infinity norm is zero. Normalization will produce NaNs.','Normalization Warning');
                end
                normalizedMatrix = data ./ denom;
            case 'row'
                denom = max(abs(data), [], 2);
                if any(denom == 0)
                    warndlg('At least one row has zero L-infinity norm. Normalization will produce NaNs.','Normalization Warning');
                end
                normalizedMatrix = bsxfun(@rdivide, data, denom);
            case 'column'
                denom = max(abs(data), [], 1);
                if any(denom == 0)
                    warndlg('At least one column has zero L-infinity norm. Normalization will produce NaNs.','Normalization Warning');
                end
                normalizedMatrix = bsxfun(@rdivide, data, denom);
            otherwise
                error('Invalid dimension. Please choose from: ''all'', ''row'', ''column''.');
        end

    case 'frobenius'
        switch lower(dimension)
            case 'all'
                denom = norm(data, 'fro');
                if denom == 0
                    warndlg('The entire matrix Frobenius norm is zero. Normalization will produce NaNs.','Normalization Warning');
                end
                normalizedMatrix = data ./ denom;
            otherwise
                error('Frobenius norm is not defined for row-wise or column-wise normalization.');
        end

    otherwise
        error('Invalid norm type. Please choose from: ''max'', ''euclidean'', ''l1'', ''l2'', ''linf'', ''frobenius''.');
end

end
