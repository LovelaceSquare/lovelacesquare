function [scaledMatrix, scalingParams] = autoscale(inputMatrix, direction)
% AUTOSCALE Perform autoscaling on a 2D data matrix.
%
% Author:         Adrián Gómez-Sánchez
% Date Created:   2021-03-18
% License:        MIT
% Reviewed by:    Lovelace's Square
% Version:        v2.0
%
% This function centers and scales the input matrix to have zero mean and 
% unit variance. Autoscaling can be performed either across columns or rows.
%
% INPUTS:
%   inputMatrix - 2D numeric matrix [nRows x nCols]. Each row can represent 
%                a sample, and each column a feature.
%   direction   - String specifying the direction of scaling:
%                 * 'column' => Scale each column independently.
%                 * 'row'    => Scale each row independently.
%                 If omitted, defaults to 'column'.
%
% OUTPUTS:
%   scaledMatrix  - Autoscaled matrix, same size as inputMatrix.
%   scalingParams - Structure containing the means and standard deviations 
%                   used for scaling:
%                   - scalingParams.mean: Mean values.
%                   - scalingParams.std:  Standard deviation values.
%
% EXAMPLES:
%   % Example 1: Scale columns of the input matrix
%   [X_scaled, params] = autoscale(X, 'column');
%
%   % Example 2: Scale rows of the input matrix
%   [X_scaled, params] = autoscale(X, 'row');
%
% DISCLAIMER:
%   Authors and Lovelace's Square are not responsible for any issues,
%   inaccuracies, or data loss arising from the use of this function.

    %% Input Validation
    if nargin < 1
        error('AUTOSCALE requires at least one input argument: inputMatrix.');
    end
    if ~isnumeric(inputMatrix) || ~ismatrix(inputMatrix)
        error('Input must be a numeric 2D matrix.');
    end
    
    if nargin < 2 || isempty(direction)
        direction = 'column'; % Default direction
    end
    
    % Validate direction
    validDirections = {'column', 'row'};
    if ~ischar(direction) && ~isstring(direction)
        error('Direction must be a string: ''column'' or ''row''.');
    end
    direction = lower(char(direction));
    if ~ismember(direction, validDirections)
        error('Invalid direction. Choose either ''column'' or ''row''.');
    end

    %% Initialize Output Variables
    scalingParams = struct('mean', [], 'std', []);

    %% Compute Scaling Parameters
    switch direction
        case 'column'
            scalingParams.mean = mean(inputMatrix, 1); % 1 x nCols
            scalingParams.std  = std(inputMatrix, 0, 1); % 1 x nCols
        case 'row'
            scalingParams.mean = mean(inputMatrix, 2); % nRows x 1
            scalingParams.std  = std(inputMatrix, 0, 2); % nRows x 1
    end

    %% Handle Zero Standard Deviations
    zeroStdMask = (scalingParams.std == 0);
    if any(zeroStdMask)
        warning('Some %s have zero standard deviation. These will not be scaled.', direction);
        scalingParams.std(zeroStdMask) = 1; % Prevent division by zero
    end

    %% Perform Autoscaling
    if strcmp(direction, 'column')
        % Utilize implicit expansion (MATLAB R2016b and later)
        scaledMatrix = (inputMatrix - scalingParams.mean) ./ scalingParams.std;
        
        % For MATLAB versions prior to R2016b, uncomment the following lines:
        % scaledMatrix = (inputMatrix - repmat(scalingParams.mean, size(inputMatrix, 1), 1)) ./ ...
        %                repmat(scalingParams.std, size(inputMatrix, 1), 1);
                
        % Optionally, set scaled values to zero where std was originally zero
        scaledMatrix(:, zeroStdMask) = 0;
        
    else % direction == 'row'
        % Utilize implicit expansion (MATLAB R2016b and later)
        scaledMatrix = (inputMatrix - scalingParams.mean) ./ scalingParams.std;
        
        % For MATLAB versions prior to R2016b, uncomment the following lines:
        % scaledMatrix = (inputMatrix - repmat(scalingParams.mean, 1, size(inputMatrix, 2))) ./ ...
        %                repmat(scalingParams.std, 1, size(inputMatrix, 2)));
                
        % Optionally, set scaled values to zero where std was originally zero
        scaledMatrix(zeroStdMask, :) = 0;
    end

end
