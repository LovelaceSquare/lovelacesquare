classdef DataValidator
% DataValidator. Validation utilities for CosmicPeakCorrection GUI.
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-16
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 2.0
%
% Validates input data matrices and correction parameters before
% processing. All methods are static -- no instantiation required.
%
% METHODS (Static):
%   validateParams - Validate derivative order, channels to remove, threshold
%   validateData   - Validate input data matrix (numeric, 2D, finite)

    methods (Static)
        function [valid, message] = validateParams(derivativeOrder, channelsToRemove, threshold, nCols)
            % VALIDATEPARAMS Validate correction parameters.
            %
            %   [valid, message] = DataValidator.validateParams(derivativeOrder, channelsToRemove, threshold, nCols)
            %
            %   Inputs:
            %       derivativeOrder  - Order of derivative
            %       channelsToRemove - Number of channels to expand on each side
            %       threshold        - Detection threshold
            %       nCols            - Number of columns in data matrix
            %
            %   Outputs:
            %       valid   - true if all parameters are valid
            %       message - Error message if invalid, empty string if valid

            valid = false;
            message = '';

            % Validate derivativeOrder
            if ~isnumeric(derivativeOrder) || ~isscalar(derivativeOrder)
                message = 'Derivative order must be a numeric scalar.';
                return;
            end
            if derivativeOrder < 1
                message = 'Derivative order must be >= 1.';
                return;
            end
            if floor(derivativeOrder) ~= derivativeOrder
                message = 'Derivative order must be an integer.';
                return;
            end
            if derivativeOrder >= nCols
                message = sprintf('Derivative order must be less than the number of channels (%d).', nCols);
                return;
            end

            % Validate channelsToRemove
            if ~isnumeric(channelsToRemove) || ~isscalar(channelsToRemove)
                message = 'Channels to remove must be a numeric scalar.';
                return;
            end
            if channelsToRemove < 0
                message = 'Channels to remove must be >= 0.';
                return;
            end
            if floor(channelsToRemove) ~= channelsToRemove
                message = 'Channels to remove must be an integer.';
                return;
            end

            % Validate threshold
            if ~isnumeric(threshold) || ~isscalar(threshold)
                message = 'Threshold must be a numeric scalar.';
                return;
            end
            if threshold <= 0
                message = 'Threshold must be > 0.';
                return;
            end
            if ~isfinite(threshold)
                message = 'Threshold must be a finite value.';
                return;
            end

            valid = true;
        end

        function [valid, message] = validateData(data)
            % VALIDATEDATA Validate input data matrix.
            %
            %   [valid, message] = DataValidator.validateData(data)
            %
            %   Inputs:
            %       data - Data matrix to validate
            %
            %   Outputs:
            %       valid   - true if data is valid
            %       message - Error message if invalid, empty string if valid

            valid = false;
            message = '';

            if isempty(data)
                message = 'Data is empty.';
                return;
            end

            if ~isnumeric(data)
                message = 'Data must be numeric.';
                return;
            end

            if ~ismatrix(data)
                message = 'Data must be a 2D matrix.';
                return;
            end

            if size(data, 2) < 3
                message = 'Data must have at least 3 channels (columns).';
                return;
            end

            if any(~isfinite(data(:)))
                message = 'Data contains non-finite values (NaN or Inf).';
                return;
            end

            valid = true;
        end
    end
end
