classdef DataValidator
    % DataValidator - Static validation methods for LAsLS baseline correction data and parameters.
    %
    % Provides input validation for spectral data matrices, algorithm parameters,
    % and interval definitions before they are passed to the LASLSCorrector.
    %
    % Authors: Adrián Gómez-Sánchez, Berta Torres-Cobos, Rodrigo Rocha de Oliveira
    % Date Created: 2024-12-16
    % License: MIT
    % Repository: https://github.com/LovelaceSquare/lovelacesquare

    % Version: 1.1

    methods (Static)
        function [isValid, msg] = validateData(data)
            % validateData - Validate that input data is a suitable numeric matrix.
            %
            % Inputs:
            %   data - The data to validate
            %
            % Outputs:
            %   isValid - Logical, true if data is valid
            %   msg     - Character vector describing the issue (empty if valid)

            isValid = true;
            msg = '';

            if isempty(data)
                isValid = false;
                msg = 'Data is empty.';
                return;
            end

            if ~isnumeric(data)
                isValid = false;
                msg = 'Data must be numeric.';
                return;
            end

            if ~ismatrix(data)
                isValid = false;
                msg = 'Data must be a 2D matrix.';
                return;
            end

            if any(isnan(data(:)))
                isValid = false;
                msg = 'Data contains NaN values. Please clean the data first.';
                return;
            end

            if any(isinf(data(:)))
                isValid = false;
                msg = 'Data contains Inf values. Please clean the data first.';
                return;
            end

            if size(data, 2) < 4
                isValid = false;
                msg = 'Data must have at least 4 spectral channels (columns).';
                return;
            end
        end

        function [isValid, msg] = validateParams(lambda, p, mu, maxIter, tolerance)
            % validateParams - Validate LASLS algorithm parameters.
            %
            % Inputs:
            %   lambda    - Global smoothness parameter (must be >= 0)
            %   p         - Global asymmetry parameter (must be in (0, 1))
            %   mu        - First-derivative penalty weight (must be >= 0)
            %   maxIter   - Maximum IRLS iterations (must be positive integer)
            %   tolerance - Convergence tolerance (must be > 0)
            %
            % Outputs:
            %   isValid - Logical, true if all parameters are valid
            %   msg     - Character vector describing the issue (empty if valid)

            isValid = true;
            msg = '';

            if ~isscalar(lambda) || ~isnumeric(lambda) || lambda < 0
                isValid = false;
                msg = 'Lambda must be a non-negative scalar.';
                return;
            end

            if ~isscalar(p) || ~isnumeric(p) || p < 0 || p > 1
                isValid = false;
                msg = 'p must be a scalar in [0, 1].';
                return;
            end

            if ~isscalar(mu) || ~isnumeric(mu) || mu < 0
                isValid = false;
                msg = 'mu must be a non-negative scalar.';
                return;
            end

            if ~isscalar(maxIter) || ~isnumeric(maxIter) || maxIter < 1 || floor(maxIter) ~= maxIter
                isValid = false;
                msg = 'Max iterations must be a positive integer.';
                return;
            end

            if ~isscalar(tolerance) || ~isnumeric(tolerance) || tolerance <= 0
                isValid = false;
                msg = 'Tolerance must be a positive scalar.';
                return;
            end
        end

        function [isValid, msg] = validateInterval(startIdx, endIdx, nChannels)
            % validateInterval - Validate an interval definition.
            %
            % Inputs:
            %   startIdx  - Start index of the interval
            %   endIdx    - End index of the interval
            %   nChannels - Total number of spectral channels
            %
            % Outputs:
            %   isValid - Logical, true if interval is valid
            %   msg     - Character vector describing the issue (empty if valid)

            isValid = true;
            msg = '';

            if ~isscalar(startIdx) || ~isnumeric(startIdx) || startIdx < 1
                isValid = false;
                msg = 'Start index must be a positive integer.';
                return;
            end

            if ~isscalar(endIdx) || ~isnumeric(endIdx) || endIdx < 1
                isValid = false;
                msg = 'End index must be a positive integer.';
                return;
            end

            if startIdx >= endIdx
                isValid = false;
                msg = 'Start index must be less than end index.';
                return;
            end

            if endIdx > nChannels
                isValid = false;
                msg = sprintf('End index (%d) exceeds number of channels (%d).', endIdx, nChannels);
                return;
            end

            if (endIdx - startIdx) < 3
                isValid = false;
                msg = 'Interval must span at least 3 channels.';
                return;
            end
        end

        function [isValid, msg] = validateIntervalParams(lambda, p)
            % validateIntervalParams - Validate per-interval lambda and p values.
            %
            % Inputs:
            %   lambda - Local lambda for the interval
            %   p      - Local p for the interval
            %
            % Outputs:
            %   isValid - Logical, true if parameters are valid
            %   msg     - Character vector describing the issue (empty if valid)

            isValid = true;
            msg = '';

            if ~isscalar(lambda) || ~isnumeric(lambda) || lambda < 0
                isValid = false;
                msg = 'Interval lambda must be a non-negative scalar.';
                return;
            end

            if ~isscalar(p) || ~isnumeric(p) || p < 0 || p > 1
                isValid = false;
                msg = 'Interval p must be a scalar in [0, 1].';
                return;
            end
        end
    end
end
