classdef DataValidator
% DataValidator. Validation utilities for AsLS_GUI.
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-16
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 2.0
%
% Validates input data matrices and correction parameters before
% processing. All methods are static.
%
% METHODS (Static):
%   validateData   - Validate input data matrix (numeric, 2D, finite)
%   validateParams - Validate lambda, p, maxIter parameters

    methods (Static)
        function [isValid, msg] = validateData(data)
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

            if any(isnan(data(:))) || any(isinf(data(:)))
                isValid = false;
                msg = 'Data contains NaN or Inf values.';
                return;
            end

            if size(data, 2) < 4
                isValid = false;
                msg = 'Data must contain at least 4 channels (columns).';
                return;
            end
        end

        function [isValid, msg] = validateParams(lambda, p, smoothEnabled, smoothLambda, maxIter)
            isValid = true;
            msg = '';

            if ~isscalar(lambda) || ~isnumeric(lambda) || ~isfinite(lambda) || lambda <= 0
                isValid = false;
                msg = 'Lambda must be a positive scalar.';
                return;
            end

            if ~isscalar(p) || ~isnumeric(p) || ~isfinite(p) || p < 0 || p > 1
                isValid = false;
                msg = 'p must be in [0, 1].';
                return;
            end

            if ~isscalar(smoothEnabled) || ~(islogical(smoothEnabled) || isnumeric(smoothEnabled))
                isValid = false;
                msg = 'Smooth toggle must be logical.';
                return;
            end

            if logical(smoothEnabled)
                if ~isscalar(smoothLambda) || ~isnumeric(smoothLambda) || ~isfinite(smoothLambda) || smoothLambda <= 0
                    isValid = false;
                    msg = 'Smooth lambda must be a positive scalar when smoothing is enabled.';
                    return;
                end
            end

            if ~isscalar(maxIter) || ~isnumeric(maxIter) || ~isfinite(maxIter) || maxIter < 1 || floor(maxIter) ~= maxIter
                isValid = false;
                msg = 'Max iterations must be a positive integer.';
                return;
            end
        end

        function [isValid, msg] = validateIntervals(intervals, nChannels)
            isValid = true;
            msg = '';

            if isempty(intervals)
                return;
            end

            if ~isstruct(intervals)
                isValid = false;
                msg = 'Intervals must be a struct array.';
                return;
            end

            requiredFields = {'start', 'end', 'lambda', 'p'};
            for k = 1:numel(requiredFields)
                if ~isfield(intervals, requiredFields{k})
                    isValid = false;
                    msg = sprintf('Interval field missing: %s', requiredFields{k});
                    return;
                end
            end

            for i = 1:numel(intervals)
                s = intervals(i).start;
                e = intervals(i).end;
                lam = intervals(i).lambda;
                p = intervals(i).p;

                if ~isscalar(s) || ~isnumeric(s) || ~isfinite(s)
                    isValid = false;
                    msg = sprintf('Interval %d has invalid start.', i);
                    return;
                end
                if ~isscalar(e) || ~isnumeric(e) || ~isfinite(e)
                    isValid = false;
                    msg = sprintf('Interval %d has invalid end.', i);
                    return;
                end
                if s < 1 || e < 1 || s > nChannels || e > nChannels
                    isValid = false;
                    msg = sprintf('Interval %d is outside channel range.', i);
                    return;
                end
                if e <= s
                    isValid = false;
                    msg = sprintf('Interval %d end must be greater than start.', i);
                    return;
                end
                if (e - s) < 1
                    isValid = false;
                    msg = sprintf('Interval %d is too narrow.', i);
                    return;
                end

                if ~isscalar(lam) || ~isnumeric(lam) || ~isfinite(lam) || lam <= 0
                    isValid = false;
                    msg = sprintf('Interval %d lambda must be positive.', i);
                    return;
                end
                if ~isscalar(p) || ~isnumeric(p) || ~isfinite(p) || p < 0 || p > 1
                    isValid = false;
                    msg = sprintf('Interval %d p must be in [0, 1].', i);
                    return;
                end
            end
        end
    end
end
