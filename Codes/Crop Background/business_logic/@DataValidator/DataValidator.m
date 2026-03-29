classdef DataValidator
% DataValidator. Static helpers that validate inputs for CropBackgroundModule.
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-16
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0

    methods (Static)

        function [valid, msg] = validate3DArray(data)
        % validate3DArray  Ensure the input is a 3-D numeric array.
        %
        %   [valid, msg] = DataValidator.validate3DArray(data)

            valid = true;
            msg   = '';

            if ~isnumeric(data)
                valid = false;
                msg   = 'Input must be a numeric array.';
                return;
            end

            if ndims(data) ~= 3 %#ok<ISMAT>
                valid = false;
                msg   = 'Input must be a 3-D array [rows x cols x channels].';
                return;
            end

            if any(size(data) == 0)
                valid = false;
                msg   = 'Input array must not have any zero-length dimension.';
                return;
            end

            if any(isnan(data(:))) || any(isinf(data(:)))
                valid = false;
                msg   = 'Input array must not contain NaN or Inf values.';
                return;
            end
        end

        function [valid, msg] = validateThresholds(minT, maxT)
        % validateThresholds  Ensure thresholds are valid numbers.
        %
        %   [valid, msg] = DataValidator.validateThresholds(minT, maxT)

            valid = true;
            msg   = '';

            if ~isnumeric(minT) || ~isscalar(minT) || isnan(minT)
                valid = false;
                msg   = 'Min threshold must be a numeric scalar.';
                return;
            end

            if ~isnumeric(maxT) || ~isscalar(maxT) || isnan(maxT)
                valid = false;
                msg   = 'Max threshold must be a numeric scalar.';
                return;
            end

            if minT > maxT
                valid = false;
                msg   = 'Min threshold must be less than or equal to max threshold.';
                return;
            end
        end

    end
end
