function [isValid, msg] = validateMaps(mdis)
% validateMaps - Validate that mdis is a proper cell array of 2D maps
%
% SYNTAX:
%   [isValid, msg] = validateMaps(mdis)
%
% INPUTS:
%   mdis - Cell array to validate
%
% OUTPUTS:
%   isValid - true if valid, false otherwise
%   msg     - Error/warning message
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-08
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0

isValid = true;
msg = '';

% Check if cell array
if ~iscell(mdis)
    isValid = false;
    msg = 'mdis must be a cell array';
    return;
end

% Check if empty
if isempty(mdis)
    isValid = false;
    msg = 'mdis is empty';
    return;
end

% Get reference size from first map (handle both mdis{i} and mdis{1,i} indexing)
if size(mdis, 1) == 1
    refSize = size(mdis{1,1});
else
    refSize = size(mdis{1});
end

if length(refSize) ~= 2
    isValid = false;
    msg = 'Maps must be 2D matrices';
    return;
end

% Check all maps (handle both indexing styles)
if size(mdis, 1) == 1
    for i = 1:length(mdis)
        if ~isnumeric(mdis{1,i})
            isValid = false;
            msg = sprintf('Map %d is not numeric', i);
            return;
        end

        if ~isequal(size(mdis{1,i}), refSize)
            isValid = false;
            msg = sprintf('Map %d has different size than map 1', i);
            return;
        end
    end
else
    for i = 1:length(mdis)
        if ~isnumeric(mdis{i})
            isValid = false;
            msg = sprintf('Map %d is not numeric', i);
            return;
        end

        if ~isequal(size(mdis{i}), refSize)
            isValid = false;
            msg = sprintf('Map %d has different size than map 1', i);
            return;
        end
    end
end

msg = sprintf('Valid: %d maps of size %dx%d', length(mdis), refSize(1), refSize(2));

end
