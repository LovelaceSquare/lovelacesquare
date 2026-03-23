function mdisResized = resizeMaps(mdis, scaleFactor)
% resizeMaps - Resize all maps by a global scaling factor
%
% SYNTAX:
%   mdisResized = resizeMaps(mdis, scaleFactor)
%
% INPUTS:
%   mdis        - Cell array of 2D maps
%   scaleFactor - Scaling factor (1 = no change, 0.5 = half size, 2 = double)
%
% OUTPUTS:
%   mdisResized - Cell array of resized maps
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-08
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0

if scaleFactor == 1
    mdisResized = mdis;
    return;
end

mdisResized = cell(size(mdis));

% Handle both mdis{i} and mdis{1,i} indexing
if size(mdis, 1) == 1
    for i = 1:length(mdis)
        mdisResized{1,i} = utils.imresize(mdis{1,i}, scaleFactor, 'bicubic');
    end
else
    for i = 1:length(mdis)
        mdisResized{i} = utils.imresize(mdis{i}, scaleFactor, 'bicubic');
    end
end

end
