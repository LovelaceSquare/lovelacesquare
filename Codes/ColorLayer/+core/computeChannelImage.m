function channelImage = computeChannelImage(mdis, channel, options)
% computeChannelImage - Compute processed single-channel image for preview
%
% SYNTAX:
%   channelImage = computeChannelImage(mdis, channel, options)
%
% INPUTS:
%   mdis    - Cell array of 2D maps
%   channel - Single channel struct with fields:
%             .mapIndex, .brightness, .saturation, .lumination
%   options - Struct with fields:
%             .perMapNormalize - true/false
%
% OUTPUTS:
%   channelImage - 2D image [0,1] double
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-08
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0

if nargin < 3
    options = struct('perMapNormalize', false);
end

% Extract map (handle both mdis{i} and mdis{1,i} indexing)
if channel.mapIndex < 1 || channel.mapIndex > length(mdis)
    if size(mdis, 1) == 1
        channelImage = zeros(size(mdis{1,1}));
    else
        channelImage = zeros(size(mdis{1}));
    end
    return;
end

if size(mdis, 1) == 1
    M = double(mdis{1, channel.mapIndex});
else
    M = double(mdis{channel.mapIndex});
end

% Per-map normalization (if enabled)
if options.perMapNormalize
    M = utils.mat2gray(M);
end

% Apply saturation clipping
if isfinite(channel.saturation)
    M(M > channel.saturation) = channel.saturation;
end

% Normalize to [0, 1] if not already done
if ~options.perMapNormalize
    M = utils.mat2gray(M);
end

% Apply brightness and lumination
M = M * channel.brightness * channel.lumination;

% Clamp to [0, 1]
channelImage = max(0, min(1, M));

end
