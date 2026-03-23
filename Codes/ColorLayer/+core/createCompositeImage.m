function compositeImage = createCompositeImage(mdis, channels, options)
% createCompositeImage - Create false-color composite from maps and channel settings
%
% SYNTAX:
%   compositeImage = createCompositeImage(mdis, channels, options)
%
% INPUTS:
%   mdis     - Cell array of 2D maps {map1, map2, ..., mapN}
%   channels - Struct array (6 elements) with fields:
%              .name       - 'Red', 'Green', 'Blue', 'Yellow', 'Cyan', 'Magenta'
%              .enabled    - true/false
%              .mapIndex   - Index into mdis (1..N)
%              .brightness - Multiplicative factor (default: 1.0)
%              .saturation - Max clipping value (default: Inf)
%              .lumination - Overall intensity scale (default: 1.0)
%              .rgb        - RGB triplet [R G B] for this channel
%   options  - Struct with fields:
%              .normalizationMode - 'None', 'Global mat2gray',
%                                   'Per-channel normalize', 'Per-map normalize'
%              .perMapNormalize   - true/false (normalize maps before compositing)
%
% OUTPUTS:
%   compositeImage - RGB image [0,1] double (H x W x 3)
%
% DESCRIPTION:
%   Creates a false-color composite by:
%   1. Optionally normalizing individual maps (if perMapNormalize = true)
%   2. For each enabled channel:
%      - Extract the selected map
%      - Apply saturation clipping
%      - Normalize to [0,1]
%      - Apply brightness and lumination
%      - Add to R, G, B channels based on RGB triplet
%   3. Apply final normalization based on normalizationMode
%
% CHANNEL COLOR MAPPING (Additive):
%   Red     [1 0 0] → R channel
%   Green   [0 1 0] → G channel
%   Blue    [0 0 1] → B channel
%   Yellow  [1 1 0] → R + G channels
%   Cyan    [0 1 1] → G + B channels
%   Magenta [1 0 1] → R + B channels
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-08
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0

%% Validate inputs
if isempty(mdis) || ~iscell(mdis)
    error('mdis must be a non-empty cell array of 2D maps');
end

if nargin < 3
    options = struct();
end

% Default options
if ~isfield(options, 'normalizationMode')
    options.normalizationMode = 'Global mat2gray';
end
if ~isfield(options, 'perMapNormalize')
    options.perMapNormalize = false;
end

%% Determine image size from first map
% Handle both mdis{i} and mdis{1,i} indexing
if size(mdis, 1) == 1
    firstMap = mdis{1,1};
else
    firstMap = mdis{1};
end
[H, W] = size(firstMap);

%% Initialize RGB channels
R = zeros(H, W);
G = zeros(H, W);
B = zeros(H, W);

%% Process each channel
for i = 1:length(channels)
    ch = channels(i);

    % Skip if disabled or invalid map index
    if ~ch.enabled || ch.mapIndex < 1 || ch.mapIndex > length(mdis)
        continue;
    end

    % Extract map (handle both indexing styles)
    if size(mdis, 1) == 1
        M = mdis{1, ch.mapIndex};
    else
        M = mdis{ch.mapIndex};
    end

    % Ensure same size
    if ~isequal(size(M), [H, W])
        warning('Map %d size mismatch. Resizing...', ch.mapIndex);
        M = utils.imresize(M, [H, W], 'bicubic');
    end

    % Convert to double
    M = double(M);

    % Per-map normalization (if enabled)
    if options.perMapNormalize
        M = utils.mat2gray(M);
    end

    % Apply saturation clipping
    if isfinite(ch.saturation)
        M(M > ch.saturation) = ch.saturation;
    end

    % Normalize to [0, 1] if not already done
    if ~options.perMapNormalize
        M = utils.mat2gray(M);
    end

    % Apply brightness and lumination
    M = M * ch.brightness * ch.lumination;

    % Add to RGB channels based on channel's RGB triplet
    R = R + M * ch.rgb(1);
    G = G + M * ch.rgb(2);
    B = B + M * ch.rgb(3);
end

%% Apply final normalization
switch options.normalizationMode
    case 'None'
        % Just clamp to [0, 1]
        R = max(0, min(1, R));
        G = max(0, min(1, G));
        B = max(0, min(1, B));

    case 'Global mat2gray'
        % Normalize entire composite as a whole
        composite = cat(3, R, G, B);
        composite = utils.mat2gray(composite);
        R = composite(:,:,1);
        G = composite(:,:,2);
        B = composite(:,:,3);

    case 'Per-channel normalize'
        % Normalize R, G, B independently
        R = utils.mat2gray(R);
        G = utils.mat2gray(G);
        B = utils.mat2gray(B);

    case 'Per-map normalize'
        % Already handled above, just clamp
        R = max(0, min(1, R));
        G = max(0, min(1, G));
        B = max(0, min(1, B));

    otherwise
        warning('Unknown normalization mode: %s. Using None.', options.normalizationMode);
        R = max(0, min(1, R));
        G = max(0, min(1, G));
        B = max(0, min(1, B));
end

%% Assemble final composite
compositeImage = cat(3, R, G, B);

end
