function imshow(img, varargin)
% imshow - Display image in axes (vanilla MATLAB implementation)
%
% SYNTAX:
%   imshow(img)
%   imshow(img, 'Parent', ax)
%
% INPUTS:
%   img - Image array (grayscale or RGB, [0,1] or uint8)
%   'Parent' - Optional parent axes
%
% DESCRIPTION:
%   Vanilla MATLAB implementation of imshow to avoid Image Processing
%   Toolbox dependency. Uses image() with proper display settings.
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-08
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0

% Parse inputs
ax = [];
for i = 1:2:length(varargin)
    if strcmpi(varargin{i}, 'Parent')
        ax = varargin{i+1};
    end
end

% Use current axes if not specified
if isempty(ax)
    ax = gca;
end

% Display image
image(ax, img);
axis(ax, 'image', 'off');
colormap(ax, gray(256));

% Set proper display range
if size(img, 3) == 1
    % Grayscale
    if isa(img, 'double') || isa(img, 'single')
        clim(ax, [0 1]);
    else
        clim(ax, [0 255]);
    end
end

end
