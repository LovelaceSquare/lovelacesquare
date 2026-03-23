function [channels, options] = loadPreset(filename)
% loadPreset - Load channel settings and options from file
%
% SYNTAX:
%   [channels, options] = loadPreset(filename)
%
% INPUTS:
%   filename - Path to preset file (.mat)
%
% OUTPUTS:
%   channels - Struct array of channel settings
%   options  - Struct with global options
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-08
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0

data = load(filename, '-mat');

if ~isfield(data, 'preset')
    error('Invalid preset file');
end

preset = data.preset;

channels = preset.channels;
options = preset.options;

end
