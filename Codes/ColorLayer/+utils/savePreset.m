function savePreset(filename, channels, options)
% savePreset - Save channel settings and options to file
%
% SYNTAX:
%   savePreset(filename, channels, options)
%
% INPUTS:
%   filename - Path to save preset (.mat file)
%   channels - Struct array of channel settings
%   options  - Struct with global options
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-08
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0

preset = struct();
preset.channels = channels;
preset.options = options;
preset.version = '1.0';
preset.dateCreated = datestr(now);

save(filename, 'preset', '-mat');

end
