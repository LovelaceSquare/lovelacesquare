%% demo_ColorLayerGUI.m
% Demo script to launch ColorLayer Editor GUI with the example composite image.
%
% Loads the bundled example image (a false-color fluorescence microscopy
% image), splits its RGB channels into separate maps, and opens the GUI.
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-08
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0

clear; clc; close all;

% Add ColorLayer to path (works from any directory)
thisFile = mfilename('fullpath');
colorLayerPath = fileparts(fileparts(thisFile));
addpath(colorLayerPath);

%% Load the example composite image
imgPath = fullfile(colorLayerPath, 'demos', 'example_composite.png');

if ~isfile(imgPath)
    error('Example image not found: %s', imgPath);
end

fprintf('Loading example image: %s\n', imgPath);
img = imread(imgPath);
fprintf('Image size: %d x %d x %d (%s)\n', size(img, 1), size(img, 2), size(img, 3), class(img));

%% Launch the GUI
% The constructor accepts H x W x 3 numeric arrays directly.
% It splits the channels into separate maps automatically.
fprintf('Launching ColorLayer GUI...\n');
ColorLayerGUI(img);
