function [croppedMatrix, retainedPixelsIdx, discardedPixelsIdx] = cropBackground(imageCube, minThresh, maxThresh)
% CROPBACKGROUND Crop background pixels from a 3D image cube by intensity range.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 2024-12-14
% License: MIT
% Reviewed by Lovelace's Square: Yes
%
% DESCRIPTION:
% This function processes a 3D image cube by summing the pixel intensities 
% across all channels to compute a global intensity map, then filters out 
% pixels whose total intensity falls outside a user-defined range specified 
% by minimum and maximum thresholds. It flattens the 2D intensity map into 
% a vector, identifies the indices of pixels that meet the threshold criteria, 
% and constructs a cropped matrix containing the multi-channel intensity 
% values for just those pixels. Additionally, it outputs both the indices 
% of the retained pixels and those that were discarded, thereby effectively 
% removing background pixels from the image data.
%
% USAGE:
%   [croppedMatrix, retainedPixelsIdx, discardedPixelsIdx] = ...
%       cropBackground(imageCube, minThresh, maxThresh)
%
% INPUTS:
%   imageCube  - 3D numeric array of size [rows x cols x channels].
%   minThresh  - Lower limit on pixel intensity (scalar).
%   maxThresh  - Upper limit on pixel intensity (scalar).
%                Pixels with global intensity < minThresh or > maxThresh are discarded.
%
% OUTPUTS:
%   croppedMatrix     - 2D numeric array of size [NretainedPixels x channels], 
%                       containing intensities of retained pixels for each channel.
%   retainedPixelsIdx - Linear indices (relative to the original rows*cols layout)
%                       of the retained pixels.
%   discardedPixelsIdx- Linear indices of the discarded pixels.
%
% EXAMPLE:
%   % Suppose you have a hyperspectral cube 'imgCube' of size [100 x 80 x 50].
%   % You want to keep pixels with total intensity between 100 and 500.
%   minT = 100; 
%   maxT = 500;
%   [croppedData, keptIdx, removedIdx] = cropBackground(imgCube, minT, maxT);
%
%   % 'croppedData' will be (Nretained x 50). 
%   % 'keptIdx' has the linear indices of those kept pixels in the original 2D plane.
%
% DISCLAIMER:
%   Authors and Lovelace's Square are not responsible for any issues, 
%   inaccuracies, or data loss arising from the use of this function.

    %% 1) Validate input
    if ndims(imageCube) ~= 3
        error('Input imageCube must be a 3D array [rows x cols x channels].');
    end
    
    [rows, cols, nChannels] = size(imageCube);

    if nargin < 3
        error(['Usage: cropBackground(imageCube, minThresh, maxThresh)\n',...
               'You must specify both minThresh and maxThresh.']);
    end
    if ~isnumeric(imageCube) || ~isnumeric(minThresh) || ~isnumeric(maxThresh)
        error('imageCube, minThresh, and maxThresh must be numeric.');
    end
    if minThresh > maxThresh
        error('minThresh must be less than or equal to maxThresh.');
    end

    %% 2) Compute global intensity by summing across channels
    %    This yields a 2D array of size [rows x cols].
    globalIntensity = sum(imageCube, 3);

    %% 3) Flatten the intensity matrix into a vector of pixel intensities
    allPixels = globalIntensity(:);  % size [rows*cols x 1]
    
    %% 4) Determine which pixels fall within the specified intensity range
    inRangeMask = (allPixels >= minThresh) & (allPixels <= maxThresh);
    retainedPixelsIdx  = find(inRangeMask);      % linear indices of retained pixels
    discardedPixelsIdx = find(~inRangeMask);     % linear indices of discarded pixels
    
    %% 5) Create the cropped matrix of retained pixel intensities
    %    We'll have (#retainedPixels) x (nChannels).
    nRetained = numel(retainedPixelsIdx);
    croppedMatrix = zeros(nRetained, nChannels, 'like', imageCube);  % preserve data type
    
    for c = 1:nChannels
        channelData = imageCube(:,:,c);
        croppedMatrix(:, c) = channelData(retainedPixelsIdx);
    end

end
