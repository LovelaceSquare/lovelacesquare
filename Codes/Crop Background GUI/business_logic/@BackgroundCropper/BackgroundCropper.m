classdef BackgroundCropper
% BackgroundCropper. Core routines for background-pixel cropping.
%
%   Provides four stateless methods:
%     computeIntensity    - sum across spectral channels
%     applyThreshold      - classify pixels as retained / discarded
%     extractPixels       - build the 2-D cropped matrix
%     preprocessIntensity - apply contrast-enhancing transforms
%
% Author: Adrian Gomez-Sanchez
% Date Created: 2026-03-16
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.1

    methods

        function globalIntensity = computeIntensity(~, imageCube)
        % computeIntensity  Sum image cube along the third dimension.
        %
        %   globalIntensity = computeIntensity(obj, imageCube)
        %
        %   Input
        %     imageCube       [rows x cols x channels] numeric array
        %
        %   Output
        %     globalIntensity [rows x cols] total intensity per pixel

            globalIntensity = sum(imageCube, 3);
        end

        function [retainedIdx, discardedIdx] = applyThreshold(~, globalIntensity, minThresh, maxThresh)
        % applyThreshold  Separate pixel indices by intensity thresholds.
        %
        %   [retainedIdx, discardedIdx] = applyThreshold(obj, globalIntensity, minThresh, maxThresh)
        %
        %   Pixels whose total intensity falls outside [minThresh, maxThresh]
        %   are marked as discarded.

            allPixels     = globalIntensity(:);
            discardedMask = (allPixels < minThresh) | (allPixels > maxThresh);
            retainedIdx   = find(~discardedMask);
            discardedIdx  = find(discardedMask);
        end

        function croppedMatrix = extractPixels(~, imageCube, retainedIdx)
        % extractPixels  Extract retained pixels into a 2-D matrix.
        %
        %   croppedMatrix = extractPixels(obj, imageCube, retainedIdx)
        %
        %   Output
        %     croppedMatrix  [nRetainedPixels x nChannels]

            [nRows, nCols, nChannels] = size(imageCube);
            reshaped      = reshape(imageCube, nRows * nCols, nChannels);
            croppedMatrix = reshaped(retainedIdx, :);
        end

        function processedIntensity = preprocessIntensity(~, rawIntensity, method)
        % preprocessIntensity  Apply a transform to enhance contrast.
        %
        %   processedIntensity = preprocessIntensity(obj, rawIntensity, method)
        %
        %   method: 'none' | 'square' | 'log' | 'sqrt'

            switch lower(method)
                case 'square'
                    processedIntensity = rawIntensity .^ 2;
                case 'log'
                    processedIntensity = log1p(abs(rawIntensity));
                case 'sqrt'
                    processedIntensity = sqrt(abs(rawIntensity));
                otherwise
                    processedIntensity = rawIntensity;
            end
        end

    end
end
