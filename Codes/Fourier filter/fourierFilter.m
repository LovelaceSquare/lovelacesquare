function [filteredData, freqVector] = fourierFilter(data, freqIntervals, dt)
% fourierFilter. Apply frequency-domain filtering to time-series signals.
%
% REFERENCES:
%   The algorithm is based on standard FFT processing techniques to 
%   selectively retain specified frequency bands.
%
% Author:  Adrián Gómez-Sánchez
% Date:    2024-12-14
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 1.0
%
% The fourierFilter function processes 2D time-series data by applying a 
% frequency-domain filter. It computes the FFT of each row (signal), and 
% retains only the frequency components within the specified ranges. The 
% inverse FFT is then used to reconstruct the filtered time-domain signals.
%
% Specifically, the function:
%   1) Computes the FFT (with zero-padding to the next power of 2) for each
%      row of the input data.
%   2) Constructs a frequency mask based on the provided frequency intervals.
%   3) Applies the mask symmetrically to retain both positive and negative 
%      frequency components, ensuring real-valued output.
%   4) Performs an inverse FFT to reconstruct the filtered signals.
%
% INPUTS:
%    data (array): A 2D numeric matrix of size [nSamples x nTimePoints] 
%       containing time-series data. Each row represents a separate signal.
%    freqIntervals (array): An N-by-2 numeric matrix specifying frequency 
%       ranges (in Hz) to retain, e.g., [0 10; 20 30] retains components 
%       between 0–10 Hz and 20–30 Hz. If empty ([]), no filtering is performed 
%       (i.e., the output equals the input).
%    dt (double): The sampling interval in seconds. Default = 1. The sampling
%       frequency is computed as fs = 1/dt.
%
% OUTPUTS:
%    filteredData (array): The filtered time-series data, having the same
%       dimensions as the input data.
%    freqVector (array): A vector containing the frequency axis (in Hz) for 
%       the positive half of the FFT (length = nFFT/2+1).
%
% EXAMPLE:
%    % Suppose 'X' is a 10x1000 matrix, where each row is a time-series signal.
%    % Retain frequencies in the ranges 0–15 Hz and 40–60 Hz, using a sampling 
%    % interval of dt = 0.01 (fs = 100 Hz):
%       intervals = [0 15; 40 60];
%       [Xfiltered, fVec] = fourierFilter(X, intervals, 0.01);
%    % Xfiltered now contains the filtered signals, and fVec is the frequency axis.
%
% DISCLAIMER:
%    Authors and Lovelace's Square are not responsible for any issues, 
%    inaccuracies, or data loss arising from the use of this function.

    % -------------------
    % 1) INPUT VALIDATION
    % -------------------
    if nargin < 2
        error('Usage: fourierFilter(data, freqIntervals, [dt]).');
    end
    if ~isnumeric(data) || ndims(data) ~= 2 || isempty(data)
        error('Input "data" must be a non-empty 2D numeric matrix [samples x timePoints].');
    end
    [nSamples, nTime] = size(data);
    
    if ~isreal(data)
        warning('Complex data detected. Using only the real part.');
        data = real(data);
    end
    
    if nargin < 3 || isempty(dt)
        dt = 1;  % default sampling interval
    end
    
    if ~isempty(freqIntervals) && (~isnumeric(freqIntervals) || size(freqIntervals,2) ~= 2)
        error('freqIntervals must be an N-by-2 numeric array or empty ([]).');
    end
    
    % Merge overlapping intervals if provided
    if ~isempty(freqIntervals)
        freqIntervals = sortrows(freqIntervals);
        freqIntervals = mergeOverlappingRanges(freqIntervals);
    end
    
    % If no intervals are specified, return the input as output.
    if isempty(freqIntervals)
        filteredData = data;
        fs = 1/dt;
        nFFT = 2^nextpow2(nTime);
        freqVector = fs/2 * linspace(0, 1, nFFT/2+1);
        return;
    end
    
    % -------------------------
    % 2) PREPARE FFT PARAMETERS
    % -------------------------
    fs = 1/dt;  % Sampling frequency
    nFFT = 2^nextpow2(nTime);  % Zero-pad to the next power of 2
    freqVector = fs/2 * linspace(0, 1, nFFT/2+1);  % Frequency axis for positive frequencies

    % Build a frequency mask for positive frequencies
    freqMask = false(1, length(freqVector));
    for i = 1:size(freqIntervals,1)
        lo = max(freqIntervals(i,1), 0);
        hi = min(freqIntervals(i,2), fs/2);
        if lo > hi
            continue;  % Skip if no valid frequency range
        end
        freqMask = freqMask | (freqVector >= lo & freqVector <= hi);
    end
    
    dcIndex = 1;           % DC component at 0 Hz
    nyquistIndex = nFFT/2+1; % Nyquist frequency (if applicable)
    allSelected = find(freqMask);
    
    % --------------
    % 3) FILTER DATA
    % --------------
    filteredData = zeros(size(data));
    for row = 1:nSamples
        % Compute the FFT for the current row with zero-padding.
        X = fft(data(row,:), nFFT);
        
        % Create an array to hold the frequency components to retain.
        keepFFT = zeros(size(X));
        
        % Retain DC component if selected.
        if any(allSelected == dcIndex)
            keepFFT(dcIndex) = X(dcIndex);
        end
        
        % Retain Nyquist frequency if nFFT is even and selected.
        if mod(nFFT,2)==0 && any(allSelected == nyquistIndex)
            keepFFT(nyquistIndex) = X(nyquistIndex);
        end
        
        % Retain positive frequency components within the specified intervals.
        corePosSel = allSelected(allSelected > dcIndex & allSelected < nyquistIndex);
        negSel = nFFT - (corePosSel - 2);  % Corresponding negative frequencies
        keepFFT(corePosSel) = X(corePosSel);
        keepFFT(negSel) = X(negSel);
        
        % Reconstruct the time-domain signal using the inverse FFT (ensuring real output).
        yRec = ifft(keepFFT, 'symmetric');
        
        % Truncate the signal to its original length.
        filteredData(row,:) = yRec(1:nTime);
    end
end

%--------------------------------------------------------------------------
function mergedRanges = mergeOverlappingRanges(ranges)
% mergeOverlappingRanges. Merge overlapping or adjacent frequency intervals.
%
% This helper function takes an N-by-2 array of frequency ranges sorted 
% in ascending order (by the first column) and merges any that overlap or 
% are adjacent.
%
% INPUT:
%    ranges (array): An N-by-2 numeric array of frequency intervals.
%
% OUTPUT:
%    mergedRanges (array): An M-by-2 array of merged frequency intervals.
%
% EXAMPLE:
%    mergeOverlappingRanges([0 10; 8 15; 20 25]) 
%    returns: [0 15; 20 25]

    if isempty(ranges)
        mergedRanges = [];
        return;
    end
    
    mergedRanges = [];
    currentStart = ranges(1,1);
    currentEnd = ranges(1,2);
    
    for i = 2:size(ranges,1)
        if ranges(i,1) <= currentEnd
            currentEnd = max(currentEnd, ranges(i,2));
        else
            mergedRanges = [mergedRanges; currentStart, currentEnd]; 
            currentStart = ranges(i,1);
            currentEnd = ranges(i,2);
        end
    end
    mergedRanges = [mergedRanges; currentStart, currentEnd];
end
