function dataBinned = binning(data, binVector, mode)
% binning Bins the data array according to specified bin sizes and mode.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 2025-04-04
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: v 2.0
%
% This function performs binning on an N-dimensional array. The binning
% process involves grouping adjacent elements specified by `binVector`
% along each dimension. Once the elements are grouped into bins, the
% function either sums or averages the values within each bin, depending
% on the selected `mode`.
%
% -- Optimized Approach --
% We minimize dimension permutations by permuting only the dimension we're
% currently binning to the leading dimension, summing/averaging, and then
% permuting back. This prevents rotating all dimensions repeatedly.
%
% INPUTS:
%    data (array): The input N-dimensional array to be binned.
%    binVector (vector): A vector specifying the bin sizes for each
%                        dimension of the input array. Each element of
%                        binVector corresponds to how many elements along
%                        that dimension will be grouped into a single bin.
%    mode (string): The mode of binning. If set to 'sum', the elements in
%                   each bin are summed. If set to 'mean', the elements
%                   are averaged.
%
% OUTPUT:
%    dataBinned (array): The resulting binned data array. Its dimensions
%                        are reduced since multiple original data points
%                        are combined into single binned values.
%
% NOTE:
%   - This function truncates the data along each dimension so that it
%     evenly fits into the specified bin sizes. If bin sizes do not neatly
%     divide the dimension sizes, the data is reduced accordingly.
%   - For large N-D arrays, this approach avoids repeated dimension
%     rotations at each step, making it typically faster in practice.
%
%
% Version History:
%    v 1.0 - Original implementation using rotating permutation method.
%    v 2.0 - Optimized version using selective dimension permutation for speed.
%
%
% EXAMPLE:
%    dataBinned = binning(rand(100, 200), [4 5], 'sum');
%
% Disclaimer:
%   Authors and Lovelace's Square are not responsible for any issues,
%   inaccuracies, or data loss arising from the use of this function.

% -------------------------------------------------------------------------
% 1) Initial checks and truncation
% -------------------------------------------------------------------------
sizes = size(data);
numDims = ndims(data);

if length(binVector) ~= numDims
    error('binVector must have the same number of dimensions as the data.');
end

% Compute truncated sizes so that each dimension is a multiple of the bin size
maxSizes             = floor(sizes ./ binVector);
AdjustedOriginalSizes = maxSizes .* binVector;

if any(AdjustedOriginalSizes == 0)
    error('The specified binning is too large for one or more dimensions.');
end
if any(binVector == 0)
    error('binVector must not contain zeros.');
end

% Truncate data so each dimension fits evenly
subCell = cell(1, numDims);
for d = 1:numDims
    subCell{d} = 1:AdjustedOriginalSizes(d);
end
data = data(subCell{:});

% -------------------------------------------------------------------------
% 2) Perform binning dimension by dimension, minimizing permutations
% -------------------------------------------------------------------------
for dimIndex = 1:numDims
    % Current dimension we are binning along is dimIndex

    % Permute 'dimIndex' to the first dimension so we can bin easily by reshape
    permOrder = 1:numDims;
    permOrder(1) = dimIndex;
    permOrder(dimIndex) = 1;
    data = permute(data, permOrder);

    % Now data's size has its "binning dimension" in dim = 1
    % Next, we bin along dimension 1 in groups of binVector(dimIndex)

    currSize = size(data);  % after permutation
    groupSize = binVector(dimIndex);

    % Reshape so that dimension 1 is grouped in blocks of 'groupSize'
    % That is: data is (groupSize x newLen x other dims)
    % We flatten everything except the first dimension
    data = reshape(data, groupSize, []);

    % Sum or mean across the first dimension
    switch mode
        case 'sum'
            data = sum(data, 1);
        case 'mean'
            data = mean(data, 1);
        otherwise
            error('Mode must be either "sum" or "mean".');
    end

    % data is now [1 x newLen], but we must restore it into the shape
    % where dimension 1 is "currSize(1)/groupSize" => newBins
    newBins = currSize(1) / groupSize;
    newShape = [newBins, currSize(2:end)];
    data = reshape(data, newShape);

    % Permute back to the original dimension arrangement
    data = ipermute(data, permOrder);

    % Update the truncated size for that dimension
    sizes(dimIndex) = sizes(dimIndex) / binVector(dimIndex);
end

% The result is fully binned now
dataBinned = data;

end
