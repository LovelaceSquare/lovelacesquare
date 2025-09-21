function [D_kernelized] = kernelize(D, num_kernels, kernel_width)
% kernelize Apply kernelization to the input matrix D
% Adapted from: Gómez-Sánchez, Adrián, et al. "Kernelizing: 
% A way to increase accuracy in trilinear decomposition analysis 
% of multiexponential signals." Analytica Chimica Acta 1273 (2023): 341545.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 2024-12-14
% License: MIT
% Reviewed by Lovelace's Square: Yes
%
% Detailed function description:
% This function takes an input matrix D and applies a set of kernels to it,
% resulting in a kernelized version of the matrix. The kernels include
% exponential-like patterns as well as impulses, and the function handles 
% both even and odd numbers of kernels. The kernels are normalized and then
% each signal is convolved with each kernel.
%
% Args:
%    D (matrix): The input data matrix to be kernelized, with dimensions
%                [num_samples, num_timepoints].
%                Each row is a signal (e.g., time series) from a sample.
%    num_kernels (int): The number of kernels to generate and apply.
%                       Must be a positive integer greater than or equal to 2.
%    kernel_width (int): The width of each kernel. Must be a positive integer 
%                        less than or equal to num_timepoints to produce 
%                        meaningful results.
%
% Returns:
%    D_kernelized (3D matrix): The kernelized data matrix with dimensions
%       [num_samples, num_kernels, num_timepoints - kernel_width + 1].
%       This is a three-way array where:
%         - dimension 1 corresponds to samples
%         - dimension 2 corresponds to different kernels
%         - dimension 3 corresponds to the time dimension after convolution
%
% Use example:
%    D = rand(10, 100); % Example data
%    num_kernels = 5;
%    kernel_width = 20;
%    D_kernelized = kernelize(D, num_kernels, kernel_width);
%
% Disclaimer:
% Authors and Lovelace's Square are not responsible for any issues, 
% inaccuracies, or data loss arising from the use of this function.


    if ~isnumeric(D) || ~ismatrix(D)
        error('kernelize:InvalidInput', ...
              'Input D must be a numeric 2D matrix.');
    end

    [num_samples, num_timepoints] = size(D);
    if num_samples < 1 || num_timepoints < 1
        error('kernelize:InvalidDimensions', ...
              'Input matrix D must have at least one sample and one timepoint.');
    end

    if ~isscalar(num_kernels) || num_kernels < 2 || floor(num_kernels) ~= num_kernels
        error('kernelize:InvalidNumKernels', ...
              'num_kernels must be an integer >= 2.');
    end

    if ~isscalar(kernel_width) || kernel_width < 1 || kernel_width > num_timepoints || floor(kernel_width) ~= kernel_width
        error('kernelize:InvalidKernelWidth', ...
              ['kernel_width must be a positive integer, and less than or equal to ', ...
               'the number of timepoints in D.']);
    end

    % If timepoints are too few to produce a valid convolution result
    output_length = num_timepoints - kernel_width + 1;
    if output_length < 1
        error('kernelize:NotEnoughTimepoints', ...
              'The specified kernel_width is too large for the given data length.');
    end

    kernels = zeros(num_kernels, kernel_width);

    % First kernel: [1 0 0 ... 0]
    kernels(1, :) = [1, zeros(1, kernel_width - 1)];

    % Last kernel: [0 0 ... 0 1]
    kernels(end, :) = [zeros(1, kernel_width - 1), 1];

    % Generate decaying exponentials for half the kernels
    half_kernels = floor((num_kernels - 2) / 2);
    for i = 2:(half_kernels + 1)
        tau = 1 / (i * 0.5);  % adjustable rate parameter
        decay = exp(-linspace(0, 1, kernel_width) / tau);
        kernels(i, :) = decay;
        kernels(num_kernels - i + 1, :) = fliplr(decay);
    end

    % If odd number of kernels, add a symmetric one in the middle
    if mod(num_kernels, 2) == 1
        mid = ceil(num_kernels / 2);
        kernels(mid, :) = exp(-abs(linspace(-1, 1, kernel_width)));
    end

    for i = 1:num_kernels
        max_val = max(kernels(i, :));
        if max_val == 0
            % If a kernel is all zeros (which is unlikely but possible), 
            % normalization would lead to division by zero.
            % We can either leave it as zeros or raise a warning.
            warning('kernelize:ZeroKernel', ...
                    'Kernel %d is all zeros. Leaving it as is.', i);
        else
            kernels(i, :) = kernels(i, :) / max_val;
        end
    end


    D_kernelized = zeros(num_samples, num_kernels, output_length);

    for i = 1:num_samples
        % Check if the data row has NaNs or Infs that might cause issues
        if any(isnan(D(i,:))) || any(isinf(D(i,:)))
            warning('kernelize:NonFiniteData', ...
                    'Sample %d contains NaN or Inf values. Results may be unreliable.', i);
        end

        for j = 1:num_kernels
            % Convolve with 'valid' option so that output length is (num_timepoints - kernel_width + 1)
            convolved = conv(D(i, :), kernels(j, :), 'valid');
            if length(convolved) ~= output_length
                % This should not happen if dimensions are correct, but a safety check is added
                error('kernelize:UnexpectedConvolutionOutput', ...
                      'Unexpected convolution output length. Check kernel and data dimensions.');
            end
            D_kernelized(i, j, :) = convolved;
        end
    end

end
