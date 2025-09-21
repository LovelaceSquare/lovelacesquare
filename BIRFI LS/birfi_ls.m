function [irf] = birfi_ls(decay, irf_size, lambda)
% BIRFI_LS. Perform IRF estimation using Tikhonov regularization.
%
% REFERENCES:
%   Gómez-Sánchez, Adrián, et al. "Blind instrument response function 
%   identification from fluorescence decays." Biophysical Reports 4.2 
%   (2024).
%
% Author: Adrián Gómez-Sánchez
% Date Created: 2024-12-16
% License: MIT
% Reviewed by Lovelace's Square: Yes
% Version: 1.0
%
% The BIRFI_LS function estimates the instrument response function (IRF) 
% from a given decay signal using a Hankel matrix approach with Tikhonov 
% regularization. This method reconstructs the IRF by solving an inverse 
% problem with a smoothing constraint controlled by the regularization 
% parameter lambda.
%
% Specifically, the IRF is obtained by solving the linear system:
%     (HankelM' * HankelM + lambda * I) * irf' = HankelM' * decay'
% where:
%   - HankelM is a structured matrix constructed from the decay data,
%   - lambda is the regularization parameter controlling smoothness,
%   - I is the identity matrix.
%
% The function returns:
%   - irf: The estimated instrument response function.
%
% INPUTS:
%   decay (array)   : A 1D numeric vector (1xD) representing the full decay.
%   irf_size (int)  : A scalar defining the size of the IRF (cutting point).
%   lambda (double) : A positive regularization parameter controlling the 
%                     smoothness of the estimated IRF. If not provided, a 
%                     default value of 100000 is used.
%
% OUTPUT:
%   irf (array): A 1D numeric vector (1xirf_size) representing the estimated IRF.
%
% EXAMPLE:
%   % Suppose 'decay_signal' is your measured decay data:
%   irf_length = 50;  % Define IRF size
%   reg_param  = 100000; % Regularization strength
%   estimated_irf = birfi_ls(decay_signal, irf_length, reg_param);
%   % 'estimated_irf' now contains the reconstructed IRF.
%
% DISCLAIMER:
%   Authors and Lovelace's Square are not responsible for any issues, 
%   inaccuracies, or data loss arising from the use of this function.

    % Validate input arguments
    if nargin < 2
        error('Usage: birfi_ls(decay, irf_size, lambda). At least decay and irf_size are required.');
    end
    if ~isvector(decay) || isempty(decay)
        error('Input "decay" must be a nonempty 1D numeric vector.');
    end
    if ~isscalar(irf_size) || irf_size <= 0
        error('Input "irf_size" must be a positive integer scalar.');
    end
    if nargin < 3
        lambda = 100000;  % Default regularization parameter
    elseif lambda <= 0
        error('The regularization parameter "lambda" must be positive.');
    end

    % Determine decay length and define cutoff
    decay = decay(:)';  % Ensure row vector
    decaysize = length(decay);
    if irf_size > decaysize
        error('irf_size cannot be greater than the length of the decay data.');
    end
    decay_cut = decay(irf_size:end);

    % Construct an auxiliary vector with zero-padding
    f2 = [zeros(1, irf_size-1), decay_cut, zeros(1, irf_size-1)];

    % Preallocate Hankel matrix
    HankelM = zeros(decaysize, irf_size);

    % Build the Hankel matrix: each row is a reversed segment of f2
    for i = 1:decaysize
        HankelM(i, :) = fliplr(f2(i : i + irf_size - 1));
    end

    % Solve the regularized system:
    %   (HankelM' * HankelM + lambda * I) * irf' = HankelM' * decay'
    irf = (HankelM' * HankelM + lambda * eye(irf_size)) \ (HankelM' * decay');

    % Convert to row vector
    irf = irf';
end

