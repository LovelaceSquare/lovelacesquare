function [irf] = birfi_ls_smoothing(decay, irf_size, penalty)
% BIRFI_LS_SMOOTHING. Perform IRF estimation with smoothing using Tikhonov regularization.
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
% The BIRFI_LS_SMOOTHING function estimates the instrument response function (IRF)
% from a given decay signal by constructing a Hankel matrix and solving the
% regularized equations. A smoothing penalty is enforced via a second-
% order finite difference operator.
%
% Specifically, the IRF is obtained by solving the linear system:
%     (HankelM' * HankelM + penalty * (D' * D)) * irf' = HankelM' * decay'
% where:
%   - HankelM is a structured matrix constructed from the decay data,
%   - penalty is the regularization parameter controlling smoothness,
%   - D is the second-order finite difference operator.
%
% The function returns:
%   - irf: The estimated instrument response function.
%
% INPUTS:
%   decay (array)    : A 1D numeric vector (1xD) representing the full decay.
%   irf_size (int)   : A scalar defining the size of the IRF.
%   penalty (double) : A positive smoothing regularization parameter that 
%                      controls the smoothness of the estimated IRF. If not provided,
%                      a default value of 1e8 is used.
%
% OUTPUT:
%   irf (array): A 1D numeric vector (1xirf_size) representing the estimated IRF.
%
% EXAMPLE:
%   % Suppose 'decay_signal' is your measured decay data:
%   irf_length = 100;      % Define IRF size
%   penaltyParam = 1e8;   % Smoothing regularization parameter
%   estimated_irf = birfi_ls_smoothing(decay_signal, irf_length, penaltyParam);
%   % 'estimated_irf' now contains the reconstructed IRF.
%
% DISCLAIMER:
%   Authors and Lovelace's Square are not responsible for any issues, 
%   inaccuracies, or data loss arising from the use of this function.

    % Set default penalty if not provided
    if nargin < 3
        penalty = 1e8;  % default smoothing regularization parameter
    end

    % Determine the length of the decay and define the cutoff.
    decaysize = length(decay);
    decay_cut = decay(irf_size:end);

    % Construct an auxiliary vector with zero-padding.
    f2 = [zeros(1, irf_size-1), decay_cut, zeros(1, irf_size-1)];

    % Build the Hankel matrix: each row is a reversed segment of f2.
    HankelM = zeros(decaysize, irf_size);
    for i = 1:decaysize
        HankelM(i, :) = fliplr(f2(i : i + irf_size - 1));
    end

    % Create the second difference operator D for the smoothing penalty.
    % D computes second differences so that D*x approximates the second derivative of x.
    D = diff(eye(irf_size), 2);

    % Solve the regularized normal equations:
    % (HankelM'*HankelM + penalty*(D'*D)) * irf' = HankelM'*decay'
    irf = (HankelM' * HankelM + penalty * (D' * D)) \ (HankelM' * decay(:));
    irf = irf';  % Convert to a row vector.
end


