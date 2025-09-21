function [Matrix] = unfoldImage(Cube)
% UNFOLDIMAGE. Flatten a 3D image cube into a 2D matrix – because chemometrics 
% models only date flat data, and Lovelace's Square is here to make your life easier.
% We want to protect your fingers from having to type the dimensions of your images 
% over and over again when performing the unfolding.
%
% Author: Adrián Gómez-Sánchez
% Date Created: 2024-12-16
% License: MIT. Do we need a License for this code? Haha
% Reviewed by Lovelace's Square: Yes
% Version: 1.0
%
% The UNFOLDIMAGE function takes a 3D image cube and flattens it into a 2D matrix 
% by merging the first two dimensions. This is done because almost all chemometrics 
% models work with matrices instead of cubes – extra dimensions are just too much hassle!
%
% INPUT:
%   Cube (array): A 3D numeric array of dimensions [x, y, z] representing an image cube.
%
% OUTPUT:
%   Matrix (array): A 2D numeric matrix of size (x*y) x z, with the spatial dimensions flattened.
%
% EXAMPLE:
%   % Suppose 'imgCube' is a 3D image array of dimensions 100x100x3:
%   unfolded = unfoldImage(imgCube);
%   % 'unfolded' will be a 10000x3 matrix – nice and flat, just like all our chemometric models like it.
%
% DISCLAIMER:
%   Authors and Lovelace's Square are not responsible for any existential crises 
%   that might arise from the difficulty of this code.

    % Retrieve the dimensions of the input cube.
    [x, y, z] = size(Cube);
    
    % Reshape the cube into a 2D matrix by merging the first two dimensions.
    Matrix = reshape(Cube, x*y, z);
end

