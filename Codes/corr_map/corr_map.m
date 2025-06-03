function corr_map(C)
% Plot the correlation map
%
% Author: Ruggero Guerrini
% Mail: gr6898@gmail.com
% Date Created: 2025-04-08
% License: MIT
%
% Description:
% Plot the correlation map printing the number in each pixel of the map
%
%
% Arguments:
%    C: Input data, spectrum. n samples x p spectral length
%
% Returns:
%    figure with the correlation map
%
% Disclaimer:
% Authors and Lovelace's Square are not responsible for any issues, inaccuracies, or data loss arising
% from the use of this function.
%

figure;
cor_ = corr(C');
imagesc(cor_);
colormap('jet');
colorbar;
axis square;
title('Correlation Map');
md = median(cor_(:));


% Uncomment these two lines if you want the range of colorbar between -1
% and 1
% clim([-1 1])
% md=0;

for i = 1:size(C,1)
    for j = 1:size(C,1)
        value = cor_(i,j);
        if value<md
                text(j, i, sprintf('%.2f', value), ...
            'HorizontalAlignment', 'center', ...
            'Color',[1 1 1], ...
            'FontWeight', 'bold');
        else
                    text(j, i, sprintf('%.2f', value), ...
            'HorizontalAlignment', 'center', ...
            'Color',[0 0 0], ...
            'FontWeight', 'bold');
        end
    end
end
