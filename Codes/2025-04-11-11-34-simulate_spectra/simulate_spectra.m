function spectra = simulate_spectra(num_spectra,num_variables,n_peaks,width,type)

% Author: Ruggero Guerrini
% Mail: gr6898@gmail.com
% Data Created: 2025-04-11
% Licence: MIT

% Version: 1.0

% Description:
% Simulate spectra with n gaussian without noise.
% Improvement will follow.
% There are two options:
% 'same'
% For each spectrum the position of the guassian is random
% 'different'
% The position of the gaussian is the same. Widht and amplitude change

% Input description:
% 'num_spectra' : number of spectra desired
% 'num_variables' : number of spectral points
% 'n_peaks': number of peaks for each spectrum desired
% 'width': width of the peaks. Scale factor of rand to estimate the width
% for each gaussian
% 'type': same or different, as described above.

% Disclaimer:
% Author and Lovelace's Square are not responsible for any issues,
% inaccuracies, or data loss arising from the use of this function.



switch type
    case 'same'
        spectra = zeros(num_spectra,num_variables);
        x = linspace(1, num_variables+1, num_variables);
        mu = randperm(num_variables-round(0.1*num_variables),n_peaks)+round(0.05*num_variables);
        sigma = width*rand(1,n_peaks);

        for j=1:num_spectra
            spectrum = zeros(size(x));
            for i = 1:length(mu)
                mu_ = mu(i);
                sigma_ = sigma(i);
                amplitude = rand*10;
                spectrum = spectrum +amplitude* exp(-((x - mu_).^2) / (2 * sigma_^2));
            end
            spectrum = spectrum / max(spectrum);
            spectra(j,:) = spectrum;
        end
    case 'different'
        spectra = zeros(num_spectra,num_variables);
        x = linspace(1, num_variables+1, num_variables);
        for j=1:num_spectra
            mu = randperm(num_variables,n_peaks);
            sigma = width*rand(1,n_peaks);
            spectrum = zeros(size(x));
            for i = 1:length(mu)
                mu_ = mu(i);
                sigma_ = sigma(i);
                amplitude = rand*10;
                spectrum = spectrum +amplitude* exp(-((x - mu_).^2) / (2 * sigma_^2));
            end
            spectrum = spectrum / max(spectrum);
            spectra(j,:) = spectrum;
        end
end
figure;
plot(spectra', 'LineWidth', 2);
xlabel('Variables');
ylabel('Intensity');
title('Spectra simualted');
grid on;