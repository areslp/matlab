function reconX = recon_data_from_autoenc(network, mappedX)
%RUN_DATA_THROUGH_AUTOENC Summary of this function goes here
%   Detailed explanation goes here

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    % Initialize some variables
    n = size(mappedX, 1);
    no_layers = length(network);
    
    % Run data through autoencoder
    activations = [mappedX ones(n, 1)];
    middle_layer = ceil(no_layers / 2);
    for i=middle_layer + 1:no_layers
        if i ~= no_layers
            activations = [1 ./ (1 + exp(-(activations * [network{i}.W; network{i}.bias_upW]))) ones(n, 1)];
        else
            activations = [activations * [network{i}.W; network{i}.bias_upW] ones(n, 1)];
        end
    end
    reconX = activations(:,1:end-1);
