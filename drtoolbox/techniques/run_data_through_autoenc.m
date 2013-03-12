function [reconX, mappedX] = run_data_through_autoenc(network, X)
%RUN_DATA_THROUGH_AUTOENC Intermediate representation and reconstruction
%
%   [reconX, mappedX] = run_data_through_autoenc(network, X)
%
% Computes intermediate representation and reconstruction of the specified
% data from the specified autoencoder.
%
%
% (C) Laurens van der Maaten
% Delft University of Technology, 2011

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    % Initialize some variables
    n = size(X, 1);
    no_layers = length(network);
    middle_layer = ceil(no_layers / 2);

    % Run data through autoencoder
    activations = [X ones(n, 1)];
    for i=1:no_layers
        if i ~= middle_layer && i ~= no_layers
            activations = [1 ./ (1 + exp(-(activations * [network{i}.W; network{i}.bias_upW]))) ones(n, 1)];
        else
            activations = [activations * [network{i}.W; network{i}.bias_upW] ones(n, 1)];
            if i == middle_layer
                mappedX = activations(:,1:end-1);
            end
        end
    end
    reconX = activations(:,1:end-1);
