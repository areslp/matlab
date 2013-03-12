function autoenc = roll_out(network)
%ROLL_OUT Rolls out an autoencoder
% 
%   network = roll_out(network)
%
% Rolls out a stack of RBMs as an autoencoder that can be finetuned using,
% e.g., backpropagation.
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    % Initialize some variables
    no_layers = length(network);
    autoenc = cell(1, 2 * no_layers);
    
    % Roll out network
    for i=1:no_layers
        autoenc{i} = network{i};
        autoenc{2 * no_layers - i + 1}.W          = network{i}.W';
        autoenc{2 * no_layers - i + 1}.bias_upW   = network{i}.bias_downW;
        autoenc{2 * no_layers - i + 1}.bias_downW = network{i}.bias_upW;
    end
