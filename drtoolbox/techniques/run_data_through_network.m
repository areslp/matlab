function mappedX = run_data_through_network(network, X)
%RUN_DATA_THROUGH_NETWORK Run data through the network
%
%   mappedX = run_data_through_network(network, X)
%
% Runs the dataset X through the parametric t-SNE embedding defined in
% network. The result is returned in mappedX.
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    % Run the data through the network
    n = size(X, 1);
    mappedX = [X ones(n, 1)];
    for i=1:length(network) - 1
        mappedX = [1 ./ (1 + exp(-(mappedX * [network{i}.W; network{i}.bias_upW]))) ones(n, 1)];
    end
    mappedX = mappedX * [network{end}.W; network{end}.bias_upW];
    