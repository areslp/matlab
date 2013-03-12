function [C, dC] = backprop_gradient(v, network, X, targets, lambda)
%BACKPROP Compute the cost gradient for CG optimization of a neural network
%
%   [C, dC] = backprop_gradient(v, network, X, targets, lambda)
%
% Compute the value of the cost function, as well as the corresponding 
% gradient for conjugate-gradient optimization of a neural network.
%
%

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

    % Deconvert the weights and store them in the network
    ind = 1;
    for i=1:no_layers
        network{i}.W        = reshape(v(ind:ind - 1 + numel(network{i}.W)),        size(network{i}.W));         ind = ind + numel(network{i}.W);
        network{i}.bias_upW = reshape(v(ind:ind - 1 + numel(network{i}.bias_upW)), size(network{i}.bias_upW));  ind = ind + numel(network{i}.bias_upW);
    end
    
    % Run the data through the network
    activations = cell(1, no_layers + 1);
    activations{1} = [X ones(n, 1)];
    for i=1:no_layers
        if i ~= middle_layer && i ~= no_layers
            activations{i + 1} = [1 ./ (1 + exp(-(activations{i} * [network{i}.W; network{i}.bias_upW]))) ones(n, 1)];
        else
            activations{i + 1} = [activations{i} * [network{i}.W; network{i}.bias_upW] ones(n, 1)];
        end
    end  

    % Compute value of cost function (= MSE)
    C = (1 / (2 * n)) .* sum(sum((activations{end}(:,1:end - 1) - targets) .^ 2)) + lambda .* sum(v .^ 2);
    
    % Only compute gradient if requested
    if nargout > 1
    
        % Compute gradients 
        dW = cell(1, no_layers);
        db = cell(1, no_layers);
        Ix = (activations{end}(:,1:end - 1) - targets) ./ n;
        for i=no_layers:-1:1

            % Compute update
            delta = activations{i}' * Ix;
            dW{i} = delta(1:end - 1,:);
            db{i} = delta(end,:);

            % Backpropagate error
            if i > 1
                if i ~= middle_layer + 1
                    Ix = (Ix * [network{i}.W; network{i}.bias_upW]') .* activations{i} .* (1 - activations{i});
                else
                    Ix = Ix * [network{i}.W; network{i}.bias_upW]';
                end
                Ix = Ix(:,1:end - 1);
            end
        end

        % Convert gradient information
        dC = zeros(numel(v), 1);
        ind = 1;
        for i=1:no_layers
            dC(ind:ind - 1 + numel(dW{i})) = dW{i}(:); ind = ind + numel(dW{i});
            dC(ind:ind - 1 + numel(db{i})) = db{i}(:); ind = ind + numel(db{i});
        end
        dC = dC + 2 .* lambda .* v;
    end
