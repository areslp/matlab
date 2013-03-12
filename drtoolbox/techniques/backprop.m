function network = backprop(network, X, T, max_iter, noise, lambda)
%BACKPROP Trains a network on a dataset using backpropagation
%
%   network = backprop(network, X, T, max_iter, noise, lambda)
%
% The function trains the specified network using backpropagation on
% dataset X with targets T for max_iter iterations. The dataset X is an NxD
% matrix, whereas the targets matrix T has size NxM. The function returns 
% the trained network in network.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology



    if ~exist('max_iter', 'var') || isempty(max_iter)
        max_iter = 10;
    end
    if ~exist('noise', 'var') || isempty(noise)
        noise = 0;
    end
    if ~exist('lambda', 'var') || isempty(lambda)
        lambda = 0;
    end

    % Initialize some variables
    n = size(X, 1);
    no_layers = length(network);
    batch_size = max(round(n / 100), 100);
    
    % Perform the backpropagation
    for iter=1:max_iter
        disp([' - iteration ' num2str(iter) ' of ' num2str(max_iter) '...']);
        
        % Loop over all batches
        index = randperm(n);
        for batch=1:batch_size:n
            
            % Select current batch
            tmpX = X(index(batch:min([batch + batch_size - 1 n])),:);
            tmpT = T(index(batch:min([batch + batch_size - 1 n])),:);
            
            % Randomly black out some of the input data
            if noise > 0
                tmpX(rand(size(tmpX)) < noise) = 0;
            end   

            % Convert the weights and store them in the network
            v = [];
            for i=1:length(network)
                v = [v; network{i}.W(:); network{i}.bias_upW(:)];
            end
            
            % Conjugate gradient minimization using 3 linesearches
%             checkgrad('backprop_gradient', v, 1e-5, network, tmpX, tmpT)
            [v, fX] = minimize(v, 'backprop_gradient', 3, network, tmpX, tmpT, lambda);
            
            % Deconvert the weights and store them in the network
            ind = 1;
            for i=1:no_layers
                network{i}.W        = reshape(v(ind:ind - 1 + numel(network{i}.W)),        size(network{i}.W));         ind = ind + numel(network{i}.W);
                network{i}.bias_upW = reshape(v(ind:ind - 1 + numel(network{i}.bias_upW)), size(network{i}.bias_upW));  ind = ind + numel(network{i}.bias_upW);
            end
            
            % Stop upon convergence
            if isempty(fX)
                reconX = run_data_through_autoenc(network, X);
                C = sum((T(:) - reconX(:)) .^ 2) ./ n;
                disp([' - final noisy MSE: ' num2str(C)]);
                return
            end
        end
        
        % Estimate the current error
        reconX = run_data_through_autoenc(network, X);
        C = sum((T(:) - reconX(:)) .^ 2) ./ n;
        disp([' - current noisy MSE: ' num2str(C)]);
    end
    