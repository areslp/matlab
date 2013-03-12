function [model, mappedX] = train_deep_autoenc(X, layers, lambda)
%TRAIN_DEEP_AUTOENC Trains a deep feedforward autoencoder on X
%
%   [model, mappedX] = train_deep_autoenc(X, layers, lambda)
%
% The function trains a deep feedforward autoencoder on
% the specified X and labels. The variable X should be a Nx1 cell
% array containing DxT matrices (where T can vary). The variable labels
% should be a Nx1 cell array containing KxT matrices. The variable lambda
% specifies the L2 regularization parameter (default = 0). The parameter
% base_eta represents the base learning rate (default = 1).
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology



    if ~exist('lambda', 'var') || isempty(lambda)
        lambda = 0;
    end
    
    % Pretrain model using stacked denoising auto-encoders
    no_layers = length(layers);
    model = cell(2 * no_layers, 1);
    mappedX = X;
    for i=1:no_layers
        noise = 0.1;
        max_iter = 30;
        [network, mappedX] = train_autoencoder(mappedX, layers(i), noise, max_iter);
        model{i}.W        = network{1}.W;
        model{i}.bias_upW = network{1}.bias_upW;
    end
    for i=1:no_layers
        model{no_layers + i}.W        = model{no_layers - i + 1}.W';
        if i ~= no_layers
            model{no_layers + i}.bias_upW = model{no_layers - i}.bias_upW;
        else
            model{no_layers + i}.bias_upW = zeros(1, size(X, 2));
        end
    end
    clear network mappedX
    
    % Compute mean squared error of initial model predictions
    reconX = run_data_through_autoenc(model, X);
    disp(['MSE of initial model: ' num2str(mean((reconX(:) - X(:)) .^ 2))]);
    
    % Finetune model using gradient descent
    noise = 0.1;
    max_iter = 30;
    model = backprop(model, X, X, max_iter, noise, lambda);
    
    % Compute mean squared error of final model predictions
    [reconX, mappedX] = run_data_through_autoenc(model, X);
    disp(['MSE of final model: ' num2str(size(X, 2) .* mean((reconX(:) - X(:)) .^ 2))]);
end
