function [network, mappedX, reconX] = train_autoencoder(X, layers, noise, max_iter)
%TRAIN_AUTOENCODER Trains an simple autoencoder
%
%   [network, mappedX, reconX] = train_encoder(X, layers, noise, max_iter)
%
% Trains up an autoencoder with the structure that is specified in layers. 
% The low-dimensional data is returned in mappedX, and the network in
% network. The variable noise indicates how much of the input data should be
% blacked out (default = 0).
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    if nargin < 2
        error('Not enough inputs.');
    end
    if isempty(layers)
        error('There should be at least one hidden layer.');
    end
    if ~exist('noise', 'var') || isempty(noise)
        noise = 0;
    end
    if ~exist('max_iter', 'var') || isempty(max_iter)
        max_iter = 50;
    end
    
    % Initialize the network
    D = size(X, 2);
    no_layers = length(layers) + 1;
    network = cell(no_layers, 1);
    network{1}.W = randn(D, layers(1)) * .0001;
    network{1}.bias_upW = zeros(1, layers(1));
    for i=2:no_layers - 1
        network{i}.W = randn(layers(i - 1), layers(i)) * .0001;
        network{i}.bias_upW = zeros(1, layers(i));
    end
    network{no_layers}.W = randn(layers(end), D) * .0001;
    network{no_layers}.bias_upW = zeros(1, D);
    reconX = run_data_through_autoenc(network, X);
    disp(['Initial MSE of reconstructions: ' num2str(mean((X(:) - reconX(:)) .^ 2))]);    
    
    % Perform backpropagation to minimize reconstruction error
    network = backprop(network, X, X, max_iter, noise);
    
    % Get representation from hidden layer
    [reconX, mappedX] = run_data_through_autoenc(network, X);
    disp(['Final MSE of reconstructions: ' num2str(mean((X(:) - reconX(:)) .^ 2))]);
    