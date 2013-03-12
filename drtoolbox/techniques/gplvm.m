function Y = gplvm(X, no_dims, sigma)
%GPLVM Gaussian Process Latent Variable Model
%
%   Y = gplvm(X, no_dims, sigma)
%
% Simple implementation of the Gaussian Process Latent Variable Model using 
% a Gaussian kernel with bandwidth alpha. The function reduces the
% dimenisonality of the dataset X to no_dims dimensions. The resulting
% low-dimensional dataset is returned in Y.
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    if ~exist('sigma', 'var') || isempty(sigma)
        sigma = 1;
    end

    % Initialize some variables
    n = size(X, 1);
    
    % Initialize solution using PCA
    disp('Preprocessing data using PCA...');
    X = bsxfun(@minus, X, mean(X, 1));
    [M, lambda] = eig(X' * X);
    [lambda, ind] = sort(diag(lambda), 'descend');
    M = M(:,ind(1:no_dims));
    Y = X * M;
    clear M lambda ind
    
    % Check gradient
    disp('Learn GPLVM model...');
    Y = minimize(Y(:), 'gplvm_grad', -500, X, sigma);
    Y = reshape(Y, [n no_dims]);
    