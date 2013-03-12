function [mappedX, mapping] = fa(X, no_dims)
%FA Perform factor analysis (FA) on dataset X
%
%   mappedX = fa(X, no_dims)
%
% Perform factor analysis on dataset X in order to reduce its
% dimensionality to no_dims dimensions. The reduced data representation is
% returned in mappedX.

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    if ~exist('no_dims', 'var')
        no_dims = 2;
    end
    
    % Make sure data is zero-mean
    mapping.mean = mean(X, 1);
    X = X - repmat(mapping.mean, [size(X, 1) 1]);

    % Initialize some variables
    X = X';
    [D, n] = size(X);
    epsilon = 1e-5;
    iter = 0;
    max_iter = 200;

    % Initialize FA model
    Sig = eye(D);               % initial variances
    A = rand(D, no_dims);       % initial linear mapping

    % Main loop
    while iter < max_iter
       
        % Update number of iterations
        iter = iter + 1;
        if rem(iter, 5) == 0
            fprintf('.');
        end

        % Perform E-step
        invC = inv(A * A' + Sig);                               % compute inverse of covariance matrix
        M = A' * invC * X;
        SC = n * (eye(no_dims) - A' * invC * A) + M * M';

        % Perform M-step
        A = (X * M') * inv(SC);
        Sig = (diag(diag(X * X' - A * M * X')) / n) + epsilon;

        % Compute log-likelihood of FA model
        newll = 0.5 * (log(det(invC)) - sum(sum((invC * X) .* X)) / n);

        % Check for convergence
        if iter ~=1 && abs(newll - ll) < epsilon
            break;
        end
        ll = newll;
    end
    
    % Apply linear mapping
    mapping.M = A;
    mappedX = X' * mapping.M;
    disp(' ');
