function [X, W, mu_X] = prewhiten(X)
%PREWHITEN Performs prewhitening of a dataset X
%
%   [X, W, mu_X] = prewhiten(X)
%
% Performs prewhitening of the dataset X. Prewhitening concentrates the main
% variance in the data in a relatively small number of dimensions, and 
% removes all first-order structure from the data. In other words, after
% the prewhitening, the covariance matrix of the data is the identity
% matrix. The function returns the subtracted data mean in mu_X, and the
% applied linear mapping in W.
% 
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    welcome;

    % Compute and apply the ZCA mapping
    mu_X = mean(X, 1);
    X = bsxfun(@minus, X, mu_X);
    mappedX = X / sqrtm(cov(X));
    if nargout > 1
        W = X \ mappedX;
    end   
    