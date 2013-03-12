function [C, dA] = mcml_grad(x, X, P)
%MCML_GRAD Computes MCML cost function and gradient 
%
%   [C, dA] = mcml_grad(x, X, P)
%
% Computes MCML cost function and gradient.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    % Decode solution
    [n, d] = size(X);
    A = reshape(x, [d d]);
    
    % Compute conditional probabilities for current solution
    D = zeros(n, n);
    for i=1:n
        diffX = bsxfun(@minus, X(i,:), X(i + 1:end,:));
        dist = sum((diffX * A) .* diffX, 2);
        D(i + 1:end, i) = dist;
        D(i, i + 1:end) = dist';
    end
    Q = exp(-D);
    Q(1:n+1:end) = 0;
    Q = bsxfun(@rdivide, Q, sum(Q, 2));
    Q = max(Q, realmin);
    
    % Compute cost function
    C = sum(P(:) .* log(P(:) ./ Q(:)));

    % Compute gradient with respect to A
    if nargin > 1
        dA = zeros(d, d);
        PQ = P - Q;
        for i=1:n
            diffX = bsxfun(@minus, X(i,:), X);
            dA = dA + bsxfun(@times, PQ(i,:), diffX') * diffX;
        end
    end
    dA = dA(:);
    