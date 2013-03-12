function [C, dC] = gplvm_grad(x, X, sigma)
%GPLVM_GRAD Gradient of the Gaussian Process Latent Variable model
%
%   [C, dC] = gplvm_grad(x, no_dims, sigma)
%
% Computes the gradient of the Gaussian Process Latent Variable model.

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    % Decode solution
    [n, d] = size(X);
    no_dims = numel(x) ./ n;
    Y = reshape(x, [n no_dims]);

    % Compute kernel matrix (in latent space)
    sum_Y = sum(Y .^ 2, 2);
    K = exp(bsxfun(@minus, bsxfun(@minus, 2 * (Y * Y'), sum_Y'), sum_Y) / (2 * sigma ^ 2));
    
    % Compute gradient with respect to kernel
%     invK = inv(K);
%     dLdK = invK * X * X' * invK - d * invK;
    tmp = (K \ X) * X';
    dLdK = (tmp - d * eye(n)) / K;

    % Compute gradient with respect to coordinates
    dC = zeros(n, no_dims);
    dLdK = K .* dLdK;
    for i=1:n
        dC(i,:) = sum(bsxfun(@times, dLdK(:,i), bsxfun(@minus, Y(i,:), Y)), 1);
    end
    dC = (-1 / sigma ^ 2) * dC(:);
    
    % Compute log-likelihood
    C = -((d * n) / 2) * log(2 * pi) - (d / 2) * log(det(K) + realmin) - .5 * trace(tmp);
