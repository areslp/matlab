function [F, dF] = nca_lin_grad(x, X, labels, no_dims, lambda)
%NCA_LIN_GRAD Computes NCA gradient on the specified dataset
%
%   [C, dC] = nca_lin_grad(x, X, labels, no_dims, lambda)
%
% Computes the linear Neighborhood Components Analysis (NCA) gradient 
% on the dataset specified through X and labels to reduce the data 
% dimensionality to no_dims dimensions. The current solutions is specified 
% in x. The function returns the current value of the cost function in C 
% and the gradient in dC.
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

    % Initialize some variables
    n = size(X, 1);
    A = reshape(x, [numel(x) / no_dims no_dims]);
    F = 0;
    dF = zeros(size(A));

    % Transform the data
    Y = X * A;
    
    % Compute conditional probabilities for current solution
    sumY = sum(Y .^ 2, 2);
    P = exp(bsxfun(@minus, bsxfun(@minus, 2 * (Y * Y'), sumY'), sumY));
    P(1:n+1:end) = 0;
    P = bsxfun(@rdivide, P, sum(P, 2));
    P = max(P, eps);

    % Compute value of cost function and gradient
    for i=1:n
        
        % Sum cost function
        inds = (labels == labels(i));
        Pi = sum(P(i, inds));
        F = F + Pi;

        % Sum gradient
        if nargout > 1
            xikA = bsxfun(@minus, Y(i,:), Y);
            xik  = bsxfun(@minus, X(i,:), X);
            dF = dF + Pi * (bsxfun(@times, xik,         P(i,:)')'     * xikA) - ...
                            bsxfun(@times, xik(inds,:), P(i, inds)')' * xikA(inds,:);
        end
    end
    
    % Include regularization term
    F = F - lambda .* sum(A(:) .^ 2) ./ numel(A);
    dF = 2 * dF - 2 * lambda .* A ./ numel(A);

    % Prepare to pass to minimize
    F = -F;
    dF = -dF(:);
