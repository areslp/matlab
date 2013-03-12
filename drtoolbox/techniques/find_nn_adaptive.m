function [D, max_k_val, no_dims] = find_nn_adaptive(X)
%FIND_NN Constructs nearest neighbor graph using adaptive nbhd selection
%
%	[D, max_k_val, no_dims] = find_nn_adaptive(X)
%
% Constructs nearest neighbor graph on the data in dataset X using adaptive
% neighborhood selection. In X, rows correspond to the observations and 
% columns to the dimensions. 
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    % First, estimate intrinsic dimensionality using the MLE estimator
    no_dims = round(intrinsic_dim(X, 'MLE'));

    % Initialize some variables
    n = size(X, 1);
    min_k = no_dims + 1;
    max_k = n;
    max_k_val = 0;

    % Perform some precomputations for speed
    X = X';
    D = sparse([], [], [], n, n, n * 12);           % Pre-allocate for speed!
    XX = sum(X .* X);
    onez = ones(1,n);
    
    % For all datapoints
    for i=1:n
        % current data point
        p = X(:,i);

        % Compute Euclidean distance to all other datapoints
        d2 = sum(p.^2).*onez + XX - 2*p'*X;

        % Sort distances
        [d2, ind] = sort(d2);
        
        % Estimate local tangent space by updating the number of neighbors k
        stop = 0;
        k = min_k;
        while ~stop && k + 1 < max_k
            
            % Update k
            k = k + 1;
            
            % Estimate local tangent space (for current value of k)
            tmpX = X(:,ind(2:k + 1)) - repmat(p, 1, k);
            lambda = svd(tmpX);
            [lambda, ind2] = sort(lambda, 'descend');
            if length(lambda) < no_dims
                break;
            end

            % Estimate T_{1}
            T = (1 / k) ^ (1 / no_dims) * sqrt(d2(k + 1));

            % Check whether stop condition is violated
            if lambda(no_dims) >= T
                stop = true;
            end
        end

        % Compute tangent space at k - 1 since k failed stop condition
        [U, lambda, M] = svd(tmpX(:,1:k-1));
        [lambda, ind2] = sort(diag(lambda), 'descend');
        U = U(:,ind2(1:no_dims))';

        % Select neighbors that correspond to the local tangent space
        stop = 0;
        while ~stop && k + 1 < max_k
            % Update k
            k = k + 1;

            % Projection of (x_{k} - x_{i}) onto tangent space
            onto = sum( (U * (X(:,ind(k + 1)) - p)).^2 );

            % Check whether stop condition is violated
            if d2(k + 1) - onto > T^2
                k = k - 1;
                stop = true;
            end
        end
        
        % Update neighborhood graph
        d2(1) = 1e-7;
        D(i, ind(1:k + 1)) = sqrt(d2(1:k + 1));
        if max_k_val < k
            max_k_val = k;
        end
    end
