function [mappedX, mapping] = nca(X, labels, no_dims, lambda)
%NCA Performs NCA on the specified dataset
%
%   [mappedX, mapping] = nca(X, labels, no_dims, lambda)
%
% Performs linear Neighborhood Components Analysis (NCA) on the 
% dataset specified through X and labels to reduce the data dimensionality 
% to no_dims dimensions. If valid_X and valid_labels are specified, the
% function uses early stopping based on NN errors.
% The function returns a embedded data in mappedX, as well as th mapping.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology



    if ~exist('no_dims', 'var') || isempty(no_dims)
        no_dims = size(X, 2);
    end
    if ~exist('lambda', 'var') || isempty(lambda)
        lambda = 0;
    end
    
    % Make sure data is zero mean
    mapping.mean = mean(X, 1);
	X = bsxfun(@minus, X, mapping.mean);

    % Initialize some variables
    max_iter = 200;
    [n, d] = size(X);
    batch_size = min(5000, n);
    no_batches = ceil(n / batch_size);
    max_iter = ceil(max_iter / no_batches);
    [lablist, foo, labels] = unique(labels);
    A = randn(d, no_dims) * .01;
    converged = false;
    iter = 0;
    
    % Main iteration loop
    while iter < max_iter && ~converged
        
        % Loop over batches
        iter = iter + 1;
        disp(['Iteration ' num2str(iter) ' of ' num2str(max_iter) '...']);
        ind = randperm(n);
        for batch=1:batch_size:n
    
            % Run NCA minimization using three linesearches
            cur_X    = double(X(ind(batch:min([batch + batch_size - 1 n])),:));
            cur_labels = labels(ind(batch:min([batch + batch_size - 1 n])));
            [A, f] = minimize(A(:), 'nca_lin_grad', 5, cur_X, cur_labels, no_dims, lambda);
            A = reshape(A, [d no_dims]);
            if isempty(f) || f(end) - f(1) > -1e-4
                disp('Converged!');
                converged = true;
            end
        end
    end
    
    % Compute embedding
    mapping.M = A;
    mappedX = X * mapping.M;
    