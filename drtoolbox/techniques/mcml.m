function mapping = mcml(X, labels, no_dims)
%MCML Performs MCML on the specified dataset
%
%   mapping = mcml(X, labels, no_dims)
%
% Performs linear Maximally Collapsing Metric Learning (MCML) on the 
% dataset specified through X and labels. The function returns a 
% Mahalanobis metrix in M. Distances through the metric are thus given by
% D(x, y) = (x - y) * A * (x - y)'.
% It the does an SVD on A gives you a (low-rank) projection of the data
% into a space where the squared Euclidean distance corresponds to the
% learned Mahalanobis metric in the original space.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


	% Make sure data is zero mean
    mapping.mean = mean(X, 1);
	X = bsxfun(@minus, X, mapping.mean);

    % Initialize some variables
    [n, d] = size(X);
    max_iter = 100;
    
    % Construct matrix with label information
    P = zeros(n, n);
    [foo, bar, labels] = unique(labels);
    for i=1:max(labels)
        ind = find(labels == i);
        P(ind, ind) = 1;
    end
    P(1:n+1:end) = 0;
    P = bsxfun(@rdivide, P, sum(P, 2));
    P = max(P, realmin);
    
    % Initialize solution
    A = randn(d, d) * .0001;
    A = (A + A') + eye(d);
    iter = 0;
    converged = false;
    
    % Run iterations
    while iter < max_iter && ~converged
        
        % Perform gradient step
        iter = iter + 1;
        disp(['Iteration ' num2str(iter) '...']);
        [A, f] = minimize(A(:), 'mcml_grad', 5, X, P);
        A = reshape(A, [d d]);
        if isempty(f) || f(end) - f(1) > -1e-4
            disp('Converged!');
            converged = true;
        end
        
        % Project A back onto the cone of PSD matrices
        [vec, val] = eig(A);
        val = diag(val);
        ind = find(val > 0);
        A = real(bsxfun(@times, val(ind)', vec(:,ind)) * vec(:,ind)');
    end 
    
    % Obtain low-dimensional representation
    [mapping.M, foo, bar] = svd(A);
    mapping.M = mapping.M(:,1:no_dims);    
