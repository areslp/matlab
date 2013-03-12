function mappedX = diffusion_maps(X, no_dims, t, sigma)
%DIFFUSION_MAPS Runs the diffusion map algorithm
%
%   mappedX = diffusion_maps(X, no_dims, t, sigma)
%
% The functions runs the diffusion map algorithm on dataset X to reduce it 
% to dimensionality no_dims. The variable sigma is the variance of the Gaussian
% used in the affinity computation (default = 1). The variable alpha
% determines the operator that is applied on the graph (default = 1).
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    % Give memory warning
    if size(X, 1) > 3000
        warning(['Due to the large number of instances (' num2str(size(X, 1)) '), diffusion maps may run out of memory.']);
    end
    
    % Normalize data
    X = double(X);
    X = X - min(X(:));
    X = X / max(X(:));

    % Compute Gaussian kernel matrix
    disp(['Compute Markov forward transition probability matrix with ' num2str(t) ' timesteps...']);
    sumX = sum(X .^ 2, 2);
    K = exp(-bsxfun(@plus, sumX, bsxfun(@plus, sumX', -2 * (X * X'))) ./ (2 .* sigma ^ 2));
    
    % Compute Markov probability matrix with t timesteps
    p = sum(K, 1)';
    K = K ./ ((p * p') .^ t);
    p = sqrt(sum(K, 1))';
    K = K ./ (p * p');
    
    % Perform economy-size SVD
    disp('Perform eigendecomposition...');
    [U, S, V] = svd(K, 0);
    U = bsxfun(@rdivide, U, U(:,1));
    mappedX = U(:,2:no_dims + 1);
    