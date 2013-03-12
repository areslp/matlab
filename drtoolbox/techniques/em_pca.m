function [mappedX, mapping] = em_pca(X, no_dims, max_iter)
%EMPCA Run an EM-based implementation of (probabilistic) PCA
%
%   [mappedX, mapping] = em_pca(X, no_dims)
%
% Performs probabilistic PCA on dataset X in order to reduce its
% dimensionality to no_dims. The dimensionality reduction is performed by
% means of an EM algorithm. The resulting low-dimensional counterpart of X
% is returned in mappedX. Information on the applied mapping (allowing for,
% e.g., out-of-sample extension) is returned in mapping.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    if ~exist('max_iter', 'var')
        max_iter = 200;
    end

    % Initialize some variables
    [n D] = size(X);                        % data dimensions
    Ez = zeros(no_dims, n);                 % expectation of latent vars
    Ezz = zeros(no_dims, no_dims, n);       % expectation of cov(z)
    Q = Inf;                                % log-likelihood
    mapping = struct;

    % Randomly initialize W and sigma
    W = rand(D, no_dims) * 2;               % factor loadings
    sigma2 = rand(1) * 2;                   % variance ^ 2
    % The covariance of the Gaussian is: C = W * W' + sigma2 * eye(D);
    
    % Make data zero-mean (possible because data mean is ML estimate for mu)
    mapping.mean = mean(X, 1);
    X = bsxfun(@minus, X, mapping.mean);
    
    % Compute data covariance and transpose data
    S = cov(X);
    X = X';
    
    % Perform EM iterations
    converged = 0;
    iter = 0;
    inW = W' * W;
    while ~converged && iter <= max_iter
            
        % Update iteration number
        iter = iter + 1;
        if rem(iter, 5) == 0
            fprintf('.');
        end
        
        % Perform E-step
        invM = inv(inW + sigma2 * eye(no_dims));
        for i=1:n
            Ez(:,i)    = invM * W' * X(:,i);         
            Ezz(:,:,i) = sigma2 * invM + Ez(:,i) * Ez(:,i)';
        end
        
        % Perform M-step (maximize mapping W)
        Wp1 = zeros(D, no_dims);
        Wp2 = zeros(no_dims, no_dims);
        for i=1:n
            Wp1 = Wp1 + X(:,i) * Ez(:,i)';
            Wp2 = Wp2 + Ezz(:,:,i);
        end
        W = Wp1 / Wp2;
        inW = W' * W;
        
        % Perform M-step (maximize discarded variance sigma)
        normX = sum(X .^ 2, 1);
        sigma2 = 0;
        for i=1:n
            sigma2 = sigma2 + (normX(i) - 2 * Ez(:,i)' * W' * X(:,i) + trace(Ezz(:,:,i) * inW));
        end
        sigma2 = (1 / (n * D)) * sigma2;
        
        % Compute likelihood of new model
        oldQ = Q;
        if iter > 1
            invC = ((1 / sigma2) * eye(D)) - ((1 / sigma2) * W * invM * W');
            detC = det(sigma2 * eye(D)) * det(eye(no_dims) + W' * ((sigma2 .^ -1) * eye(D)) * W);
            Q = (-n / 2) * (D * log(2 * pi) + log(detC) + trace(invC * S));
        end
        
        % Stop condition to detect convergence
        if abs(oldQ - Q) < 1e-3
            converged = 1;
        end
    end
    
    % Compute mapped data
    disp(' ');
    mapping.M = (inW \ W')';
    mapping.sigma2 = sigma2;
    mappedX = X' * mapping.M;    
    