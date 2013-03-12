function [mappedX, mapping] = laplacian_eigen(X, no_dims, k, sigma, eig_impl)
%LAPLACIAN_EIGEN Performs non-linear dimensionality reduction using Laplacian Eigenmaps
%
%   [mappedX, mapping] = laplacian_eigen(X, no_dims, k, sigma, eig_impl)
%
% Performs non-linear dimensionality reduction using Laplacian Eigenmaps.
% The data is in matrix X, in which the rows are the observations and the
% columns the dimensions. The variable dim indicates the preferred amount
% of dimensions to retain (default = 2). The variable k is the number of 
% neighbours in the graph (default = 12).
% The reduced data is returned in the matrix mappedX.
%
%

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
    if ~exist('k', 'var')
        k = 12;
    end
	if ~exist('sigma', 'var')
		sigma = 1;
    end
    if ~exist('eig_impl', 'var')
        eig_impl = 'Matlab';
    end
    
    % Construct neighborhood graph
    disp('Constructing neighborhood graph...');
    if size(X, 1) < 4000
        G = L2_distance(X', X');
        
        % Compute neighbourhood graph
        [tmp, ind] = sort(G); 
        for i=1:size(G, 1)
            G(i, ind((2 + k):end, i)) = 0; 
        end
        G = sparse(double(G));
        G = max(G, G');             % Make sure distance matrix is symmetric
    else
        G = find_nn(X, k);
    end
    G = G .^ 2;
    G = G ./ max(max(G));
    
    % Only embed largest connected component of the neighborhood graph
    blocks = components(G)';
    count = zeros(1, max(blocks));
    for i=1:max(blocks)
        count(i) = length(find(blocks == i));
    end
    [count, block_no] = max(count);
    conn_comp = find(blocks == block_no);    
    G = G(conn_comp, conn_comp);
    
    % Compute weights (W = G)
    disp('Computing weight matrices...');
    
    % Compute Gaussian kernel (heat kernel-based weights)
    G(G ~= 0) = exp(-G(G ~= 0) / (2 * sigma ^ 2));
        
    % Construct diagonal weight matrix
    D = diag(sum(G, 2));
    
    % Compute Laplacian
    L = D - G;
    L(isnan(L)) = 0; D(isnan(D)) = 0;
	L(isinf(L)) = 0; D(isinf(D)) = 0;
    
    % Construct eigenmaps (solve Ly = lambda*Dy)
    disp('Constructing Eigenmaps...');
	tol = 0;
    if strcmp(eig_impl, 'JDQR')
        options.Disp = 0;
        options.LSolver = 'bicgstab';
        [mappedX, lambda] = jdqz(L, D, no_dims + 1, tol, options);			% only need bottom (no_dims + 1) eigenvectors
    else
        options.disp = 0;
        options.isreal = 1;
        options.issym = 1;
        [mappedX, lambda] = eigs(L, D, no_dims + 1, tol, options);			% only need bottom (no_dims + 1) eigenvectors
    end
    
    % Sort eigenvectors in ascending order
    lambda = diag(lambda);
    [lambda, ind] = sort(lambda, 'ascend');
    lambda = lambda(2:no_dims + 1);
    
    % Final embedding
	mappedX = mappedX(:,ind(2:no_dims + 1));

    % Store data for out-of-sample extension
    mapping.K = G;
    mapping.vec = mappedX;
    mapping.val = lambda;
    mapping.X = X(conn_comp,:);
    mapping.sigma = sigma;
    mapping.k = k;
    mapping.conn_comp = conn_comp;
    
    