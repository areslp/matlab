function [mappedX, mapping] = kernel_pca(X, no_dims, varargin)
%KERNEL_PCA Perform the kernel PCA algorithm
%
%   [mappedX, mapping] = kernel_pca(X, no_dims)
%   [mappedX, mapping] = kernel_pca(X, no_dims, kernel)
%   [mappedX, mapping] = kernel_pca(X, no_dims, kernel, param1)
%   [mappedX, mapping] = kernel_pca(X, no_dims, kernel, param1, param2)
%
% The function runs kernel PCA on a set of datapoints X. The variable
% no_dims sets the number of dimensions of the feature points in the 
% embedded feature space (no_dims >= 1, default = 2). 
% For no_dims, you can also specify a number between 0 and 1, determining 
% the amount of variance you want to retain in the PCA step.
% The value of kernel determines the used kernel. Possible values are 'linear',
% 'gauss', 'poly', 'subsets', or 'princ_angles' (default = 'gauss'). For
% more info on setting the parameters of the kernel function, type HELP
% GRAM.
% The function returns the locations of the embedded trainingdata in 
% mappedX. Furthermore, it returns information on the mapping in mapping.
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
    kernel = 'gauss';
    param1 = 1;
	param2 = 3;
    if nargin > 2
		kernel = varargin{1};
		if length(varargin) > 1 & strcmp(class(varargin{2}), 'double'), param1 = varargin{2}; end
		if length(varargin) > 2 & strcmp(class(varargin{3}), 'double'), param2 = varargin{3}; end
    end
    
    % Store the number of training and test points
    ell = size(X, 1);

    if size(X, 1) < 2000

        % Compute Gram matrix for training points
        disp('Computing kernel matrix...'); 
        K = gram(X, X, kernel, param1, param2);

        % Normalize kernel matrix K
        mapping.column_sums = sum(K) / ell;                       % column sums
        mapping.total_sum   = sum(mapping.column_sums) / ell;     % total sum
        J = ones(ell, 1) * mapping.column_sums;                   % column sums (in matrix)
        K = K - J - J';
        K = K + mapping.total_sum;
 
        % Compute first no_dims eigenvectors and store these in V, store corresponding eigenvalues in L
        disp('Eigenanalysis of kernel matrix...');
        K(isnan(K)) = 0;
        K(isinf(K)) = 0;
        [V, L] = eig(K);
    else
        % Compute column sums (for out-of-sample extension)
        mapping.column_sums = kernel_function([], X', 1, kernel, param1, param2, 'ColumnSums') / ell;
        mapping.total_sum   = sum(mapping.column_sums) / ell;
        
        % Perform eigenanalysis of kernel matrix without explicitly
        % computing it
        disp('Eigenanalysis of kernel matrix (using slower but memory-conservative implementation)...');
        options.disp = 0;
        options.isreal = 1;
        options.issym = 1;
        [V, L] = eigs(@(v)kernel_function(v, X', 1, kernel, param1, param2, 'Normal'), size(X, 1), no_dims, 'LM', options);
        disp(' ');
    end
    
    % Sort eigenvalues and eigenvectors in descending order
    [L, ind] = sort(diag(L), 'descend');
    L = L(1:no_dims);
	V = V(:,ind(1:no_dims));
    
    % Compute inverse of eigenvalues matrix L
	disp('Computing final embedding...');
    invL = diag(1 ./ L);
    
    % Compute square root of eigenvalues matrix L
    sqrtL = diag(sqrt(L));
    
    % Compute inverse of square root of eigenvalues matrix L
    invsqrtL = diag(1 ./ diag(sqrtL));
    
    % Compute the new embedded points for both K and Ktest-data
    mappedX = sqrtL * V';                     % = invsqrtL * V'* K
    
    % Set feature vectors in original format
    mappedX = mappedX';
    
    % Store information for out-of-sample extension
    mapping.X = X;
    mapping.V = V;
    mapping.invsqrtL = invsqrtL;
    mapping.kernel = kernel;
    mapping.param1 = param1;
    mapping.param2 = param2;
    