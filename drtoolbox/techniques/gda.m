function mappedX = gda(X, Y, no_dims, varargin)
%GDA Perform Generalized Discriminant Analysis
%
%	mappedX = gda(X, Y, no_dims)
%	mappedX = gda(X, Y, no_dims, kernel)
%	mappedX = gda(X, Y, no_dims, kernel, param1)
%	mappedX = gda(X, Y, no_dims, kernel, param1, param2)
%
% Perform Generalized Discriminant Analysis. GDA or Kernel LDA is the 
% nonlinear generalization of LDA by means of the kernel trick. X is the
% data on which to perform GDA, Y are the corresponding labels.
% The value of kernel determines the used kernel. Possible values are 'linear',
% 'gauss', 'poly', 'subsets', or 'princ_angles' (default = 'gauss'). For
% more info on setting the parameters of the kernel function, type HELP
% GRAM.
% The function returns the locations of the embedded trainingdata in 
% mappedX.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    % Process inputs
    if ~exist('no_dims', 'var')
        no_dims = 2;
    end
    kernel = 'gauss';
    param1 = 1;
	param2 = 0;
    if length(varargin) > 0 & strcmp(class(varargin{1}), 'char'), kernel = varargin{1}; end 
	if length(varargin) > 1 & strcmp(class(varargin{2}), 'double'), param1 = varargin{2}; end
	if length(varargin) > 2 & strcmp(class(varargin{3}), 'double'), param2 = varargin{3}; end
    
	% Make sure labels are nice
	[foo, bar, Y] = unique(Y, 'rows');

	% Get dimensions
	[n, dim] = size(X);
	nclass = max(Y);

	% Sort data according to labels
	[foo, ind] = sort(Y);
	Y = Y(ind);
	X = X(ind,:);

	% Compute kernel matrix
    disp('Computing kernel matrix...');
    K = gram(X, X, kernel, param1, param2);

    % Compute centering matrix
    ell = size(X, 1);
    D = sum(K) / ell;
    E = sum(D) / ell;
    J = ones(ell, 1) * D;
    K = K - J - J' + E * ones(ell, ell);

    % Perform eigenvector decomposition of kernel matrix (Kc = P * gamma * P')
    disp('Performing eigendecomposition of kernel matrix...');
    K(isnan(K)) = 0;
    K(isinf(K)) = 0;
    [P, gamma] = eig(K);

	if size(P, 2) < n
		error('Singularities in kernel matrix prevent solution.');
    end
	
	% Sort eigenvalues and vectors in descending order
	[gamma, ind] = sort(diag(gamma), 'descend');
	P = P(:,ind);

	% Remove eigenvectors with relatively small value
	minEigv = max(gamma) / 1e5;
	ind = find(gamma > minEigv);
	P = P(:,ind);
	gamma = gamma(ind);
	rankK = length(ind);
    
	% Recompute kernel matrix
    K = P * diag(gamma) * P';

	% Construct diagonal block matrix W
	W = [];
	for i=1:nclass
		num_data_class = length(find(Y == i));
		W = blkdiag(W, ones(num_data_class) / num_data_class);
	end  

	% Determine target dimensionality of data 
	old_no_dims = no_dims;
    no_dims = min([no_dims rankK nclass]);
    if old_no_dims > no_dims
        warning(['Target dimensionality reduced to ' num2str(no_dims) '.']);
    end

	% Perform eigendecomposition of matrix (P' * W * P)
	disp('Performing GDA eigenanalysis...');
    [Beta, lambda] = eig(P' * W * P);
	lambda = diag(lambda);
	
	% Sort eigenvalues and eigenvectors in descending order
	[lambda, ind] = sort(lambda, 'descend');
	Beta = Beta(:,ind(1:no_dims));

	% Compute final embedding mappedX
	mappedX = P * diag(1 ./ gamma) * Beta;

	% Normalize embedding
	for i=1:no_dims
		mappedX(:,i) = mappedX(:,i) / sqrt(mappedX(:,i)' * K * mappedX(:,i));
	end
