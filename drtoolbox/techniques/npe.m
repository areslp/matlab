function [mappedX, mapping] = npe(X, no_dims, k, eig_impl)
%NPE Perform the Neighborhood Preserving Embedding algorithm
%
%       [mappedX, mapping] = npe(X, no_dims, k)
%       [mappedX, mapping] = npe(X, no_dims, k, eig_impl)
% 
% Runs the Neighborhood Preserving Embedding algorithm on dataset X to 
% reduce it to dimensionality no_dims. The number of neighbors that is used
% by LPP is specified by k (default = 12).
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    if size(X, 2) > size(X, 1)
        error('Number of samples should be higher than number of dimensions.');
    end
    if ~exist('no_dims', 'var')
        no_dims = 2; 
    end
    if ~exist('k', 'var')
        k = 12;
    end
    if ~exist('eig_impl', 'var')
        eig_impl = 'Matlab';
    end
    
    % Get dimensionality and number of dimensions
    [n, d] = size(X);
    mapping.mean = mean(X, 1);

    % Compute pairwise distances and find nearest neighbours (vectorized implementation)
    disp('Finding nearest neighbors...');    
    [distance, neighborhood] = find_nn(X, k);
    max_k = size(neighborhood, 2);
    if nargout > 1
        mapping.nbhd = distance;
    end
    X = X';
    neighborhood = neighborhood';
        
    % Find reconstruction weights for all points by solving the MSE problem 
    % of reconstructing a point from each neighbours. A used constraint is 
    % that the sum of the reconstruction weights for a point should be 1.
    disp('Compute reconstruction weights...');
    if k > d 
        tol = 1e-5;
    else
        tol = 0;
    end

    % Construct reconstruction weight matrix
    W = zeros(max_k, n);
    for i=1:n
        nbhd = neighborhood(:,i);
        if ischar(k)
           nbhd = nbhd(nbhd ~= 0);
        end
        kt = numel(nbhd);
        z = X(:,nbhd) - repmat(X(:,i), 1, kt);                  % Shift point to origin
        C = z' * z;												% Compute local covariance
        C = C + eye(kt, kt) * tol * trace(C);					% Regularization of covariance (if K > D)
        wi = C \ ones(kt, 1);                                   % Solve linear system
        wi = wi / sum(wi);                                      % Make sure that sum is 1
        W(:,i) = [wi; nan(max_k - kt, 1)];
    end

    % Now that we have the reconstruction weights matrix, we define the 
    % sparse cost matrix M = (I-W)'*(I-W).
    M = sparse(1:n, 1:n, ones(1, n), n, n, 4 * max_k * n);
    for i=1:n
       w = W(:,i);
       ww(~isnan(w)) = 0;
       j = neighborhood(:,i);
       w = w(j ~= 0);
       j = j(j ~= 0);
       M(i, j) = M(i, j) - w';
       M(j, i) = M(j, i) - w;
       M(j, j) = M(j, j) + w * w';
    end
	
	% For sparse datasets, we might end up with NaNs or Infs in M. We just set them to zero for now...
	M(isnan(M)) = 0;
	M(isinf(M)) = 0;

    % Compute XWX and XX and make sure these are symmetric
    X = X';
    WP = X' * M * X;
    DP = X' * X;
    DP = (DP + DP') / 2;
    WP = (WP + WP') / 2;

    % Solve generalized eigenproblem
    if size(X, 1) > 1500 && no_dims < (size(X, 1) / 10)
        if strcmp(eig_impl, 'JDQR')
            options.Disp = 0;
            options.LSolver = 'bicgstab';
            [eigvector, eigvalue] = jdqz(WP, DP, no_dims, 'SA', options);
        else
            options.disp = 0;
            options.issym = 1;
            options.isreal = 0;
            [eigvector, eigvalue] = eigs(WP, DP, no_dims, 'SA', options);
        end
    else
        [eigvector, eigvalue] = eig(WP, DP);
    end
    
    % Sort eigenvalues in descending order and get smallest eigenvectors
    [eigvalue, ind] = sort(diag(eigvalue), 'ascend');
    eigvector = eigvector(:,ind(1:no_dims));
    
    % Compute final linear basis and map data
    mappedX = X * eigvector;
    mapping.M = eigvector;
