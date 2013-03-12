function [mappedX, mapping] = lpp(X, no_dims, k, sigma, eig_impl)
%LPP Perform linearity preserving projection
%
%   [mappedX, mapping] = lpp(X, no_dims, k, sigma, eig_impl)
%
% Perform the Linearity Preserving Projection on dataset X to reduce it to 
% dimensionality no_dims. The number of neighbors that is used by LPP is
% specified by k (default = 12). The variable sigma determines the
% bandwidth of the Gaussian kernel (default = 1).
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

    % Compute XDX and XLX and make sure these are symmetric
    disp('Computing low-dimensional embedding...');
    DP = X' * D * X;
    LP = X' * L * X;
    DP = (DP + DP') / 2;
    LP = (LP + LP') / 2;

    % Perform eigenanalysis of generalized eigenproblem (as in LEM)
    if size(X, 1) > 200 && no_dims < (size(X, 1) / 2)
        if strcmp(eig_impl, 'JDQR')
            options.Disp = 0;
            options.LSolver = 'bicgstab';
            [eigvector, eigvalue] = jdqz(LP, DP, no_dims, 'SA', options);
        else
            options.disp = 0;
            options.issym = 1;
            options.isreal = 1;
            [eigvector, eigvalue] = eigs(LP, DP, no_dims, 'SA', options);
        end
    else
        [eigvector, eigvalue] = eig(LP, DP);
    end
    
    % Sort eigenvalues in descending order and get smallest eigenvectors
    [eigvalue, ind] = sort(diag(eigvalue), 'ascend');
    eigvector = eigvector(:,ind(1:no_dims));
    
    % Compute final linear basis and map data
    mappedX = X * eigvector;
    mapping.M = eigvector;
    mapping.mean = mean(X, 1);
