function [mappedX, mapping] = fastmvu(X, no_dims, k, finetune, eig_impl)
%FAST_MVU Runs the Fast Maximum Variance Unfolding algorithm
%
%   [mappedX, mapping] = fastmvu(X, no_dims, k, finetune)
%
% Computes a low dimensional embedding of data points using the Fast
% Maximum Variance Unfolding algorithm. The data is specified in an NxD 
% data matrix X. The lowdimensional representation is returned in mappedX.
% Information on the mapping is returned in the mapping struct.
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
    if ~exist('finetune', 'var')
        finetune = true;
    end
    if ~exist('eig_impl', 'var')
        eig_impl = 'Matlab';
    end

    % Initialize some parameters
    maxiter = 10000;                    % maximum number of iterations
    eta = 1e-5 ;                        % regularization parameter
    initial_dims = 2 * no_dims;         % dimensionality of first guess
    
    % Compute sparse distance matrix D on dataset X
    disp('Constructing neighborhood graph...');
    D = find_nn(X, k);
    D(D ~= 0) = 1;
    
    % Select largest connected component
    blocks = components(D)';
    count = zeros(1, max(blocks));
    for i=1:max(blocks)
        count(i) = length(find(blocks == i));
    end
    [count, block_no] = max(count);
    conn_comp = find(blocks == block_no);
    D = D(conn_comp,:);
    D = D(:,conn_comp);
    mapping.conn_comp = conn_comp;
    mapping.D = D;

    % Perform eigendecomposition of the graph Laplacian
    disp('Perform eigendecomposition of graph Laplacian...');
    DD = diag(sum(D, 2));                       % degree matrix
    L = DD - D;                                 % graph Laplacian
    L(isnan(L)) = 0;
	L(isinf(L)) = 0;
    if strcmp(eig_impl, 'JDQR')
        options.Disp = 0;
        options.LSolver = 'bicgstab';
        [laplX, lambda] = jdqr(L, initial_dims + 1, 'SR', options);          % only need bottom (initial_dims + 1) eigenvectors
    else
        options.disp = 0;
        options.isreal = 1;
        options.issym = 1;
        [laplX, lambda] = eigs(L, initial_dims + 1, 'SR', options);			 % only need bottom (initial_dims + 1) eigenvectors
    end
    [lambda, ind] = sort(diag(lambda), 'ascend');
    laplX = real(laplX(:,ind(2:initial_dims + 1)));
    clear DD L

    % Maximize sum of pairwise distances while retaining distances inside
    % neighborhood graph distances. I.e. perform SDP optimization, starting
    % with eigendecomposition of the Laplacian. The constraints of the 
    % optimization are formed by the upper triangle of D.
    disp('Perform semi-definite programming...');
    disp('CSDP OUTPUT =============================================================================');
    try 
        [LL, mappedX, L, newV, idx] = sdecca2(laplX', triu(D), eta, 0);
    catch
        error('Error while performing SDP. Maybe the binaries of CSDP are not suitable for your platform.');
    end
    disp('=========================================================================================');
    
    % Perform initial conjugate gradient search (in initial-dimensional space)
    if finetune
        disp('Finetune initial solution using conjugate gradient descent...');
        mappedX = hillclimber2c(D, mappedX, 'maxiter', maxiter, 'eta', eta);
    end

    % Perform PCA to remove noise and further reduce dimensionality
    disp('Perform PCA to obtain final solution...');
    [mappedX, mapping2] = compute_mapping(mappedX', 'PCA', no_dims);

    % Finetune final solution
    if finetune
        mappedX = mappedX';
        disp('Finetune final solution using conjugate gradient descent...');
        mappedX = hillclimber2c(D, mappedX, 'maxiter', maxiter, 'eta', 0);
        mappedX = mappedX';
    end
    
    % Save data for the out-of-sample extension
    mapping.k = k;
    mapping.L = L;
    mapping.X = X;
    mapping.newV = newV;
    mapping.idx = idx;
    mapping.vec = laplX;
    mapping.pca_map = mapping2;
    mapping.pca_map.name = 'PCA';
    mapping.finetune = finetune;
    