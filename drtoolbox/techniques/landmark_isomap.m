function [mappedX, mapping] = landmark_isomap(X, no_dims, k, percentage)
%ISOMAP Runs the Isomap algorithm
%
%   [mappedX, mapping] = landmark_isomap(X, no_dims, k, percentage); 
%
% The functions runs the Landmark Isomap algorithm on dataset X to reduce the
% dimensionality of the dataset to no_dims. The number of neighbors used in
% the compuations is set by k (default = 12). The variable percentage has to be
% between 0 and 1, and determines the number of landmarks that is used (default = 0.1).
%
% If the neighborhood graph that is constructed is not completely
% connected, only the largest connected component is embedded. The indices
% of this component are returned in mapping.conn_comp.
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
	if ~exist('percentage', 'var')
		percentage = 0.2;
    end

    % Construct neighborhood graph
    disp('Constructing neighborhood graph...'); 
    D = find_nn(X, k);
    
    % Select largest connected component
    blocks = components(D)';
    count = zeros(1, max(blocks));
    for i=1:max(blocks)
        count(i) = length(find(blocks == i));
    end
    [count, block_no] = max(count);
    conn_comp = find(blocks == block_no);    
    D = D(conn_comp, conn_comp);
    mapping.D = D;
    n = size(D, 1);

    % Compute shortest paths
    disp('Computing shortest paths...');
    landmarks = randperm(n);
	landmarks = landmarks(1:round(percentage * n));
    nl = length(landmarks); 
    D = dijkstra(D, landmarks);
	D = full(D)';
    mapping.DD = D;
    mapping.landmarks = landmarks;

	% Do not embed in more dimensions than (nl - 1)
	if no_dims > nl - 1
		no_dims = nl - 1;
		warning(['Target dimensionality reduced to ' num2str(no_dims) '...']);
	end
    
    % Performing MDS using eigenvector implementation
    disp('Constructing low-dimensional embedding...'); 
    D = D .^ 2;
    subB = -.5 .* (bsxfun(@minus, bsxfun(@minus, D, sum(D, 2) ./ nl), sum(D, 1) ./ n) + sum(D(:)) ./ (n * nl));
    subB2 = subB' * subB;
	subB2(isnan(subB2)) = 0;
	subB2(isinf(subB2)) = 0;
    [alpha, beta] = eig(subB2);
    val = beta .^ (1 / 2);
    mapping.invVal = inv(val);
    vec = subB * alpha * mapping.invVal;
	if size(vec, 2) < no_dims
		no_dims = size(vec, 2);
		warning(['Target dimensionality reduced to ' num2str(no_dims) '...']);
	end
	
    % Computing final embedding
    [val, ind] = sort(real(diag(val)), 'descend'); 
    vec = vec(:,ind(1:no_dims));
    val = val(1:no_dims);
    mappedX = real(bsxfun(@times, vec, sqrt(val)'));
    
    % Store data for out-of-sample extension
    mapping.conn_comp = conn_comp;
    mapping.k = k;
    mapping.X = X(conn_comp,:);
    mapping.alpha = alpha;
    mapping.beta = beta;
    mapping.no_dims = no_dims;
    
    