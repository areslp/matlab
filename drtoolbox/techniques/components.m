function blocks = components(A)
%COMPONENTS Finds connected components in a graph defined by a adjacency matrix
%
%   blocks = components(A)
%
% Finds connected components in a graph defined by the adjacency matrix A.
% The function outputs an n-vector of integers 1:k in blocks, meaning that
% A has k components. The vector blocks labels the vertices of A according 
% to component.
% If the adjacency matrix A is undirected (i.e. symmetric), the blocks are 
% its connected components. If the adjacency matrix A is directed (i.e. 
% unsymmetric), the blocks are its strongly connected components.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    % Check size of adjacency matrix
    [n, m] = size(A);
    if n ~= m, error ('Adjacency matrix must be square'), end;

    % Compute Dulmage-Mendelsohn permutation on A
    if ~all(diag(A)) 
        [foo, p, bar, r] = dmperm(A | speye(size(A)));
    else
        [foo, p, bar, r] = dmperm(A);  
    end

    % Compute sizes and number of clusters
    sizes = diff(r);
    k = length(sizes);

    % Now compute the array blocks
    blocks = zeros(1, n);
    blocks(r(1:k)) = ones(1, k);
    blocks = cumsum(blocks);

    % Permute blocks so it maps vertices of A to components
    blocks(p) = blocks;
