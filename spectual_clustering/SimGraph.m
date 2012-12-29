function [ W ] = SimGraph( M, Type, Param, varparam )
%SIMGRAPH Returns adjacency matrix for M
%   If M is a square n-by-n matrix, SimGraph returns the
%   adjacency matrix for the similarity graph specified by
%   Type and using the parameter Param, where M is assumed to
%   be the distance matrix.
%   If M is not square, it is assumed to be a d-by-n matrix
%   containing data points. Then, SimGraph also returns the
%   adjacency matrix, but computes it more efficiently.
%
%   The last parameter varparam can either be a scalar, defining
%   the size of neighborhoods when a Gaussian similarity function
%   is applied, or it can be a reference to a function that will
%   be used for weighting the similarity graph (use '@' when
%   calling the function).
%
%   'M' - Matrix either a distance matrix or data points
%   'Type' - Defines type of similarity graph to be used
%      1 - Full Similarity Graph
%      2 - Epsilon Similarity Graph
%      3 - kNearest Neighbors Graph
%      4 - Mutual kNeares Neighbors Graph
%   'Param' - In case of an epsilon graph, this defines epsilon.
%      In case of a kNearest neighbors graph, this defines k.
%   'varparam' (optional) - Either a scalar, defining the size of
%      neighborhoods or a function, calculating the similarities
%      from the distance matrix (has to take arbitrary matrices
%      as the only input and return a matrix of the same size).
%
%   Author: Ingo Buerk
%   Year  : 2011/2012
%   Bachelor Thesis

% if matrix is not square, we assume it to be the data point
% matrix and, using sparse matrices, compute the adjacency
% matrix parallel to calculating the pairwise distances, in
% order to save memory

% catch illegal inputs
if ~any(Type == (1:4))
    error('Unknown Similarity Graph type.');
end

% number of data points
n = size(M, 2);

% in case of the full similarity graph, there is no need for
% efficient sparse techniques
if (size(M, 1) ~= n) && Type == 1
    
    % calculate distance matrix
    W = squareform(pdist(M'));
    
    % clear variable from memory
    clear M;
    
    % however, spectral algorithms require similarity
    % matrices rather than distance matrices, so if
    % sigma was given, we apply a Gaussian similarity
    % function
    if nargin == 4 && ~isempty(varparam)
       W = exp(-W.^2 ./ (2*varparam^2)); 
    end
    
elseif (size(M, 1) ~= n) && Type ~= 1
    
    % in cases of epsilon or k-nearest similarity graphs, we
    % save all nonzero elements and their corresponding
    % indices in vectors first, thus making creating the sparse
    % matrix more efficient
    
    % k-nearest similarity graphs
    if Type == 3 || Type == 4
        
        % preallocate memory
        indi = zeros(1, Param * n);
        indj = zeros(1, Param * n);
        inds = zeros(1, Param * n);

        % loop through each data point, calculate the distance
        % to all other data points, sort them and store the
        % k smallest to create the distance matrix
        for ii = 1:n
            dist = distEuclidean(repmat(M(:, ii), 1, n), M);
            
            [s, O] = sort(dist, 'ascend');
            
            indi(1, (ii-1)*Param+1:ii*Param) = ii;
            indj(1, (ii-1)*Param+1:ii*Param) = O(1:Param);
            inds(1, (ii-1)*Param+1:ii*Param) = s(1:Param);
        end
        
        % clear variable from memory
        clear M;
        
        % create the sparse distance matrix for the directed
        % k-nearest neighbors
        W = sparse(indi, indj, inds, n, n);
        
        % clear variables from memory
        clear indi indj inds dist s O;
        
        % from the directed distance matrix, we can get the
        % undirected normal or mutual k-nearest neighbors
        % distance matrix quite simply     
        if Type == 3
            W = max(W, W');
        else
            W = min(W, W');
        end
        
        % since the algorithms take similarity graphs, we need to
        % weight the edges of the graph constructed above with
        % similarity values corresponding to the distances.
        
        % if no argument was given (or 0), we construct an
        % unweighted graph
        if nargin < 4 || isempty(varparam) || varparam == 0
            W = (W ~= 0);
        
        % if the value is a scalar number, we apply a Gaussian
        % neighborhood function
        elseif nargin == 4 && isnumeric(varparam)
            W = spfun(@(W) (simGaussian(W, varparam)), W);
        
        % otherwise, we assume varparam to be a function and call
        % it to calculate the similarity values
        else
            W = spfun(varparam, W);
        end
        
    elseif Type == 2
        
        % in case of the epsilon similarity graph, we can not
        % preallocate any memory, since we don't know how many
        % nonzero elements we are actually going to have
        indi = [];
        indj = [];
        inds = [];
        
        % go through all data points
        for ii = 1:n
            % calculate the pairwise distances between the
            % i-th data point an all data points, which is
            % exactly the i-th column of the distance matrix
            dist = distEuclidean(repmat(M(:, ii), 1, n), M);
            
            % epsilon graphs are not going to be weighted,
            % therefore we only have to set the entries
            % to 0 or 1
            dist = (dist < Param);
            
            lastind     = size(indi, 2);
            count       = nnz(dist);
            [~, col]    = find(dist);
            
            % save indices and values
            indi(1, lastind+1:lastind+count) = ii;
            indj(1, lastind+1:lastind+count) = col;
            inds(1, lastind+1:lastind+count) = 1;
        end
        
        % clear variable from memory
        clear M;
        
        % create the sparse adjacency matrix for the epsilon
        % graph
        W = sparse(indi, indj, inds, n, n);
        
        % clear variables from memory
        clear indi indj inds dist lastind count col v;
    
    end

else
    
    % if the matrix given was square, we assume it is the
    % distance matrix and therefore only calculate the
    % adjacency matrix for the similarity graph
    
    switch Type
        case 1
            % since all vertices will be connected, just use
            % the distance matrix itself
            W = M;
            
            % clear variable from memory
            clear M;
            
            % however, spectral algorithms require similarity
            % matrices rather than distance matrices, so if
            % sigma was given, we apply a Gaussian similarity
            % function
            if nargin == 4 && ~isempty(varparam)
                W = exp(-W.^2 ./ (2*varparam^2));
            end
        case 2
            % epsilon neighborhood similarity graphs are
            % unweighted, therefore we can just use a logical
            % sparse matrix
            W = sparse(M < Param);
            
            % clear variable from memory
            clear M;
        case {3, 4}
            % preallocate memory for nonzero elements and their
            % corresponding indices
            indi = zeros(1, Param * n);
            indj = zeros(1, Param * n);
            inds = zeros(1, Param * n);
            
            % go through each column and pick the first k 
            % neighbors
            for ii = 1:n
                [s, O] = sort(M(:, ii), 'ascend');
                
                indi(1, (ii-1)*Param+1:ii*Param) = ii;
                indj(1, (ii-1)*Param+1:ii*Param) = O(1:Param);
                inds(1, (ii-1)*Param+1:ii*Param) = s(1:Param);
            end
            
            % clear variable from memory
            clear M;
            
            % create directed neighbor graph
            W = sparse(indi, indj, inds, n, n);
            
            % clear variables from memory
            clear indi indj inds s O;
            
            % create normal or mutual kNearest neighbors graph
            if Type == 3
                W = max(W, W');
            else
                W = min(W, W');
            end
            
            % since the algorithms take similarity graphs, we 
            % need to weight the edges of the graph constructed 
            %above with similarity values corresponding to the 
            %distances.
            
            % if no argument was given (or 0), we construct an
            % unweighted graph
            if nargin < 4 || isempty(varparam) || varparam == 0
                W = (W ~= 0);
                
            % if the value is a scalar number, we apply a 
            %Gaussian neighborhood function
            elseif nargin == 4 && isnumeric(varparam)
                W = spfun(@(W) (simGaussian(W, varparam)), W);
                
            % otherwise, we assume varparam to be a function and 
            % call it to calculate the similarity values
            else
                W = spfun(varparam, W);
            end
    end
    
end


end

