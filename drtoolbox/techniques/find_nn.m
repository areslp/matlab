function [D, ni] = find_nn(X, k)
%FIND_NN Finds k nearest neigbors for all datapoints in the dataset
%
%	[D, ni] = find_nn(X, k)
%
% Finds the k nearest neighbors for all datapoints in the dataset X.
% In X, rows correspond to the observations and columns to the
% dimensions. The value of k is the number of neighbors that is
% stored. The function returns a sparse distance matrix D, in which
% only the distances to the k nearest neighbors are stored. For
% equal datapoints, the distance is set to a tolerance value.
% The method is relatively slow, but has a memory requirement of O(nk).
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


	if ~exist('k', 'var') || isempty(k)
		k = 12;
    end
    
    % Perform adaptive neighborhood selection if desired
    if ischar(k)
        [D, max_k] = find_nn_adaptive(X);
        ni = zeros(size(X, 1), max_k);
        for i=1:size(X, 1)
            tmp = find(D(i,:) ~= 0);
            tmp(tmp == i) = [];
            tmp = [tmp(2:end) zeros(1, max_k - length(tmp) + 1)];
            ni(i,:) = tmp;
        end
    
    % Perform normal neighborhood selection
    else

        % Memory conservative implementation
        k = k + 1;
        if size(X, 1) > 2000
            X = X';
            n = size(X, 2);
            D = zeros(n, k);
            XX = sum(X .^ 2, 1);
            onez = ones(1,n);
            if nargout > 1, ni = zeros(n, k, 'uint16'); end
            for i=1:n
                p = X(:,i);
                xx = sum(p .^ 2);
                xX = p' * X;
                d = bsxfun(@plus, XX - 2 * xX, xx);
                [d, ind] = sort(d);
                d = sqrt(d(1:k));
                ind = ind(1:k);
                d(d == 0) = 1e-7;
                D(i,:) = d;
                ni(i,:) = ind;
            end
            D = sparse(repmat((1:size(ni, 1))', [1 size(ni, 2)]), double(ni(:)), double(D(:)), size(ni, 1), size(ni, 1));

        % Faster implementation
        else
            n = size(X, 1);
			D = L2_distance(X', X');
            [foo, ind] = sort(D, 2);
            flat = repmat((1:n)', 1, n - k) + n * ind(:,k+1:end) - n;
            D(flat(:)) = 0;
            D(1:n+1:end) = 1e-7;
            D = sparse(double(D));

            if nargout > 1, ni = uint16(ind(:,1:k)); end
        end
        if nargout > 1
            ni = ni(:,2:end);
        end
    end