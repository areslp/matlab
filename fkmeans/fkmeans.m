function [label, centroid, dis] = fkmeans(X, k, options)
% FKMEANS Fast K-means with optional weighting and careful initialization.
% [L, C, D] = FKMEANS(X, k) partitions the vectors in the n-by-p matrix X
% into k (or, rarely, fewer) clusters by applying the well known batch
% K-means algorithm. Rows of X correspond to points, columns correspond to
% variables. The output k-by-p matrix C contains the cluster centroids. The
% n-element output column vector L contains the cluster label of each
% point. The k-element output column vector D contains the residual cluster
% distortions as measured by total squared distance of cluster members from
% the centroid.
%
% FKMEANS(X, C0) where C0 is a k-by-p matrix uses the rows of C0 as the
% initial centroids instead of choosing them randomly from X.
%
% FKMEANS(X, k, options) allows optional parameter name/value pairs to 
% be specified. Parameters are:
%
%   'weight' - n-by-1 weight vector used to adjust centroid and distortion
%              calculations. Weights should be positive.
%   'careful' - binary option that determines whether "careful seeding"
%               as recommended by Arthur and Vassilvitskii is used when
%               choosing initial centroids. This option should be used
%               with care because numerical experiments suggest it may
%               be counter-productive when the data is noisy.
%
% Notes
% (1) The careful seeding procedure chooses the first centroid at random
% from X, and each successive centroid from the remaining points according
% to the categorical distribution with selection probabilities proportional
% to the point's minimum squared Euclidean distance from the already chosen
% centroids. This tends to spread the points out more evenly, and, if the
% data is made of k well separated clusters, is likely to choose an initial
% centroid from each cluster. This can speed convergence and reduce the
% likelihood of getting a bad solution [1]. However, in experiments where
% 5% uniformly distributed noise data was added to such naturally clustered
% data the results were frequently worse then when centroids were chosen at
% random.
% (2) If, as is possible, a cluster is empty at the end of an iteration,
% then there may be fewer than k clusters returned. In practice this seems
% to happen very rarely.
% (3) Unlike the Mathworks KMEANS this implementation does not perform a
% final, slow, phase of incremental K-means ('onlinephase') that guarantees
% convergence to a local minimum. 
%
% References
% [1] "k-means++: The Advantages of Careful Seeding", by David Arthur and
% Sergei Vassilvitskii, SODA 2007.

n = size(X,1);

% option defaults
weight = 0; % uniform unit weighting
careful = 0;% random initialization

if nargin == 3
    if isfield(options, 'weight')
        weight = options.weight;
    end
    if isfield(options,'careful')
        careful = options.careful;
    end
end

% If initial centroids not supplied, choose them
if isscalar(k)
    % centroids not specified
    if careful
        k = spreadseeds(X, k);
    else
        k = X(randsample(size(X,1),k),:);
    end
end

% generate initial labeling of points
[~,label] = max(bsxfun(@minus,k*X',0.5*sum(k.^2,2)));
k = size(k,1);

last = 0;

if ~weight
    % code defactoring for speed
    while any(label ~= last)
        % remove empty clusters
        [~,~,label] = unique(label);
        % transform label into indicator matrix
        ind = sparse(label,1:n,1,k,n,n);
        % compute centroid of each cluster
        centroid = (spdiags(1./sum(ind,2),0,k,k)*ind)*X;
        % compute distance of every point to each centroid
        distances = bsxfun(@minus,centroid*X',0.5*sum(centroid.^2,2));
        % assign points to their nearest centroid
        last = label;
        [~,label] = max(distances);
    end
    dis = ind*(sum(X.^2,2) - 2*max(distances)');
else
    while any(label ~= last)
        % remove empty clusters
        [~,~,label] = unique(label);
        % transform label into indicator matrix
        ind = sparse(label,1:n,weight,k,n,n);
        % compute centroid of each cluster
        centroid = (spdiags(1./sum(ind,2),0,k,k)*ind)*X;
        % compute distance of every point to each centroid
        distances = bsxfun(@minus,centroid*X',0.5*sum(centroid.^2,2));
        % assign points to their nearest centroid
        last = label;
        [~,label] = max(distances);
    end
    dis = ind*(sum(X.^2,2) - 2*max(distances)');
end
label = label';

% Code below this line reused from the file exchange submission K-means++
% (http://www.mathworks.com/matlabcentral/fileexchange/28901-k-means) in
% accordance with the license:
% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
% Copyright (c) 2010, Michael Chen
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% 
function D = sqrdistance(A, B)
% Square Euclidean distances between all sample pairs
% A:  n1 x d data matrix
% B:  n2 x d data matrix
% WB: n2 x 1 weights for matrix B
% D: n2 x n1 pairwise square distance matrix
%    D(i,j) is the squared distance between A(i,:) and B(j,:)
% Written by Michael Chen (sth4nth@gmail.com). July 2009.
n1 = size(A,1); n2 = size(B,2);
m = (sum(A,1)+sum(B,1))/(n1+n2);
A = bsxfun(@minus,A,m);
B = bsxfun(@minus,B,m);
D = full((-2)*(A*B'));
D = bsxfun(@plus,D,full(sum(B.^2,2))');
D = bsxfun(@plus,D,full(sum(A.^2,2)))';
end

function [S, idx] = spreadseeds(X, k)
% X: n x d data matrix
% k: number of seeds
% reference: k-means++: the advantages of careful seeding.
% by David Arthur and Sergei Vassilvitskii
% Adapted from softseeds written by Mo Chen (mochen@ie.cuhk.edu.hk), 
% March 2009.
[n,d] = size(X);
idx = zeros(k,1);
S = zeros(k,d);
D = inf(n,1);
idx(1) = ceil(n.*rand);
S(1,:) = X(idx(1),:);
for i = 2:k
    D = min(D,sqrdistance(S(i-1,:),X));
    idx(i) = find(cumsum(D)/sum(D)>rand,1);
    S(i,:) = X(idx(i),:);
end
end

end
