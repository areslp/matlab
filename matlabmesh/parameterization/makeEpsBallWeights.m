function [ weights ] = makeEpsBallWeights( P, mode, use_eps, options, do_normalize )
% [ weights ] = makeEpsBallWeights( P, mode, use_eps, options, do_normalize )
%  compute weights for epsilon-ball neighbours
%   [mode] determines weight type. valid settings:
%     'uniform'   - uniform weights   [Tutte]
%     'invdist'   - inverse distance weights
%     'invdist2'  - squared inverse distance weights
%     'gaussian'  - heat-kernel weight e^(-dist^2/4t) 
%          options.t = size        
%   [use_eps] epsilon-ball radius
%   [options] type-specific parameters
%   [do_normalize] set to 1 to normalize weights (default 0)
%   [Ryan Schmidt  rms@dgp.toronto.edu  07/2009]

% compute NxN distance matrix (ick!)
[N,dim] = size(P);
X = P';
X2 = sum(X.^2,1);
distance = sqrt( repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X );

weights = sparse(N,N);

% estimate edge length using some random points
nbrs = unique(ceil(rand(50,1)*N));
sampled = distance(nbrs,:);
sampled2 = sort(sampled,2);
nearestk = sampled2(:,2:7);
edgelen_avg = mean(nearestk(:));
edgelen_max = mean(max(nearestk,[],2));

if ~exist('do_normalize','var')
    do_normalize = 0;
end
if ~exist('use_eps','var')
    use_eps = edgelen_max;
end
if ~exist('options','var')
    options = [];
end
if ~isfield(options,'t')
    options.t = use_eps/100;
end


if strcmp(mode, 'uniform')
    for v = 1:N
        j = find(distance(:,v) < use_eps);
        j = j(j~=v);
        weights(v,j) = 1;
    end
elseif strcmp(mode, 'invdist')
    for v = 1:N
        j = find(distance(:,v) < use_eps);
        j = j(j~=v);
        dists = vmag( vadd(P(j,:), -P(v,:)) );
        weights(v, j) = 1 ./ ( dists + eps );
    end
elseif strcmp(mode, 'invdist2')
    for v = 1:N
        j = find(distance(:,v) < use_eps);
        j = j(j~=v);
        dists = vmag2( vadd(P(j,:), -P(v,:)) );
        weights(v, j) = 1 ./ ( dists + eps );
    end
elseif strcmp(mode, 'gaussian')
    for v = 1:N
        j = find(distance(:,v) < use_eps);
        j = j(j~=v);
        dists2 = vmag2( vadd(P(j,:), -P(v,:)) );
        weights(v,j) = exp( -dists2 / (4*options.t) );
    end   
else
    fprintf('[makeEpsBallWeights] FAILURE: unknown weight mode %s\n', mode);
end


if do_normalize
    for v = 1:N
        j = find(weights(v,:)~=0);
        weights(v, j) = weights(v, j) / sum(weights(v, j));  
    end
end

end