function [ weights ] = makeKNbrWeights( P, mode, use_k, options, do_normalize )
% [ weights ] = makeKNbrWeights( P, mode, use_k, options, do_normalize )
%  compute weights for epsilon-ball neighbours
%   [mode] determines weight type. valid settings:
%     'uniform'   - uniform weights   [Tutte]
%     'invdist'   - inverse distance weights
%     'invdist2'  - squared inverse distance weights
%     'gaussian'  - heat-kernel weight e^(-dist^2/4t) 
%          options.t = size        
%   [use_k] k-nbrhood size
%   [options] type-specific parameters
%   [do_normalize] set to 1 to normalize weights (default 0)
%   [Ryan Schmidt  rms@dgp.toronto.edu  09/2009]

% compute NxN distance matrix (ick!)
[N,dim] = size(P);
X = P';
X2 = sum(X.^2,1);
distance = sqrt( repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X );
[sorted,index] = sort(distance);
all_nbrs = index(2:(1+use_k),:)';


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
    options.t = edgelen_avg;
end


if strcmp(mode, 'uniform')
    for v = 1:N
        j = all_nbrs(v,:);
        weights(v,j) = 1;
    end
elseif strcmp(mode, 'invdist')
    for v = 1:N
        j = all_nbrs(v,:);
        dists = vmag( vadd(P(j,:), -P(v,:)) );
        weights(v, j) = 1 ./ ( dists + eps );
    end
elseif strcmp(mode, 'invdist2')
    for v = 1:N
        j = all_nbrs(v,:);
        dists = vmag2( vadd(P(j,:), -P(v,:)) );
        weights(v, j) = 1 ./ ( dists + eps );
    end
elseif strcmp(mode, 'gaussian')
    for v = 1:N
        j = all_nbrs(v,:);
        dists2 = vmag2( vadd(P(j,:), -P(v,:)) );
        weights(v,j) = exp( -dists2 / (4*options.t) );
    end  
elseif strcmp(mode, 'optimal3')
    for v = 1:N
        reg_tol = 10e-3;
        nbrs = all_nbrs(v,:);
        z = vadd( P(nbrs,:), -P(v,:) );   % shift to origin
        C = z * z';                                 % local covariance
        C = C + eye(size(C,1))*reg_tol*trace(C);       % regularization (K > dim)
        weights(v, nbrs) = C \ ones(size(C,1),1);        % solve Cw = 1

        % enforce sum-to-1
        weights(v, nbrs) = weights(v, nbrs) / sum(weights(v,nbrs));   % enforce sum(w) = 1        
    end      
else
    fprintf('[makeKNbrWeights] FAILURE: unknown weight mode %s\n', mode);
end


if do_normalize
    for v = 1:N
        j = find(weights(v,:)~=0);
        weights(v, j) = weights(v, j) / sum(weights(v, j));  
    end
end

end