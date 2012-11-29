function [ W ] = makeOneRingWeights( mesh, mode, do_normalize, area_type, options )
% [ W ] = makeOneRingWeights( mesh, mode, normalize, area_type )
%  compute one-right weights for vertex neighbours
%   [mode] determines weight type. valid settings:
%     'uniform'   - uniform weights   [Tutte]
%     'invdist'   - inverse distance weights
%     'invdist2'  - squared inverse distance weights
%     'dcp'       - discrete conformal (cotan) weights [Desbrun02]
%     'cotan'     - same as above
%     'dap'       - discrete authalic weights          [Desbrun02]
%     'meanvalue' - mean-value weights                 [Floater04]
%     'optimal2'  - optimal weights in flattened one-rings
%     'optimal3'  - optimal weights in 3D one-ring
%   [do_normalize] set to 1 to normalize weights (default 0)
%   [area_type] per-vertex area type (rows multiplied by 1/Ai)
%     'uniform'     - Ai = 1   (default)
%     'onering'     - one-ring area
%     'mixed'       - [Meyer02] mixed voronoi weights (see vertexArea())
%   [Ryan Schmidt  rms@dgp.toronto.edu  07/2009]


if ~exist('do_normalize','var')
    do_normalize = 0;
end
if ~exist('area_type','var')
    area_type = 'uniform';
end
if ~ exist('options', 'var')
    options = struct('normalize',0); 
end
if ~ isfield(options,'regtol')
    options.regtol = 10e-3;
end


n = size(mesh.v,1);
W = mesh.e;
A = vertexArea(mesh,[],area_type);

if strcmp(mode, 'uniform')
    for v = 1:n
        [i,j] = find(mesh.e(v,:));
        W(v,j) = 1;
    end
elseif strcmp(mode, 'invdist')
    for v = 1:n
        [i,j] = find(mesh.e(v,:));
        dists = vmag( vadd(mesh.v(j,:), -mesh.v(v,:)) );
        W(v, j) = 1 ./ ( dists + eps );
    end
elseif strcmp(mode, 'invdist2')
    for v = 1:n
        [i,j] = find(mesh.e(v,:));
        dists = vmag2( vadd(mesh.v(j,:), -mesh.v(v,:)) );
        W(v, j) = 1 ./ ( dists + eps );
    end
elseif strcmp(mode, 'gaussian')
    edgelen_avg = edgeStats(mesh.v, mesh.e);
    t = edgelen_avg;
    for v = 1:n
        [i,j] = find(mesh.e(v,:));
        dists2 = vmag2( vadd(mesh.v(j,:), -mesh.v(v,:)) );
        W(v,j) = exp( -dists2 / 4*t );
    end
elseif strcmp(mode, 'optimal2')
    for v = 1:n
        [ringv,vpos] = oneringExpMap(mesh,v);
        ncount = numel(ringv);
        C = vpos * vpos';                         % local covariance
        %C = C + eye(ncount)*trace(C)/100;
        C = C + eye(ncount)*options.regtol*trace(C);   % regularization (K > dim)
        W(v, ringv) = C \ ones(ncount,1);        % solve Cw = 1
        
        %vposw = C \ ones(ncount,1);
        %vposw = vposw / sum(vposw);
        %rpos = sum(repmat(vposw,1,3) .* vpos);
        %fprintf('recons err at v%d is %f\n', v, vmag(rpos));

        % enforce sum-to-1
        W(v, ringv) = W(v, ringv) / sum(W(v,ringv));   % enforce sum(w) = 1
    end
elseif strcmp(mode, 'optimal3')
    for v = 1:n
        [i,j] = find(mesh.e(v,:));
        vpos = vadd(mesh.v(j,:), -mesh.v(v,:));
        ncount = numel(j);
        C = vpos * vpos';                         % local covariance
        C = C + eye(ncount)*options.regtol*trace(C);   % regularization (K > dim)
        W(v, j) = C \ ones(ncount,1);        % solve Cw = 1

        % enforce sum-to-1
        W(v, j) = W(v, j) / sum(W(v,j));   % enforce sum(w) = 1
    end
elseif strcmp(mode, 'dcp')
    W = cotanWeights(mesh);
elseif strcmp(mode, 'cotan')
    W = cotanWeights(mesh);
elseif strcmp(mode, 'dap')
    W = cotanWeights(mesh, [], 1);
elseif strcmp(mode, 'meanvalue')
    W = meanvalueWeights(mesh);
end

for v = 1:n
    W(v,:) = W(v,:) * (1/A(v));
end

if do_normalize
    for v = 1:n
        [i,j] = find(mesh.e(v,:));
        W(v, j) = W(v, j) / sum(W(v, j));  
    end
end

end