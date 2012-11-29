function [ uv, geo_dists, options_cache ] = expmap( points, normals, vtx, options )
% [ uv, geo_dists, options_cache ] = expmap( points, normals, vtx, options)
% compute discrete exponential
% points, normals are lists of points and normals
% vtx is vertex you want expmap around (arbitrary points not supported)
% options is structure of parameters
%  options.graphnbrs : size of graph nbrhoods that Dijkstra will run on (k-NN)
%  options.stopcrit  : early stopping criteria
%    'none'    - run to completion
%    'k'       - K points have been parameterized
%    'maxdist' - maximum (dijkstra) geodesic distance
%  options.stopparam : parameter for stopping criteria
%    ( count for 'k', approx geodesic distance for 'maxdist' )
%  options.upwindavg : use upwind averaging (more robust - default is 1)
%
% return values:
% uv is per-point uv value
%    - uv's of points outside stopcrit are (inf,inf)
% geo_dists is approx dijkstra geodesic distance from vtx to each point
%    - values outside stopcrit are inf
% options_cache is copy of options with cached values 
%      (pass as argument next time to get better performance)
%
%   [Ryan Schmidt  rms@dgp.toronto.edu  03/2009]

do_profiling = 0;


if ~ exist('options', 'var')
    options = struct('graphnbrs',8,'stopcrit','none','stopparam',inf); 
end
if ~ isfield(options,'graphnbrs')
     options.graphnbrs = 8;
end
if ~ isfield(options,'stopcrit')
     options.nbrtype = 'none'; 
end
if ~ isfield(options,'stopparam')
    if strcmp( options.stopcrit, 'k' )
        options.stopparam = inf;
    else
        error('default setting for stopparam with stopcrit=maxdist is not implemented');
    end
end
if ~ isfield(options,'upwindavg')
    options.upwindavg = 1;
end

% initialize cache if it is not available
have_cache = 0;
if isfield(options,'cache')
    have_cache = 1;
else
    options.cache = struct('G',[],'tan1',[],'tan2',[]);
end


% number of per-vertex neighbours in dijkstra graph
Kgraph = options.graphnbrs;
[dim,N] = size(points');


% make local neighbour distance matrix (does double-duty as connectivity graph)
if do_profiling tic; end
if ~ have_cache
    options.cache.G = findKNN(points,Kgraph);
%    X = points';
%    X2 = sum(X.^2,1);
%    distance = sqrt( repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X );
%    options.cache.G = spalloc(N, N, Kgraph * N);
%    [sorted,index] = sort(distance);
%    for jj = 1:N
%        nbrs = index(2:(1+Kgraph),jj);
%       nbrdists = sorted(2:(1+Kgraph),jj);
%        options.cache.G(jj,nbrs) = nbrdists;
%   end
   
    % [RMS] symmetrize G matrix (should we do this?)
    options.cache.G = max(options.cache.G, options.cache.G');
end
if do_profiling fprintf('%f - make distance matrix\n', toc); end


% early-out termination criteria
stop_dist = inf;
stop_count = N;
if strcmp(options.stopcrit,'k')
    stop_count = min(N, options.stopparam + 1);
%    fprintf('[expmap.m] [TODO] need to set reasonable stop distance threshold...');
    stop_dist = inf;
elseif strcmp(options.stopcrit,'maxdist')
    stop_dist = options.stopparam;
end

% compute distances using Dijkstra's algorithm
if do_profiling tic; end
dist_thresh = 2*stop_dist;    
[geo_dists,Nearest] = dijkstra_maxdist(options.cache.G, vtx, dist_thresh);
geo_dists = geo_dists';
if do_profiling fprintf('%f - dijkstra\n', toc); end

if do_profiling tic; end
% determine embedding order
[Didx,idx] = sort(geo_dists);

% initialize tangent frames
if ~ have_cache
    options.cache.tan1 = zeros(size(points));
    options.cache.tan2 = zeros(size(points));
    for i = 1:N
        [e1,e2] = tangentFrame( normals(i,:) );
        options.cache.tan1(i,:) = e1;
        options.cache.tan2(i,:) = e2;
    end
end

% all uv's are initialized to inf (invalid value)
uv = ones(N,2) * inf;

% initialize seed point
si = idx(1);
uv(si,:) = [0,0];
ns = normals(si,:);
t1s = options.cache.tan1(si,:);
if do_profiling fprintf('%f - tangent frames, etc\n', toc); end

% [TODO] write mex version of this? 
% embed each vertex in order of increasing dijkstra distance
if do_profiling tic; end
for ii = 2:stop_count
    i = idx(ii);
    
    if geo_dists(i) > stop_dist
        break;
    end
    
    if ~ options.upwindavg
  
        % nearest upwind nbr
        ui = Nearest(i);
        nu = normals(ui,:);
        t1u = options.cache.tan1(ui,:);

        % project nbr vector into tangent frame at upwind nbr (preserve length)
        v = points(i,:) - points(ui,:);
        len = vmag(v);
        t2u = options.cache.tan2(ui,:);
        localuv = [ vdot(v,t1u), vdot(v,t2u) ];
        localuv = len * normalize(localuv);

        % [RMS] this is backwards....try going the other way?
        %    (also comments are wrong)

        % rotate from frame at unbr into seed frame
        nrot = valign( ns, nu );
        rt1s = mvmul(nrot,t1s);
        % rotate around normal at seed point to align tangent frames
        [trot,axis,angle] = valign(rt1s, t1u);
        if vdot(axis,nu) < 0
            angle = -angle;
        end

        % construct 2D rotation in tangent space at upwind nbr
        cosA = cos(angle);
        sinA = sin(angle);
        frameRot = [cosA, sinA; -sinA, cosA];

        % apply rotation
        localuv = mvmul(frameRot,localuv);
        uv_u = localuv(1);
        uv_v = localuv(2);

        % accumulate vector in tangent space
        uv(i,:) = uv(ui,:) + [uv_u, uv_v];
        
    else
        % ui = set of upwind nbrs connected to i that are already fixed
        unbrs = find(options.cache.G(i,:)~=0);
        unbrs = unbrs(uv(unbrs,1)~=Inf);
        weightsum = 0;
        sumuv = [0,0];
        for jj = 1:numel(unbrs)
            ui = unbrs(jj);
            weight = 1 / ( vmag2(points(i,:) - points(ui,:)) + eps );
            uv_jj = estimateUV( i, ui, si, points, normals, uv, options.cache.tan1, options.cache.tan2 );
            sumuv = sumuv + weight * uv_jj;
            weightsum = weightsum + weight;
        end
        uv(i,:) = sumuv / weightsum;
        
    end
    
end
if do_profiling fprintf('%f - iteration\n', toc); end

% [TODO] can remove this after we start passing stopcrit to dijkstra algo
geo_dists( uv(:,1) == inf ) = inf;

options_cache = options;

end






function [uv_est] = estimateUV(i, ui, si, points, normals, uv, tan1, tan2 )

    ns = normals(si,:);
    t1s = tan1(si,:);
    
    nu = normals(ui,:);    
    t1u = tan1(ui,:);
    t2u = tan2(ui,:);
    
    v = points(i,:) - points(ui,:);
    len = vmag(v);
    localuv = [ vdot(v,t1u), vdot(v,t2u) ];
    localuv = len * normalize(localuv);

    % [RMS] this is backwards....try going the other way?
    %    (also comments are wrong)

    % rotate from frame at unbr into seed frame
    nrot = valign( ns, nu );
    rt1s = mvmul(nrot,t1s);
    % rotate around normal at seed point to align tangent frames
    [trot,axis,angle] = valign(rt1s, t1u);
    if vdot(axis,nu) < 0
        angle = -angle;
    end

    % construct 2D rotation in tangent space at upwind nbr
    cosA = cos(angle);
    sinA = sin(angle);
    frameRot = [cosA, sinA; -sinA, cosA];

    % apply rotation
    localuv = mvmul(frameRot,localuv);
    uv_u = localuv(1);
    uv_v = localuv(2);

    % accumulate vector in tangent space
    uv_est = uv(ui,:) + [uv_u, uv_v];
end
