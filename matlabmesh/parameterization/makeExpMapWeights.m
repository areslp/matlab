function [ W, err2, vtx_u, vtx_v ] = makeExpMapWeights( pointset, options )
% [ W, err2, vtx_u, vtx_v ] = makeExpMapWeights( pointset, options )
%   compute expmap weights for each point with associated normal n in pointset
% options.weightmode : determines weight type. valid settings:
%    'uniform'  - uniform weights
%    'deluniform' - uniform weights for connected neighbours in 
%                   tangent-space delaunay triangulation
%    'invdist'  - inverse distance weights
%    'optimal2' - optimal weights from uvs (solve linear system)
%    'optimal3' - optimal weights from points (solve linear system)
%    'gaussian2' - gaussian heat-kernel weights from uvs (from Laplacian Eigemaps)
%    'gaussian3' - gaussian heat-kernel weights from points (from Laplacian Eigemaps)
% options.nbrtype : determines type of nbrhood
%    'k'       - K-nearest-geodesic neighbourhood
%    'geoball' - geodesic-distance ball
% options.nbrsize : size parameter for nbrhood
%    ( count for 'k', radius for 'eball' )
% options.vnbrsize : per-vertex nbrhood size (ignored by default)
%
% options.heath    : h/t parameter for heat kernel gaussians
% options.regtol   : regularization tolerance for optimal weights 
%                       (default 1e-3)
%
% options.normalize : normalize weights so they sum to 1  
%                       (default true for optimal weights, false for rest)
%
% options.silent   : set to 1 to get no output (default 0)
%   [Ryan Schmidt  rms@dgp.toronto.edu  03/2009]

N = size(pointset.v,1);

if ~ exist('options', 'var')
    options = struct('weightmode','optimal2','nbrtype','k','nbrsize',15); 
end
if ~ isfield(options,'weightmode')
     options.weightmode = 'optimal2'; 
end
if ~ isfield(options,'nbrtype')
     options.nbrtype = 'k'; 
end
if ~ isfield(options,'nbrsize')
    if strcmp( options.nbrtype, 'k' )
        options.nbrsize = 15;
    else
        error('default setting for nbrsize with nbrtype=geoball is not implemented');
    end
end
if ~ isfield(options,'regtol')
    options.regtol = 10e-3;
end
if ~ isfield(options,'normalize')
    if strcmp(options.weightmode,'optimal2') | strcmp(options.weightmode,'optimal3') | strcmp(options.weightmode,'optimal2b')
        options.normalize = 1;
    else
        options.normalize = 0;
    end
end
be_silent = 0;
if isfield(options,'silent')
    be_silent = options.silent;
end
if strcmp(options.weightmode,'gaussian2') | strcmp(options.weightmode,'gaussian3')
   if ~isfield(options,'heath') & strcmp(options.nbrtype,'geoball')
       options.heath = options.nbrsize / 3;
   elseif ~isfield(options,'heath')
       error('must either pass geoball size or radius size .heath for heat-kernel weights');
   end
end


vnbrsize = ones(N,1) * options.nbrsize;
if isfield(options,'vnbrsize') && numel(options.vnbrsize) == N
    vnbrsize = options.vnbrsize;
end


use_fast_expmap = 1;

fprintf('[makeExpMapWeights] computing expmaps\n');

% construct (nbrs,uv) for each vertex
E = sparse(N,N);        % [TODO] better pre-allocation for speed...
if use_fast_expmap
    enable_W_cache = 0;
    fprintf('   using expmapCL.exe\n');

    fake_mesh = struct('v',pointset.v, 'n', pointset.n, 'u',[], 'fidx',[], 'vidx',1:N);

    if strcmp(options.nbrtype,'k')
        [uu, vv] = fastExpMaps(fake_mesh, 'k', options.nbrsize);
    else
        [uu, vv] = fastExpMaps(fake_mesh, 'h', options.nbrsize);
    end
    

    geo_dists = sqrt(uu.^2+vv.^2);
    U = uu;
    V = vv;
    %E(find(geo_dists(:)~=0)) = 1;
    Ecount = nnz(geo_dists);
    
else
    fprintf('   using expmap.m\n');
    
    % use cached version if possible (opens much faster)
    enable_W_cache = globalConfig('enableWeightMatrixCaching');
    hash = meshHash(pointset.v, points.n);
    modestring = ['makeExpMapWeights','_',num2str(hash),'_',options.weightmode,'_',options.nbrtype,'_',num2str(round(abs(log(options.nbrsize))*100)),'_',num2str(round(abs(log(options.regtol))*100))];
    cachepath = [globalConfig('cachePath'), modestring];
    if enable_W_cache & exist(cachepath)
        VARS = load(cachepath, '-mat');
        options2 = VARS.options;
        W = VARS.W;
        err2 = VARS.err2;
        vtx_u = VARS.vtx_u;
        vtx_v = VARS.vtx_v;
        % do version checking?
        fprintf('[makeExpMapWeights] loaded cached weights %s\n', cachepath);    
        return;
    end


    % construct options settings for expmap.m
    nbrsize_fudge = 1;      
    if strcmp(options.nbrtype,'k')
        expmap_options.stopcrit = 'k';
    else
        expmap_options.stopcrit = 'maxdist';
        nbrsize_fudge = 1.1;    % deal w/ roundoff problems, etc
    end

    for v = 1:N
        if ~be_silent & mod(v,100) == 0
            fprintf('%d / %d finished (%3.2f%%)\n', v, N, (v/N)*100);
        end

        nbrsize = vnbrsize(v);
        expmap_options.stopparam = nbrsize * nbrsize_fudge;
        
        % save and re-use cache after first run, but don't overwrite each time
        if v == 1
            [uv,geo_dists,expmap_options] = expmap(pointset.v, points.n, v, expmap_options);
        else
            [uv,geo_dists] = expmap(pointset.v, points.n, v, expmap_options);
        end

        % find neighbours for point v
        nbrs = find(geo_dists < Inf);
        nbrs = nbrs(nbrs ~= v);

        %E(v,nbrs) = 1;
        U(v,nbrs) = uv(nbrs,1);
        V(v,nbrs) = uv(nbrs,2);
    end
    Ecount = nnz(U);
end


fprintf('[makeExpMapWeights] computing weights\n');

% return UV coordinates
vtx_u = U;
vtx_v = V;

% regularization factor for optimal weights
reg_tol = options.regtol;

% now construct weights for each vertex

% precompute strcmp statements (helps speed a lot)
use_k_nbrs = strcmp(options.nbrtype , 'k');
weight_mode = 0;
if strcmp(options.weightmode , 'uniform')
    weight_mode = 1;
elseif strcmp(options.weightmode, 'deluniform')
    weight_mode = 2;
elseif strcmp(options.weightmode , 'invdist')
    weight_mode = 3;
elseif strcmp(options.weightmode, 'gaussian2')
    weight_mode = 4;
elseif strcmp(options.weightmode, 'gaussian3')
    weight_mode = 5;
elseif strcmp(options.weightmode , 'optimal2')
    weight_mode = 6;
elseif strcmp(options.weightmode , 'optimal3')
    weight_mode = 7;
elseif strcmp(options.weightmode, 'optimal2b')
    weight_mode = 8;
end

sum_matrix_time = 0;
sum_setup_time = 0;
all_nbr_dists2 = U.^2+V.^2;
err2 = zeros(N,1);
%W = sparse([],[],[],N,N,nnz(E));
W = sparse([],[],[],N,N,Ecount);
for v = 1:N 

    tic;
    if ~be_silent & mod(v,500) == 0
        fprintf('    %d / %d finished (%3.2f%%)\n', v, N, (v/N)*100);
    end
    
    % find uv's for this point
    nbrsize = vnbrsize(v);
    if use_k_nbrs
        nbrs = find(all_nbr_dists2(v,:)>0);
        dists2 = full(all_nbr_dists2(v,nbrs));
        [vsort,isort] = sort(dists2);
        if numel(nbrs) > nbrsize
            nbrs = nbrs(isort(1:nbrsize));
        end
        uv = full([U(v,nbrs)', V(v,nbrs)']);
    else
        nbrs =  find(all_nbr_dists2(v,:) > 0 & all_nbr_dists2(v,:) < nbrsize*nbrsize)';
        uv = full([U(v,nbrs)', V(v,nbrs)']);
    end
    ncount = numel(nbrs);
    
    % compute weights for neighbours
    setuptime = toc;
    sum_setup_time = sum_setup_time + setuptime;
    
    if weight_mode == 1          % 'uniform'
        W(v, nbrs) = 1;
        
    elseif weight_mode == 2      % 'deluniform'
        TRI = delaunay([0, uv(:,1)], [0,uv(:,2)]);
        [rr,cc]=find(TRI==1);
        tris = TRI(rr,:);
        dnbrs = unique(tris(:));
        W(v,nbrs(dnbrs)) = 1;
        
    elseif weight_mode == 3      % 'invdist'
        dists = vmag(uv);
        W(v, nbrs) = 1 ./ ( dists + eps );
        
    elseif weight_mode == 4     %  'gaussian2'
        dists2 = vmag2( uv );
        W(v,nbrs) = exp( -dists2 / (4*options.heath) );

    elseif weight_mode == 5     %  'gaussian3'
        dists3 = vmag2( vadd( pointset.v(nbrs,:), -pointset.v(v,:) ) );
        W(v,nbrs) = exp( -dists3 / (4*options.heath) );
        
    elseif weight_mode == 6     % 'optimal2')
        tic;
        % based on code from Roweis & Saul's lle.m
        z = uv;                             % points already at origin (no need to shift)
        C = z * z';                         % local covariance
        C = C + eye(ncount)*reg_tol*trace(C);   % regularization (K > dim)
        W(v, nbrs) = C \ ones(ncount,1);        % solve Cw = 1

        solvetime = toc;
        sum_matrix_time = sum_matrix_time+solvetime;
        
    elseif weight_mode == 7        % 'optimal3'
        % based on code from Roweis & Saul's lle.m
        z = vadd( pointset.v(nbrs,:), -pointset.v(v,:) );   % shift to origin
        C = z * z';                                 % local covariance
        C = C + eye(ncount)*reg_tol*trace(C);       % regularization (K > dim)
        W(v, nbrs) = C \ ones(ncount,1);        % solve Cw = 1

    elseif weight_mode == 8      %  'optimal2b'

        z = uv;                                     % shift to origin
        C = z * z';                                 % local covariance
        C = C + eye(ncount)*sqrt(reg_tol)*trace(C);       % regularization (K > dim)
        
        %row_weights = full(options.areas(nbrs));
        
        M = [C*2, -ones(ncount,1); ones(1,ncount), 0];
        if 0
            X = M \ [zeros(ncount,1) ; 1];
        else
            % weighted-weights term
            wk = 0.001;  %0.001;
            Mn = size(M,1);
            HM = wk*2*eye(Mn,Mn);
            HM(Mn,Mn) = 0;
            M = M + HM;
            dists = vmag2(z);
            h = 1-exp(-dists/max(dists));
            X = M \ [-2*h*wk ; 1];
        end
        
        W(v,nbrs) = X(1:ncount);
    end
    
    if options.normalize
        W(v, nbrs) = W(v, nbrs) / sum(W(v,nbrs));   % enforce sum(w) = 1
    end
        
    
    % reconstruction error
    newu = vdot( uv', W(v,nbrs) )';
    err2(v) = vmag2(newu);
    
end

if ~be_silent fprintf('      total solve time: %f    matlab setup time %f\n', sum_matrix_time, sum_setup_time); end;

% cache
if enable_W_cache
    save( cachepath, 'options', 'W', 'err2', 'vtx_u', 'vtx_v', '-mat' );
end 







