function [ mesh ] = readMesh( filename, options )
% [ mesh ] = readMesh( filename, options )
%  Currently only supports .OBJ files.
%  options specifies flags for reading
%   'n' - generate normals if not contained in file
%   'C' - ignore cache
%   'B' - no boundary
%
%   [Ryan Schmidt  rms@dgp.toronto.edu  09/2008]

if ~exist('options','var')
    options = '';
end

enable_mesh_cache = globalConfig('enableMeshCaching');
if ~isempty(findstr(options, 'C'))
    enable_mesh_cache = 0;
end

% check that filename exists
meshpath = which(filename);
if ~exist('meshpath')
    error(['File ',filename,' does not exist in matlab path\n']);
end

% use cached version if possible (opens much faster)
cachepath = [meshpath, '.meshcache'];
if enable_mesh_cache & exist(cachepath)
    S = load(cachepath, '-mat');
    mesh = S.mesh;
    if mesh.version == globalConfig('meshVersion');
        fprintf('[readMesh] loaded cached mesh %s\n', meshpath);
        return;
    else
        delete(cachepath);
        clear mesh;
    end
end
    

fid = fopen(filename, 'r');

v = [];
f = [];

read_n = [];        % these are the normals we have read, indexed by order
n = [];             % same normals rewritten to be in vertex order

u = [];

has_normals = 0;
has_texture = 0;

C = textscan(fid, '%s %s %s %s', 1);
while ~feof(fid)
    if strcmp(C{1}{1}, 'v')
        v1 = str2double(C{2}{1});  v2 = str2double(C{3}{1});  v3 = str2double(C{4}{1});
        v = cat(1, v, [v1,v2,v3]);

    elseif strcmp(C{1}{1}, 'vn')
        n1 = str2double(C{2}{1});  n2 = str2double(C{3}{1});  n3 = str2double(C{4}{1});
        read_n = cat(1, read_n, [n1,n2,n3]);
        has_normals = 1;

    elseif strcmp(C{1}{1}, 'vt')
        u1 = str2double(C{2}{1});  u2 = str2double(C{3}{1});
        u = cat(1, u, [u1,u2]);
        has_texture = 1;
        
    elseif strcmp(C{1}{1}, 'f')
        fi1 = regexp(C{2}{1}, '\d+', 'match');
        fi2 = regexp(C{3}{1}, '\d+', 'match');
        fi3 = regexp(C{4}{1}, '\d+', 'match');
        
        v1 = str2num(fi1{1});
        v2 = str2num(fi2{1});
        v3 = str2num(fi3{1});
        
        % rewrite normals in vertex-index order
        if has_normals
            nindex = 2 + has_texture;
            n1 = str2num(fi1{nindex});
            n2 = str2num(fi2{nindex});
            n3 = str2num(fi3{nindex});
            n(v1,:) = read_n(n1,:);
            n(v2,:) = read_n(n2,:);
            n(v3,:) = read_n(n3,:);
        end
        
        % [TODO] rewrite texture coords in vertex-index order
        
        f = cat(1, f, [v1,v2,v3]);

    end
    C = textscan(fid, '%s %s %s %s', 1);
end


makemesh_options = [];
if ~isempty(findstr(options, 'B'))
    makemesh_options.noboundary = 1;
end
mesh = makeMesh(v,f,n,u, makemesh_options);

if ~isempty(findstr(options, 'n')) && numel(mesh.n) == 0
    mesh.n = estimateNormal(mesh);
end

if enable_mesh_cache
    save( cachepath, 'mesh', '-mat' );
end

fclose(fid);


end

