function [ pointset ] = readPoints( filename, options )
% [ pointset ] = readPoints( filename, options )
%  Currently only supports .OBJ-style vertex lists
%  options specifies flags for reading
%   'n' - generate normals if not contained in file
%
%   [Ryan Schmidt  rms@dgp.toronto.edu  09/2009]

if ~exist('options','var')
    options = '';
end

enable_cache = globalConfig('enableMeshCaching');

% check that filename exists
filepath = which(filename);
if ~exist('filepath')
    error(['File ',filename,' does not exist in matlab path\n']);
end

% use cached version if possible (opens much faster)
cachepath = [filepath, '.ptscache'];
if enable_cache & exist(cachepath)
    S = load(cachepath, '-mat');
    pointset = S.pointset;
    if pointset.version == globalConfig('pointsetVersion');
        fprintf('[readMesh] loaded cached data %s\n', filepath);
        return;
    else
        delete(cachepath);
        clear pointset;
    end
end
    

fid = fopen(filename, 'r');

v = [];
n = [];
u = [];
c = [];

has_normals = 0;
has_texture = 0;

C = textscan(fid, '%s %s %s %s', 1);
while ~feof(fid)
    if strcmp(C{1}{1}, 'v')
        v1 = str2double(C{2}{1});  v2 = str2double(C{3}{1});  v3 = str2double(C{4}{1});
        v = cat(1, v, [v1,v2,v3]);

    elseif strcmp(C{1}{1}, 'vn')
        n1 = str2double(C{2}{1});  n2 = str2double(C{3}{1});  n3 = str2double(C{4}{1});
        n = cat(1, n, [n1,n2,n3]);

    elseif strcmp(C{1}{1}, 'vt')
        u1 = str2double(C{2}{1});  u2 = str2double(C{3}{1});
        u = cat(1, u, [u1,u2]);
        
    elseif strcmp(C{1}{1}, 'vc')
        c1 = str2double(C{2}{1});  c2 = str2double(C{3}{1});  c3 = str2double(C{4}{1});
        c = cat(1, c, [c1,c2,c3]);
        
    end
    C = textscan(fid, '%s %s %s %s', 1);
end

if size(n,1) ~= size(v,1)
    n = [];
end
if size(u,1) ~= size(v,1)
    u = [];
end
if size(c,1) ~= size(v,1)
    c = [];
end

pointset = makePointSet(v,n,u,c);

% extract largest connected component
G = (pointset.e>0);
blocks = components(G)';
count = zeros(1, max(blocks));
for i=1:max(blocks)
    count(i) = length(find(blocks == i));
end
[count, block_no] = max(count);
conn_comp = find(blocks == block_no);
pointset = makePointSet( pointset.v(conn_comp,:), pointset.n(conn_comp,:) );


%if ~isempty(findstr(options, 'n')) && numel(pointset.n) == 0
%    pointset.n = estimateNormal(pointset);
%end

if enable_cache
    save( cachepath, 'pointset', '-mat' );
end

fclose(fid);

end

