function [ vTris, vVerts ] = onering( mesh, nVert, mode )
%ONERING returns list of tris and verts in one-ring of vertex
%   mode == 'ccw' : sort one-ring 

vVerts = find(mesh.e(nVert,:)~=0)';

vTris = cat(1, full(mesh.te(nVert,vVerts))', full(mesh.te(vVerts,nVert)));
vTris = unique(vTris);
vTris(find(vTris==0)) = [];

if ~exist('mode','var')
    return;
end
if strcmp(mode,'ccw') && ~isempty(vVerts)
    % need to set vVerts(1) to 'first' boundary vert.
    % [TODO] unfortunately this code doesn't do that...only
    %  sets to *a* boundary vert...
    %  !!! also fails for 'strip' parts of meshes, which have more than 2
    %  bdry nbrs !!!
    if mesh.isboundaryv(nVert)
        isb = mesh.isboundaryv(vVerts);
        swapi = find(isb~=0);
        tmp = vVerts(1); vVerts(1) = vVerts(swapi(1));  vVerts(swapi(1)) = tmp;
    end
    
    curv = vVerts(1);
    vSorted = [curv];  
    rest = vVerts(2:end);
    tnbrs = [mesh.te(nVert,curv), mesh.te(curv,nVert)];
    vnbrs = [0,0];
    for j = 1:2
        if tnbrs(j) ~= 0 
            vnbrs(j) = tripick2(mesh.f(tnbrs(j),:),nVert,curv);
        end
    end
%    vnbrs = [tripick2(mesh.f(tnbrs(1),:),nVert,curv), tripick2(mesh.f(tnbrs(2),:),nVert,curv) ];
    prev = curv;
    usev = 1;  if vnbrs(usev) == 0; usev = 2; end;
    curv = vnbrs(usev);
    tSorted = [tnbrs(usev)];
    
    while ~isempty(rest)
        vSorted = [vSorted,curv];
        rest = rest(rest~=curv);
        
        tnbrs = [mesh.te(nVert,curv), mesh.te(curv,nVert)];
        if tnbrs(1) == 0 | tnbrs(2) == 0
            break;  % hit boundary - need to stop...
        end  
        vnbrs = [0,0];
        for j = 1:2
            if tnbrs(j) ~= 0 
                vnbrs(j) = tripick2(mesh.f(tnbrs(j),:),nVert,curv);
            end
        end        
        %vnbrs = [tripick2(mesh.f(tnbrs(1),:),nVert,curv), tripick2(mesh.f(tnbrs(2),:),nVert,curv) ];
        if vnbrs(1) == prev
            prev = curv;
            curv = vnbrs(2);
            tSorted = [tSorted, tnbrs(2)];
        elseif vnbrs(2) == prev
            prev = curv;
            curv = vnbrs(1);
            tSorted = [tSorted, tnbrs(1)];
        else
            error('can we get here?');
        end
    end
    
    if numel(tSorted)~=numel(vTris) || numel(vSorted)~=numel(vVerts)
        warning('onering(): could not sort!');
        return;
    end
    vTris = tSorted;
    vVerts = vSorted;
    
    
end
