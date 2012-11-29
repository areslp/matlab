%% halfedge tests

mesh = readMesh('patch.obj');


% test that he->vertex mapping and flip edges are right
[ii,jj] = find(mesh.e~=0);
for ei = 1:size(i,1);
    i = ii(ei);  j = jj(ei);
    
    l1 = mesh.he(i,j);
    l2 = mesh.he(j,i);
    
    if (l1 ~= 0 && mesh.hev(l1) ~= i) || (l2 ~= 0 && mesh.hev(l2) ~= j)
        fprintf('hev is broken? ei is %d\n', ei);
    end
    
    l1f = mesh.heflip(l1);
    l2f = mesh.heflip(l2);
    
    if l1f ~= l2 || l2f ~= l1
        fprintf('heflip is broken? ei is %d\n', ei);
    end
end

% check that vertex->he mapping is working 
for vi = 1:size(mesh.vidx)
    hestart = mesh.veh(vi);
    henext = mesh.henext(hestart);
    henext2 = mesh.henext(henext);
    henext3 = mesh.henext(henext2);
    if henext3 ~= hestart
%        fprintf('face %d %d %d (%d)\n', mesh.hev(hestart), mesh.hev(henext), mesh.hev(henext2), mesh.hev(henext3) );
        error('veh or henext is broken? vi %d\n', vi);
    end
end
    
% check that face->he mapping is working
for fi = 1:size(mesh.fidx);
    f = mesh.f(fi,:);
    he1 = mesh.feh(fi);
    he2 = mesh.henext(he1);
    he3 = mesh.henext(he2);
    f2 = [mesh.hev(he1), mesh.hev(he2), mesh.hev(he3)];
    
    % this doesn't actually guarantee anything, but it should work most of the time
    if numel(unique([f,f2])) ~= 3
%        fprintf('face %d - %d %d %d ||| %d %d %d\n',fi, f(1), f(2), f(3), mesh.hev(he1), mesh.hev(he2), mesh.hev(he3));
        error('feh or henext is broken? fi %d\n', fi);
    end
end


% check that ordered onering iterations are working
for vi = 1:size(mesh.vidx)
    if mesh.isboundaryv(vi)
        continue;
    end
    
    hecur = mesh.veh(vi);
    onering_he = mesh.hev( mesh.heflip(hecur) );
    hecur = mesh.henext( mesh.heflip( hecur ) );
    curv = mesh.hev( mesh.heflip(hecur) );
    while curv ~= onering_he(1);
        onering_he = [onering_he, curv];
        hecur = mesh.henext( mesh.heflip( hecur ) );
        curv = mesh.hev( mesh.heflip(hecur) );        
    end

    [onering_t, onering_v] = onering(mesh,vi);
    
    if numel(unique([onering_he,onering_v'])) ~= numel(onering_v)
        [onering_he;onering_v']
        error('onering check failed - vi %d\n', vi);
    end   
end

