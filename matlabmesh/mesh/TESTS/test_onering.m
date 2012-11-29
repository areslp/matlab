mesh = readMesh('sphere.obj');

nVert = 5;
[t,v] = onering(mesh, nVert, 'ccw');

dangle = (2*pi) / numel(v);

[ringmesh,vmap,tmap] = subMesh(mesh,t');
ringmesh.v(vmap(vmap(:,2)==nVert),:) = [0,0,0];
for i = 1:numel(v)
    vi = vmap(vmap(:,2)==v(i));
    ringmesh.v(vi,:) = [cos(i*dangle),sin(i*dangle),0];
end
plotMesh(ringmesh,'fevi');



vcenter = mesh.v(nVert,:);
e1 = vadd(mesh.v(v,:),-vcenter);
e2 = vadd(mesh.v(circshift(v,[0,-1]),:),-vcenter);
angles = vangle(e1,e2);
anglediff = (2*pi-sum(angles))/numel(v);
vangles = cumsum(angles+anglediff);
vangles = vangles-vangles(1);
vr = vmag(e1);
newv = [vr.*cos(vangles),vr.*sin(vangles),zeros(numel(v),1)];

[ringmesh,vmap,tmap] = subMesh(mesh,t');
figure;plotMesh(ringmesh);
ringmesh.v(vmap(vmap(:,2)==nVert),:) = [0,0,0];
for i = 1:numel(v)
    vi = vmap(vmap(:,2)==v(i));
    ringmesh.v(vi,:) = newv(i,:);
end
figure;plotMesh(ringmesh);

% validate oneringExpMap function, which effectively tests
% that onering() returns correctly-sorted nbrhood
% [TODO] checks often fail at boundary vertices...
for vi = 1:numel(mesh.vidx)
    [ringv,posv, ringt] = oneringExpMap(mesh,vi);
    
    % check that vertex ordering was preserved in first triangle
    [v2,v3] = tripick(mesh.f(ringt(1),:),vi);
    if ringv(1) ~= v2 | ringv(2) ~= v3
        fprintf('vertex ordering changed at vertex %d??\n', vi);
    end
    
    % check that flattened areas are positive
    %  (doesn't actually check last vertex...)
    for i = 1:numel(ringv)-1
        areai = triarea2s([0,0,0],posv(i,:),posv(i+1,:));
        if areai < 0
            fprintf('vertex %d has negative area in flattening');
        end
    end
    
    
    tmp = mesh;
    tmp.v(vi,:) = [0,0,0];
    tmp.v(ringv,:) = posv;

    % check that angle sum around sorted ringv is 2pi
    vcenter = tmp.v(vi,:);
    e1 = vadd(tmp.v(ringv,:),-vcenter);
    e2 = vadd(tmp.v(circshift(ringv,[0,-1]),:),-vcenter);
    angles = vangle(e1,e2);    
    if abs(2*pi-sum(angles)) > sqrt(eps)
        fprintf('vert %d: angle sum is %f\n', vi, sum(angles));
    end
    
    % check that angle sum around triangles is 2pi
    anglesum = 0;
    tris = oneringf(mesh,vi);
    for ti = 1:numel(tris)
        tri = mesh.f(tris(ti),:);
        [n1,n2] = tripick(tri,vi);
        anglesum = anglesum + vangle(tmp.v(n1,:)-tmp.v(vi,:), tmp.v(n2,:)-tmp.v(vi,:));
    end
    if abs(2*pi-anglesum) > sqrt(eps)
        fprintf('vert %d: angle sum is %f\n', vi, 2*pi-anglesum);
    end        
end
    