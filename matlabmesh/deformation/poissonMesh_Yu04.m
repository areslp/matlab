function [ output_mesh ] = poissonMesh_Yu04( mesh, boundary_cons, target_tris )
%POISSONMESH_YU04 Summary of this function goes here
%   Detailed explanation goes here
%   boundary_cons:      rows of [vtx_i, x, y, z]
%
% [RMS TODO] replace target_tris with a list of transformation
%    matrices for original tris. Move (broken) transformation estimation
%    to another file. Also add in triangle area weights (from paper)

output_mesh = mesh;
nVerts = numel(mesh.vidx);
orig_tris = [ mesh.v(mesh.f(:,1),:), mesh.v(mesh.f(:,2),:), mesh.v(mesh.f(:,3),:) ];

A = mesh.e(mesh.vidx,:);
A(A~=0) = 0;

% construct A matrix
for i = 1:nVerts
    vi = mesh.vidx(i);
    
    if mesh.isboundaryv(vi) == 1
        A(vi,vi) = 1;
        continue;
    end
    
    vTris = onering(mesh, vi);
    for oi = 1:numel(vTris)
        k = vTris(oi);
        f = mesh.f(k,:);
        Tk = orig_tris(k,:);
        
        gradBik = gradB(f, vi, Tk);
        
        for j = 1:3
            gradBjk = gradB(f, f(j), Tk);
            A(vi,f(j)) = A(vi,f(j)) + vdot(gradBik, gradBjk);
        end
    end
end


% construct B matrix based on W field. W field is defined
% by applying transformations to gradients *and vertices* of
% original triangles
B = zeros(nVerts,3);
for i = 1:nVerts
    vi = mesh.vidx(i);
    vTris = onering(mesh, vi);

    if mesh.isboundaryv(vi) == 1
        B(vi,:) = mesh.v(vi,:);
        continue;
    end
    
    for oi = 1:numel(vTris)
        k = vTris(oi);
        f = mesh.f(k,:);
        Tk = orig_tris(k,:);

        Tknew = target_tris(k,:);
        mat = estimateTransform(Tk, Tknew);
        
        % transform Tk to Tknew centroid/plane (based on estimated rotation)
        v1 = Tk(1:3);  v2 = Tk(4:6);  v3 = Tk(7:9);
        c1 = (v1 + v2 + v3) / 3;
        c2 = (Tknew(1:3) + Tknew(4:6) + Tknew(7:9)) / 3;
        v1 = (mat * (v1 - c1)')' + c2;
        v2 = (mat * (v2 - c1)')' + c2;
        v3 = (mat * (v3 - c1)')' + c2;
        Tknew = [v1,v2,v3];      
        
        gradBik = gradB(f, vi, Tk);
                
        x = Tknew([1,4,7]);  y = Tknew([2,5,8]);   z = Tknew([3,6,9]);
        for j = 1:3
             gradBjk = (mat * gradB(f, f(j), Tk)')';
             B(vi,1) = B(vi,1) + x(j) * vdot( gradBik, gradBjk );
             B(vi,2) = B(vi,2) + y(j) * vdot( gradBik, gradBjk );
             B(vi,3) = B(vi,3) + z(j) * vdot( gradBik, gradBjk );
        end
    end
    
end

% set boundary constraints
nBoundary = size(boundary_cons,1);
for i = 1:nBoundary;
    bi = boundary_cons(i,1);
    bp = boundary_cons(i,2:4);
    B(bi,:) = bp;
end

output_mesh.v(:,1) = A \ B(:,1);
output_mesh.v(:,2) = A \ B(:,2);
output_mesh.v(:,3) = A \ B(:,3);


% this plots all the transformed triangles (but is very slow!)
% newplot;
% hold all
% for ti = 1:numel(mesh.fidx)
%     Tk = orig_tris(ti,:);
%     Tknew = target_tris(ti,:); 
%     mat = estimateTransform(Tk, Tknew);    
%     v1 = Tk(1:3);  v2 = Tk(4:6);  v3 = Tk(7:9);
%     c1 = (v1 + v2 + v3) / 3;
%     c2 = (Tknew(1:3) + Tknew(4:6) + Tknew(7:9)) / 3;
%     v1 = (mat * (v1 - c1)')' + c2;
%     v2 = (mat * (v2 - c1)')' + c2;
%     v3 = (mat * (v3 - c1)')' + c2;
%     T = [v1,v2,v3];       
%     patch(T([1,4,7]),T([2,5,8]),T([3,6,9]),'red')    
% end
% hold off;
% drawnow;




function [a] = area(tri)
v1 = tri(1:3);  v2 = tri(4:6);  v3 = tri(7:9);
a = triarea(v1,v2,v3);



% estimate transformation that has been applied to a triangle
% [RMS] this does not work very well..
function [mat] = estimateTransform(tri, dest)

% construct coordinate frame of dest tri
d1 = dest(1:3);  d2 = dest(4:6);  d3 = dest(7:9);
n2 = ncross(d2-d1,d3-d1);
e1 = normalize(d2-d1);
e2 = ncross(e1,n2);
m2 = [e1;e2;n2]';

% construct coordinate frame of original tri
o1 = tri(1:3);  o2 = tri(4:6);  o3 = tri(7:9);
n1 = ncross(o2-o1,o3-o1);
e1 = normalize(o2-o1);
e2 = ncross(e1,n1);
m1 = [e1;e2;n1]';

% this is matrix that takes original tri into plane of dest tri
mat = m2 * m1';

% subtract off centroids
c = (d1 + d2 + d3)/3;
d1 = normalize(d1 - c);  d2 = normalize(d2 - c);  d3 = normalize(d3 - c);
c = (o1 + o2 + o3)/3;
o1 = normalize(o1 - c);  o2 = normalize(o2 - c);  o3 = normalize(o3 - c);

% transform original vectors into new frame
o1 = (mat * o1')';  o2 = (mat * o2')';  o3 = (mat * o3')';

% compute average in-frame rotation angle for each vertex
t1 = angle(o1,d1);  t2 = angle(o2,d2);  t3 = angle(o3,d3);
r = -(t1+t2+t3)/3;
mat = axisrot(n2, r) * m2 * m1';









% gradient of triangle linear basis function for vertex i
% tri = [ i1, i2, i3 ]  (indices of triangle verts)
% Tk = row of vertices of triangle   [ v1 , v2 , v3 ]
function [ d ] = gradB(tri, i, Tk)

% get position of vertex i
vi = find(tri == i);
v1 = Tk((1:3)*vi);

% get positions of opposing verts of tri
otheri = find(tri ~= i);
t1 = otheri(1);  t2 = otheri(2);
v2 = Tk((1:3)*t1);
v3 = Tk((1:3)*t2);

% gradient direction is perp to opposite edge, lying in face
n = ncross(v2 - v1, v3 - v1);
e = v3 - v2;
d = ncross(e,n);   % this is gradient direction

% make sure direction points towards vertex i
midpoint = 0.5 * ( v2 + v3 );
if vdot(d, v1-midpoint) < 0
    d = -d;
end

% scale vector so that B(v1) = 1
l = vdot(v1-v2, d);
d = d / l;

% scale by area (?)
%A = triarea(v1,v2,v3);
%v = v / (2*A);



