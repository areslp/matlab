function [ A ] = faceArea( mesh, faces )
% [ A ] = faceArea( mesh, faces )
%  compute area of mesh faces
%  If faces empty or undefined, compute areas for entire mesh 
%  !!! currently assumes faces are triangles !!!

if ~ exist('faces', 'var') || numel(faces) == 0
    faces = mesh.fidx;
end
n = numel(faces);

A = zeros(n,1);
for i = 1:n
    f = mesh.f(i,:);
    fv = mesh.v(f,:);
    A(i) = triarea(fv(1,:), fv(2,:), fv(3,:));
end

end

