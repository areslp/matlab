function [ N ] = faceNormal( mesh, faces )
% [ N ] = faceNormal( mesh, faces )
%  compute normals of mesh faces
%  If faces empty or undefined, compute areas for entire mesh 
%  !!! currently assumes faces are triangles !!!

if ~ exist('faces', 'var') || numel(faces) == 0
    faces = mesh.fidx;
end
n = numel(faces);

N = zeros(n,3);
for i = 1:n
    f = mesh.f(i,:);
    fv = mesh.v(f,:);
    N(i,:) = trinormal(fv(1,:), fv(2,:), fv(3,:));
end

end