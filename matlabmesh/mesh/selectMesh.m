function [ idx ] = selectMesh( mesh, f, type )
%SELECTMESH select components of mesh - returns indices of selected elements
%   f = implicit function that takes arguments x,y,z - return value < 0 is inside
%   type = selection type (either 'v' or 'f')

numV = numel(mesh.vidx);
signs = zeros(numV,1);
for i = 1:numV
    signs(i) = f( mesh.v(i,1), mesh.v(i,2), mesh.v(i,3) );
end

idx = [];
if strcmp(type, 'f')
    for i = 1:numel(mesh.fidx)
        f = mesh.f( i, : );
        inside = 1;
        for j = 1:numel(f)
            inside = inside && signs(f(j)) < 0;
        end
        if inside
            idx = [idx;i];
        end
    end            
      
else
    idx = mesh.vidx( find(signs < 0) );
end


end
