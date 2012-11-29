function [ vVertices ] = oneringv( mesh, nVertex )
% [ vVertices ] = oneringv( mesh, nVertex )
%    return indices of one-ring vertices at nVertex

vVertices = find(mesh.e(nVertex,:)~=0)';

end
