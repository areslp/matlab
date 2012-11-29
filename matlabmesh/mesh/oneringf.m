function [ vTris ] = oneringf( mesh, nVertex )
% [vTris] = oneringf( mesh, nVertex )  
%   return indices of one-ring faces at nVertex

alltris = [ nonzeros(mesh.te(nVertex,:)), nonzeros(mesh.te(:,nVertex)) ];
vTris = unique(alltris);

end
