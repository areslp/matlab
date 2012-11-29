function [ W ] = meanvalueWeights( mesh, vertices )
% W = meanvalueWeights(mesh, vertices)
%   Compute mean-value weights 
%   - if vertices is undefined/empty, compute for entire mesh
%
%   [Ryan Schmidt  rms@dgp.toronto.edu  09/2009]


if ~ exist('vertices', 'var') || numel(vertices) == 0
    vertices = mesh.vidx;
end

n = numel(vertices);

W = mesh.e(vertices,:);
W(W~=0)=-1;

% [TODO] can we vectorize this? loop(if(if(loop))) ack!
% [TODO] conformal weights are symmetric...can we take advantage?

for i = 1:n
    qi = vertices(i);
    ov = oneringv(mesh,qi);
    for j = 1:numel(ov)
        qj = ov(j);
        
        % find opposing vertices (2 -> interior, 1 -> boundary)
        faces = [mesh.te(qi,qj), mesh.te(qj,qi)];
        faces = faces(faces~=0);
        verts = zeros(numel(faces),1);
        for k = 1:numel(faces)
            f = mesh.f(faces(k),:);
            verts(k) = f(f~=qi & f~=qj);
        end

        % sum up cotangents
        sumAB = 0;
        for k = 1:numel(verts)
            qk = verts(k);
            v1 = mesh.v(qk,:)-mesh.v(qi,:);
            v2 = mesh.v(qj,:)-mesh.v(qi,:);
            alphak = vangle(v1,v2);
            sumAB = sumAB + tan(alphak/2);
        end
        W(qi,qj) = sumAB / vmag2(mesh.v(qi,:)-mesh.v(qj,:));
    end
end



end
