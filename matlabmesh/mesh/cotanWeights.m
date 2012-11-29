function [ W ] = cotanWeights( mesh, vertices, authalic, areaWeighted )
% W = cotanWeights(mesh, vertices, authalic, areaWeighted)
%   Compute cotangent weights 
%    vertices : vertex set (uses mesh.vidx if empty or undefined)
%    authalic : if = 1, compute authalic (area-preserving) weights instead [Desbrun02]
%    areaWeighted: if = 1, divide each cotan coeff by triangle area ( see [Mullen08] )
%   
%   - if authalic is nonzero, compute authalic cotan weights
%       instead of conformal cotan weights (see [Desbrun02])
%   
%   [Ryan Schmidt  rms@dgp.toronto.edu  07/2009]


if ~ exist('vertices', 'var') || numel(vertices) == 0
    vertices = mesh.vidx;
end
if ~exist('authalic', 'var')
    authalic = 0;
end
if ~exist('areaWeighted', 'var')
    areaWeighted = 0;
end
n = numel(vertices);

W = mesh.e(vertices,:);
W(W~=0)=-1;

% [TODO] can we vectorize this? loop(if(if(loop))) ack!
% [TODO] conformal weights are symmetric...can we take advantage?

if areaWeighted
    faceAreas = faceArea(mesh);
else
    faceAreas = ones(numel(mesh.fidx),1);
end

for i = 1:n
    qi = vertices(i);
    ov = oneringv(mesh,qi);
    for j = 1:numel(ov)
        qj = ov(j);
        
        % find opposing vertices (2 -> interior, 1 -> boundary)
        faces = [mesh.te(qi,qj), mesh.te(qj,qi)];
        faces = faces(faces~=0);
        verts = zeros(numel(faces),1);
        vertfaces = zeros(numel(faces),1);
        for k = 1:numel(faces)
            f = mesh.f(faces(k),:);
            verts(k) = f(f~=qi & f~=qj);
            vertfaces(k) = faces(k);
        end

        % sum up cotangents
        sumAB = 0;
        if authalic
            for k = 1:numel(verts)
                qo = verts(k);
                v1 = mesh.v(qi,:)-mesh.v(qj,:);
                v2 = mesh.v(qo,:)-mesh.v(qj,:);
                sumAB = sumAB + vcot(v1,v2);
            end               
            sumAB = sumAB / vmag2( mesh.v(qi,:) - mesh.v(qj,:) );
        else
            for k = 1:numel(verts)
                qo = verts(k);
                v1 = mesh.v(qi,:)-mesh.v(qo,:);
                v2 = mesh.v(qj,:)-mesh.v(qo,:);
                sumAB = sumAB + vcot(v1,v2) / faceAreas(vertfaces(k));
            end
        end
        W(qi,qj) = sumAB;
        
    end
end



end
