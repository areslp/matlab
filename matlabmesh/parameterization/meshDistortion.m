function [ Dangle, Darea, L2, Linf, Cangle, Carea ] = meshDistortion( mesh, ignoreBoundary, perVertex )
%MESHDISTORTION compute distortion over mesh
%  [TODO] check that these are actually working. Only L2 and Linf
%    give values of 1 if we check distortion wrt same planar mesh...

if ~exist('ignoreBoundary','var')
    ignoreBoundary = 0;
end
if ~exist('perVertex', 'var')
    perVertex = 0;
end

use_faces = mesh.f;
if ignoreBoundary
    face_roi = mesh.fidx( mesh.isboundaryf == 0 );
    use_faces = mesh.f( face_roi, :);
end


% construct scaled u that has same area as mesh
areaV = meshArea(mesh.v, use_faces);
areaU = meshArea(mesh.u, use_faces);
centroid = sum(mesh.u) / size(mesh.u,1);
scale = sqrt(areaV / areaU);
scaled_u = vadd(mesh.u, -centroid) * scale;

[Dangle,Darea] = triDistortion(mesh.v, scaled_u, use_faces);
[face_L2,face_Linf,TareaV,TareaU] = triStretch(mesh.v, scaled_u, use_faces);

sumAreaV = sum(TareaV);
sumAreaU = sum(TareaU);
areaScale = sumAreaU / sumAreaV;

Dangle = sum( Dangle .* TareaV ) / sumAreaV;
Darea = sum( Darea .* TareaV ) / sumAreaV;
face_L2 = face_L2 .* face_L2 .* TareaV;

if perVertex
    face_L2 = sqrt( face_L2 / sumAreaV ) * sqrt(areaScale);
    face_Linf = face_Linf * sqrt(areaScale);
    L2 = zeros(size(mesh.vidx));
    for i = 1:numel(mesh.vidx)
        fi = oneringf(mesh,i);
        L2(i) = sum(face_L2(fi)) / numel(fi);
        Linf(i) = max(face_Linf(fi));
    end
else
    L2 = sqrt( sum(face_L2) / sumAreaV ) * sqrt(areaScale);
    Linf = max(face_Linf) * sqrt(areaScale);
end



% compute Dirichlet and Chi Energy for parameterization from [Desbrun02]
%  (does not include boundary triangles)

weightsAngle = cotanWeights(mesh);
[i,j] = find(weightsAngle);
vdists = weightsAngle;
vdists(sub2ind(size(weightsAngle),i,j)) = vmag2(mesh.u(i,:) - mesh.u(j,:));
%vdists( mesh.vidx(mesh.isboundaryv==1), : ) = 0;
if perVertex
    Cangle = full( sum(weightsAngle .* vdists) )';
else
    Cangle = full( sum(sum(weightsAngle .* vdists)) );
end

weightsArea = cotanWeights(mesh, [], 1);
if perVertex
    Carea = full( sum(weightsArea .* vdists) )';
else
    Carea = full( sum(sum(weightsArea .* vdists)) );
end



end
