function [ submesh, vertmap, facemap ] = subMesh( mesh, face_roi )
%[ submesh, vertmap, facemap ] = subMesh( mesh, face_roi )
%   pick out submesh using face ROI

newverts = [];


newf = mesh.f(face_roi,:);
vert_roi = unique( [newf(:,1); newf(:,2); newf(:,3)] );
newv = mesh.v(vert_roi,:);
newn = [];
if ( ~ isempty(mesh.n) )
    newn = mesh.n(vert_roi,:);
end
newu = [];
if ( ~ isempty(mesh.u) )
    newu = mesh.u(vert_roi,:);
end

vertmap = [1:numel(vert_roi); vert_roi']';
facemap = [1:numel(face_roi); face_roi']';

% rewrite face indices
for i = 1:numel(face_roi)
    f = newf(i,:);
    for j = 1:numel(f);
        f(j) = vertmap(find(vertmap(:,2) == f(j)), 1);
    end
    newf(i,:) = f;
end

submesh = makeMesh(newv, newf, newn, newu);

end
