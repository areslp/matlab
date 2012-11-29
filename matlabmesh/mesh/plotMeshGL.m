function [ ] = plotMeshGL( mesh, mode, vC )
% plotMeshGL( mesh, mode, vC )
%   plot a mesh using external OpenGL viewer
%     (mode flags currently not implemented)
%
%
%   [Ryan Schmidt  rms@dgp.toronto.edu  09/2009]


% write mesh to temp file
meshpath = [globalConfig('cachePath'), '_plotMeshGL_mesh_temp.obj'];
writeMesh(mesh, meshpath);

% generate expmaps
view3D = which('view3D.exe');
command = [view3D, ' ', meshpath];

dos([command, ' &']);

