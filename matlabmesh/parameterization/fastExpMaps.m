function [ vU, vV ] = fastExpMaps( mesh, nbrtype, nbrsize )
%UNTITLED1 Summary of this function goes here
%   Detailed explanation goes here

% write mesh to temp file
meshpath = [globalConfig('cachePath'), '_fastExpMaps_mesh_temp.obj'];
writeMesh(mesh, meshpath);

% generate expmaps
expmapCL = which('expmapCL.exe');
temppath = [globalConfig('cachePath'),'_fastExpMaps_data_temp.txt'];
if strcmp(nbrtype,'g') | strcmp(nbrtype,'k')
    command = [expmapCL,' ',meshpath,' ',nbrtype,' ',num2str(nbrsize),' > ',temppath];
elseif strcmp(nbrtype,'h')
    command = [expmapCL,' ',meshpath,' ',nbrtype,' ',num2str(nbrsize),' 6 > ',temppath];
else
    error(['[fastExpMaps] nbrtype ',nbrtype,' is not a valid option (only use g,k,h)']);
end
    
tic;
[status, result] = system(command);
emtime = toc;
if status < 0
    error(['expmapCL failed! msg was \n',result]);
else 
    fprintf('%s',result);
end

%read them back in
fprintf('[fastExpMaps] importing data...\n');
data = importdata(temppath);
fprintf('[fastExpMaps] parsing expmaps...\n');
N = numel(mesh.vidx);
i = data(:,1);
j = data(:,2);
vU = sparse(i,j, data(:,3),N,N);
vV = sparse(i,j, data(:,4),N,N);

end
