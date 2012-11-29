

mesh = readMesh('plane.obj');
vmin = find(mesh.v(:,1) == min(mesh.v(:,1)));
vmax = find(mesh.v(:,1) == max(mesh.v(:,1)));
vpos = [vmin;vmax];
consP = [vpos, mesh.v(vpos,:), ones(numel(vpos),1)*5];

weightF = 5;
matfix = eye(3,3);
matrot = axisrot( [0,0,1], -pi*0.5 );
consF_fix = [vmin, repmat(reshape(matfix',1,9),numel(vmin),1), ones(numel(vmin),1)*weightF];
consF_rot = [vmax, repmat(reshape(matrot',1,9),numel(vmax),1), ones(numel(vmax),1)*weightF];
consF = [consF_fix;consF_rot];


deformed_mesh = deformRotInvCoords(mesh,consF,consP);
%deformed_mesh.v = mesh.v;
figure; hold all;
plotMesh(deformed_mesh,'fen');
scatter3( consP(:,2), consP(:,3), consP(:,4), 'rs', 'filled' );
scatter3( deformed_mesh.v(consF(:,1),1), deformed_mesh.v(consF(:,1),2), deformed_mesh.v(consF(:,1),3), 'bo', 'filled' );
hold off;








mesh = readMesh('armadilloman_tiny.obj');
mesh.n = estimateNormal(mesh, [], 'faceavg_area');
nose = 73;

%vfix = [4,5,11,12,  16,17,18,24,  61,75 ]';
vfix = [4,5,11,12,  16,17,18,24 ]';
consP_fix = [vfix, mesh.v(vfix,:), ones(numel(vfix),1)*10];
consP_move = [];
consP_move = [nose, mesh.v(nose,:)+[0,-90,30], 3];
consP = [consP_fix;consP_move];

weightF = 5;
matfix = eye(3,3);
consF_fix = [vfix, repmat(reshape(matfix',1,9),numel(vfix),1), ones(numel(vfix),1)*weightF];
noserot = axisrot( [1,0,0], pi*0.5 );
consF_nose = [];
consF_nose = [nose, reshape(noserot',1,9), weightF];
consF = [consF_fix;consF_nose];

deformed_mesh = deformRotInvCoords(mesh,consF,consP);
%deformed_mesh.v = mesh.v;
figure; hold all;
plotMesh(deformed_mesh,'fen');
plotMesh(mesh,'en');
scatter3( consP(:,2), consP(:,3), consP(:,4), 'rs', 'filled' );
scatter3( deformed_mesh.v(consF(:,1),1), deformed_mesh.v(consF(:,1),2), deformed_mesh.v(consF(:,1),3), 'bo', 'filled' );
hold off;



















mesh = readMesh('plane.obj');

% construct boundary constraints 
% have to constrain all boundary points (Dirichlet condition)
vboundary = mesh.vidx(mesh.isboundaryv~=0);
consB = [vboundary, mesh.v(vboundary,:), ones(numel(vboundary),1)*10];

% construct rotation for each boundary vertex
% (here we are rotating outwards, to form a bubble)
center = sum(mesh.v) / numel(mesh.vidx);
maxdist = max(vmag(vadd(mesh.v,-center)));
rotn = [];
numB = numel(vboundary);
consF = [ vboundary, zeros(numB,9), ones(numB,1)*10 ];
for ii = 1:numB
    i = consB(ii,1);
    vec = mesh.v(i,:) - center;
    [vec,dist] = normalize(vec);
    rotaxis = [-vec(3),0,vec(1)];
    rotmat = axisrot(rotaxis, -(pi/2) * (dist/maxdist) );
    consF(ii,2:10) = reshape(rotmat',1,9); 
end

% compute and plot deformation
[rinv_mesh,rotations] = deformRotInvCoords(mesh,consF,consB);
figure; hold all;
plotMesh(rinv_mesh,'fen');
scatter3( consB(:,2), consB(:,3), consB(:,4), 'rs', 'filled' );
hold off;

% compute and plot deformation
vboundary = mesh.vidx(mesh.isboundaryv~=0);
consB = [vboundary, mesh.v(vboundary,:), ones(numel(vboundary),1)*1000];
consH = [];

deformed_mesh = deformLaplacian(mesh, rotations, [consB;consH] );
figure; hold all;
plotMesh(deformed_mesh,'fen');
scatter3( consB(:,2), consB(:,3), consB(:,4), 'rs', 'filled' );
if ~isempty(consH)
    scatter3( consH(:,2), consH(:,3), consH(:,4), 'bo', 'filled' );
end
hold off;

