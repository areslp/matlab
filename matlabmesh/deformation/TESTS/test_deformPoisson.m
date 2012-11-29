

mesh = readMesh('plane.obj');

% construct boundary constraints 
% have to constrain all boundary points (Dirichlet condition)
vboundary = mesh.vidx(mesh.isboundaryv~=0);
consB = [vboundary, mesh.v(vboundary,:), ones(numel(vboundary),1)*10];

% construct handle constraints
% (shifting center vertex up in Y)
center = sum(mesh.v) / numel(mesh.vidx);
[val,idx] = min( vmag2(vadd(mesh.v,-center)) );
consH = [idx, mesh.v(idx,:)+[0,0.13,0], 10];

% construct rotation for each mesh triangle
% (here we are rotating outwards, to form a bubble)
rotations = [];
center = sum(mesh.v) / numel(mesh.vidx);
maxdist = max(vmag(vadd(mesh.v,-center)));
rotvals = [];
rotn = [];
for ti = 1:numel(mesh.fidx)
    face = mesh.f(ti,:);
    tricenter = sum(mesh.v(face,:)) / 3;
    trivec = tricenter - center;
    [trivec,dist] = normalize(trivec);
    rotaxis = [-trivec(3),0,trivec(1)];
    rotmat = axisrot(rotaxis, -(pi/2) * (dist/maxdist) );
    rotations = [rotations; reshape(rotmat,1,9) ];
    rotvals = [rotvals;(dist/maxdist)];
end

% compute and plot deformation
deformed_mesh = deformPoisson(mesh,rotations, [consB;consH] );
figure; hold all;
plotMesh(deformed_mesh,'fen');
scatter3( consB(:,2), consB(:,3), consB(:,4), 'rs', 'filled' );
if ~isempty(consH)
    scatter3( consH(:,2), consH(:,3), consH(:,4), 'bo', 'filled' );
end
hold off;




