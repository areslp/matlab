% read a test mesh
mesh = readMesh('patch.obj');
mesh = readMesh('patch2.obj', 'n');

mesh = readMesh('bump_ref.obj');
mesh = clipEars(mesh);

plotMesh(mesh, 'vefb');


% set options.display = 1;  to see visual debug output
% set options.verbose = 0;  to get silent operation

% embed w/ Isomap.m using full distance matrix based on K-nbrs
clear options;
X = mesh.v';
Kiso = 8;
options.dims = 2;
options.display = 0; 
D = Isomap_L2_distance(X, X, 1); 
[Yiso,Riso,Eiso] = Isomap(D, 'k', Kiso, options);
Y = Yiso.coords{1};
mesh.u = Y';

plotMesh(mesh, 'UefbO');


% embed w/ IsomapII.m using full distance matrix based on K-nbrs
clear options;
X = mesh.v';
Kiso = 8;
options.dims = 2;
options.display = 0; 
D = Isomap_L2_distance(X, X, 1); 
[Yiso,Riso,Eiso] = IsomapII(D, 'k', Kiso, options);
Y = Yiso.coords{1};
mesh.u = Y';

plotMesh(mesh, 'UefbO');


% embed w/ IsomapII.m using distance function
clear options;
global Isomap_X;    % this global is required by Isomap_dfun (used below)
Isomap_X = mesh.v';
Kiso = 8;
options.dims = 2;
options.display = 0; 
[Yiso,Riso,Eiso] = IsomapII('Isomap_dfun', 'k', Kiso, options);
Y = Yiso.coords{1};
clear global Isomap_X;
mesh.u = Y';

plotMesh(mesh, 'UefbO');



% embed w/ IsomapII.m using distance function and sparse landmark points
clear options;
global Isomap_X;    % this global is required by Isomap_dfun (used below)
Isomap_X = mesh.v';
Kiso = 8;
options.dims = 2;
options.display = 0; 
nVerts = size(mesh.v,1);
nLandmarks = 10 * ceil(log(nVerts));
options.landmarks = ceil(rand(100,1)*nVerts);
[Yiso,Riso,Eiso] = IsomapII('Isomap_dfun', 'k', Kiso, options);
Y = Yiso.coords{1};
clear global Isomap_X;
mesh.u = Y';

plotMesh(mesh, 'UefbO');




% C-Isomap
% embed w/ IsomapII.m using full distance matrix based on K-nbrs
clear options;
X = mesh.v';
Kiso = 8;
options.dims = 2;
options.display = 0; 
options.conformal = 1;
D = Isomap_L2_distance(X, X, 1); 
[Yiso,Riso,Eiso] = IsomapII(D, 'k', Kiso, options);
Y = Yiso.coords{1};
mesh.u = Y';

plotMesh(mesh, 'UefbO');



