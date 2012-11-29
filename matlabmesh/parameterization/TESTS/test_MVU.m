% read a test mesh
mesh = readMesh('patch.obj');
mesh = readMesh('patch2.obj', 'n');

plotMesh(mesh, 'vefb');


% embed mesh vertices w/ MVU (Maximum Variance Unfolding)
%  [WARNING] this algorithm is very slow even for quite small datasets.
%     I recommend starting MaxIters at something small (10-20) and
%     increasing in small steps. The rightmost column in the output
%     is the error, once it stops changing significantly you can
%     probably stop. In the result output what you want to see
%     is the line 'Partial Success: SDP solved with reduced accuracy'.
%     That seems to be as good as it is going to get...   -RMS

X = mesh.v';
Knbrs = 5;         % larger values here make optimization much slower
MaxIters = 100;

pars.solver=0;  % CSDP
pars.factor = 0;  pars.slack = 3;  % only permit shrinking of distances
pars.maxiter = MaxIters; % max iterations

Dis=mvu_distance(X);  
[Ynbr,D]=mvu(Dis, Knbrs, pars); 
mesh.u = Ynbr(1:2, :)';


plotMesh(mesh, 'UefbO');



