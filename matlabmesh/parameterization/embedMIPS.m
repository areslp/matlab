function [ uvmesh ] = embedMIPS( mesh, max_iter )
%EMBEDMIPS computes free-boundary MIPS parameterization of mesh
%   [Ryan Schmidt  rms@dgp.toronto.edu  09/2008]
%   Described in 'MIPS: An Efficient Global Parameterization Method' 
%   by Hormann & Griener)

% embed mesh inside circle as initialization
mymesh = mesh;
n = size(mymesh.n, 1);
boundaryUV = embedBoundary( mesh, 'circle' );
weights = makeOneRingWeights(mesh, 'uniform');
uvmesh = mesh;
uvmesh.u = embedInterior(mesh, boundaryUV, weights);
uvmesh.v = [uvmesh.u(:,1), uvmesh.u(:,2), zeros(n,1)];
uvmesh = orientMesh(uvmesh, 'ccw2');

% make sure that face orientations are the same,
%  and then compute interior triangle angles
mymesh.f = uvmesh.f;
tAngles = getAnglesV(mymesh);

uv = uvmesh.u;
uv_lastCheck = uv;
EPrev = MIPSEnergy( uv, uvmesh.f, tAngles );

% iterative MIPS optimization as described in 
% 'Using Most Isometric Parameterizations for Remeshing Polygonal
% Surfaces', Labsik, Hormann, & Greiner
if ~exist('max_iter'), max_iter = 100; end;
tol = 0.0001;
for ni = 1:max_iter
    
    % do tolerance check every K iterations. Also, 
    % sometimes we get a triangle flip (the optimizer
    % allows this because it has negative energy, which is 'good')
    % When that happens, reset uv to last known good values
    % (randomized pick should prevent problem from repeating)
    if mod(ni,10) == 0
        ENew = MIPSEnergy( uv, uvmesh.f, tAngles );
        fprintf('iter %6d : MIPSEnergy %8.4f  (delta %8.4f)\n', ni, ENew, ENew-EPrev);
        if ( ENew < 0 )
            uv = uv_lastCheck;
        elseif ( EPrev - ENew < tol )
            break;
        else
            uv_lastCheck = uv;
        end
    end
   
   %pick random vert to optimize
   vi = 1 + round(rand(1,1) * (n-1));
   vTris = onering( uvmesh, vi );
   nTris = size(vTris,1);
   
   % re-order face vertices so that first vtx is vi
   faces = uvmesh.f(vTris,:);
   for ti = 1:nTris
       t = faces(ti,:);
       i = find(t==vi);
       if ( i == 2 ) 
           t = [ t(2), t(3), t(1) ];
       elseif ( i == 3 )
           t = [ t(3), t(1), t(2) ];
       end
       faces(ti,:) = t;
   end
   
   % optimize this vertex wrt it's fixed one-ring
%    options = optimset('Display','iter', 'FunValCheck', 'on', 'GradObj', 'on');
    options = optimset('Display','notify', 'FunValCheck', 'on', 'GradObj', 'on');
    [minuv, minE] = fminunc( @(uv_i) MIPSEnergyVGrad_vec(uv_i, vi, mymesh, uv, faces) , uv(vi,:)', options );
    uv(vi,:) = minuv;

%   {'after: ', MIPSEnergyV(mymesh, uv, faces), MIPSEnergy(uv, faces, tAngles) }
 
   
end

uvmesh.u = uv;
uvmesh.v = [uvmesh.u(:,1), uvmesh.u(:,2), zeros(n,1)];




% compute MIPS energy for triangles in faces
function [ E ] = MIPSEnergyV( mesh3D, uv, faces )
E = 0;
nTris = size(faces,1);
for ti = 1:nTris   
    t = faces(ti,:);
    vA = uv(t(1),:);  vB = uv(t(2),:);   vC = uv(t(3),:);
    eBC = vB-vC;  eBA = vB-vA;  eCA = vC-vA;  eCB = vC-vB;
    [cotA,cotB,cotC] = getAngles( mesh3D, t(1), t(2), t(3) );
    num = cotA*dot(eBC,eBC) + cotB*dot(eBA,eBA) + cotC*dot(eCA,eCA);
    denom = cross2(eBA,eCA);
    E = E + num/denom;
end
function [ E ] = MIPSEnergyV_vec( uv_i, vi, mesh3D, uv, faces )
tmpuv = uv;
tmpuv(vi,:) = uv_i;
E = MIPSEnergyV( mesh3D, tmpuv, faces );
   
% compute gradient of MIPS energy
% This function assumes that only the first vertex in each
% face (which must be the same) is free. It computes the
% gradient of that vertex wrt it's fixed boundary ring
function [ dUV, E ] = MIPSEnergyVGrad( mesh3D, uv, faces )
E = 0;
dUV = [0,0];
nTris = size(faces,1);
for ti = 1:nTris
    t = faces(ti,:);
    vA = uv(t(1),:);  vB = uv(t(2),:);   vC = uv(t(3),:);
    eBC = vB-vC;  eBA = vB-vA;  eCA = vC-vA;  eCB = vC-vB;
    [cotA,cotB,cotC] = getAngles( mesh3D, t(1), t(2), t(3) );

    % compute numerator & denominator of energy
    num = cotA*dot(eBC,eBC) + cotB*dot(eBA,eBA) + cotC*dot(eCA,eCA);
    denom = cross2(eBA,eCA);

    % compute derivative wrt uvA
    dA = -2*(cotB*eBA + cotC*eCA)/denom - perpdot(eCB)*(num/denom^2);
    
    dUV = dUV + dA;
    E = E + num/denom;
end
function [ E, dUV ] = MIPSEnergyVGrad_vec( uv_i, vi, mesh3D, uv, faces )
tmpuv = uv;
tmpuv(vi,:) = uv_i;
[dUV, E] = MIPSEnergyVGrad( mesh3D, tmpuv, faces );
dUV = dUV';



% compute face angles for face [t1,t2,t3]
function [ cotA, cotB, cotC ] = getAngles( mesh, t1, t2, t3 )
    vA = mesh.v(t1,:);  vB = mesh.v(t2,:);   vC = mesh.v(t3,:);
    cosAlpha = dot( vB-vA, vC-vA );
    cosBeta = dot( vA-vC, vB-vC );
    cosGamma = dot( vA-vB, vC-vB );
    cotA = cot( acos( clamp(cosAlpha,-1,1) ) );
    cotB = cot( acos( clamp(cosBeta,-1,1) ) );
    cotC = cot( acos( clamp(cosGamma,-1,1) ) );



% compute all face angles in mesh
function [ tAnglesV ] = getAnglesV( mesh )
nt = size(mesh.f,1);
tAnglesV = [];
for ti = 1:nt
    t = mesh.f(ti,:);
    [cA,cB,cC] = getAngles( mesh, t(1), t(2), t(3) );
    tAnglesV(ti,:) = [cA,cB,cC];
end    
    
% compute MIPS energy for entire mesh
function [ E ] = MIPSEnergy( uv, faces, tAngles )
E = 0;
nt = size(faces,1);
for ti = 1:nt
    t = faces(ti,:);
    vA = uv(t(1),:);  vB = uv(t(2),:);   vC = uv(t(3),:);
    eBC = vB-vC;  eBA = vB-vA;  eCA = vC-vA;  eCB = vC-vB;
    cotA = tAngles(ti,1);   cotB = tAngles(ti,2);  cotC = tAngles(ti,3);
    num = cotA*dot(eBC,eBC) + cotB*dot(eBA,eBA) + cotC*dot(eCA,eCA);
    denom = cross2(eBA,eCA);
    E = E + num/denom;
end



