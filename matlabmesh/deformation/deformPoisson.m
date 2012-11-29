function [ deformed_mesh ] = deformPoisson( mesh, triTrans, consP )
%[ deformed_mesh ] = deformPoisson( mesh, triTrans, consP )
%
%   triTrans:     per-triangle transformation matrices
%                   - Tx9 matrix with rows  [x11,x12,x12,x21,...,x33]
%   consP:        rows of [vtx_i, x, y, z, weight_i]
%     (NOTE: hard constraints currently used, so weight is ignored...)

N = numel(mesh.vidx);
T = numel(mesh.fidx);

% construct 3TxN linear differential operator G, such that G*f = grad(f)
G = sparse([],[],[],3*T,N,3*T*3);
for ti = 1:T
    face = mesh.f(ti,:);
    row = 3*(ti-1)+1;
    
    % gradients of barycentric functions
    [g1,g2,g3] = trigrad( mesh.v(face(1),:), mesh.v(face(2),:), mesh.v(face(3),:) );    
    grad = [g1;g2;g3];
    for k = 0:2
        G(row+k,face) = grad(:,k+1);
    end
end

% compute gradients for each vertex coordinate function   (sanity check)
% gx/gy/gz are stacked vectors [f1_x,f1_y,f1_z,f2_x,...]
% gx = G*mesh.v(:,1);
% gy = G*mesh.v(:,2);
% gz = G*mesh.v(:,3);
% gxyz = [gx,gy,gz]

% compute gradients of transformed triangles
gxyz = zeros(3*T,3);
for i = 1:T
    trans = reshape( triTrans(i,:), 3,3 );
    triv = mesh.v(mesh.f(i,:),:);
    triv = mvmul(trans,triv);
    [g1,g2,g3] = trigrad(triv(1,:),triv(2,:),triv(3,:));
    
    r = 3*(i-1)+1;
    for k = 1:3
        gxyz(r:r+2,k) = triv(1,k)*g1 + triv(2,k)*g2 + triv(3,k)*g3;
    end
end

% construct and solve system (G'*M*G)x = (G'*M)grad(x)

%  (M is weight/"mass" matrix)
%M = speye(3*T);
fA = faceArea(mesh);
M = sparse(1:3*T,1:3*T, reshape([fA,fA,fA]', 3*T, 1), 3*T,3*T);

A = G' * M * G;

% matrix above is same as cotan Laplacian
%A = -makeOneRingWeights(mesh,'cotan',0) / 2;     % divide by 2 is important here...
%A(1:N+1:N*N) = -sum(A,2); 

RHS = G' * M * gxyz;

% solve for new positions w/ hard constraints
consi = consP(:,1);
consv = consP(:,2:4);
X = hardConstrainSolve(A,RHS,consi,consv);


deformed_mesh = mesh;
deformed_mesh.v = [X(:,1), X(:,2), X(:,3)];
deformed_mesh.n = estimateNormal(deformed_mesh);

end
