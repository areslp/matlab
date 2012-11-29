function [g1, g2, g3] = trigrad( p1, p2, p3 )
% [g1, g2, g3] = trigrad( p1, p2, p3 )
%   g1/g2/g3 are barycentric basis function gradients

v = [p1;p2;p3];
face = [1,2,3];
grad = zeros(3,3);

% construct gradient using geometry
% grad for i is perp to edge jk and tri normal, length is dist from i to jk 
for t = 1:3
    i = face(t);
    [j,k] = tripick(face, i);
    eji = v(i,:) - v(j,:);
    ejk = normalize( v(k,:) - v(j,:) );
    pproj = v(j,:) + vdot(eji,ejk)*ejk;   % projection of eji onto ejk
    vperp = v(i,:) - pproj;
    grad(t,:) = vperp / vmag2(vperp);       % direction * 1/length
end

% [RMS] this is the matrix form from the survey (after eq14)
%     matp = [ v(face(1),:) - v(face(3),:);  ...  
%              v(face(2),:) - v(face(3),:);   ...
%              faceNormal(mesh,ti) ];
%     matrhs = [ [1,0,-1];[0,1,-1];[0,0,0] ];
%     grad = (matp \ matrhs)';

g1 = grad(1,:);
g2 = grad(2,:);
g3 = grad(3,:);

end
