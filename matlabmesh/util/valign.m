function [ rotmat, axis, theta ] = valign( vfrom, vto )
%[ rotmat, axis, theta ] = valign( vfrom, vto )
%   generate matrix that rotates vector vfrom into vector vto

n1 = vfrom;
n2 = vto;
costheta = dot(n1,n2) / (norm(n1)*norm(n2));
axis = cross(n1,n2);
if ( norm(axis) > eps )
    axis = axis / norm(axis);
    theta = acos( clamp(costheta,-1,1) );
    rotmat = axisrot( axis, theta );
elseif costheta < 0
    [p1,p2] = tangentFrame(n1);
    axis = p1;
    theta = pi;
    rotmat = axisrot(axis, theta);
else
    axis = [0,0,1];
    theta = 0;
    rotmat = eye(3);
end


end
