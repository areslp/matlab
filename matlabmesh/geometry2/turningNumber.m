function [ num ] = turningNumber( v )
% [ num ] = turningNumber( v )
%   compute turning number of closed polygon
%   v is Nx2 list of points (first/last not duplicated)

N = size(v,1);
anglesum = 0;
for k1 = 1:N;
    k2 = mod(k1,N) + 1;
    k3 = mod(k1+1,N) + 1;
    e1 = normalize(vadd( v(k2,:), -v(k1,:) ));
    e2 = normalize(vadd( v(k3,:), -v(k2,:) ));
    asign = sign(cross2(e1,e2));
    anglesum = anglesum + asign*vangle(e1,e2);
end
num = round(anglesum / (2*pi));



end
