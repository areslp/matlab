function [ A ] = triarea( p1, p2, p3 )
% [ A ] = triarea( p1, p2, p3 )
%   area of triangle defined by 3 points
    u = p2 - p1;
    v = p3 - p1;
    
    % 3D formula
    %A = 0.5 * vmag( cross(u,v) );

    %   n-D formula from [Meyer02] eq 13
    A = 0.5 * sqrt( dot(u,u)*dot(v,v) - dot(u,v)^2 );
end
