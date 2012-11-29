function [ N ] = trinormal( p1, p2, p3 )
% [ N ] = trinormal( p1, p2, p3 )
    N = ncross( p2-p1, p3-p1 );
end
