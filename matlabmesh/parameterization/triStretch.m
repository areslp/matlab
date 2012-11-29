function [ L2, Linf, TareaV, TareaU ] = triStretch( v, u, t )
%TRISTRETCH compute per-triangle stretch metrics 
% Based on formuls in Sander et al 2002 "Texture Mapping Progressive Meshes"
%   
%   v - 3D vertex positions
%   u - 2D vertex positions
%   t - triangle indices

nv = size(v,1);
nt = size(t,1);

L2 = zeros(nt,1);
Linf = zeros(nt,1);
TareaV = zeros(nt,1);
TareaU = zeros(nt,1);

for i = 1:nt
    tu = u( t(i,:), : );
    p1 = tu(1,:);  p2 = tu(2,:);  p3 = tu(3,:);
    s1 = p1(1,1);   t1 = p1(1,2);
    s2 = p2(1,1);   t2 = p2(1,2);
    s3 = p3(1,1);   t3 = p3(1,2);
    
    tv = v( t(i,:), : );
    q1 = tv(1,:);  q2 = tv(2,:);  q3 = tv(3,:);

    A = ( (s2-s1)*(t3-t1) - (s3-s1)*(t2-t1) ) / 2;
    Ss = ( q1*(t2-t3) + q2*(t3-t1) + q3*(t1-t2) ) / (2*A);
    St = ( q1*(s3-s2) + q2*(s1-s3) + q3*(s2-s1) ) / (2*A);

    a = dot(Ss,Ss);
    b = dot(Ss,St);
    c = dot(St,St);
    
    det = sqrt( (a-c)^2 + 4*b^2 );
    sMax = sqrt( 0.5 * (a + c + det) );     % max singular value
    sMin = sqrt( 0.5 * (a + c - det) );     % min singular value
    
    TareaV(i) = triarea( q1, q2, q3 );
    TareaU(i) = A;
    L2(i) = sqrt(0.5 * (a+c));
    Linf(i) = sMax;
end


end
