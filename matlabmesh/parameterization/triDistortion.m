function [ Eangle, Earea, Tarea ] = triDistortion( v, u, t )
%MEASURETRIDISTORTION measure angle and area distortion of a triangle

nt = size(t,1);

Eangle = zeros(nt,1);
Earea = zeros(nt,1);
Tarea = zeros(nt,1);

for i = 1:nt
    tv = v( t(i,:), : );
    tu = u( t(i,:), : );

    vA = tv(1,:);   vB = tv(2,:);   vC = tv(3,:);
    cosAlpha = dot( vB-vA, vC-vA );
    cosBeta = dot( vA-vC, vB-vC );
    cosGamma = dot( vA-vB, vC-vB );
    cotA = cot( acos( clamp(cosAlpha,-1,1) ) );
    cotB = cot( acos( clamp(cosBeta,-1,1) ) );
    cotC = cot( acos( clamp(cosGamma,-1,1) ) );    
    areav = triarea( vA, vB, vC );
    
    vA = tu(1,:);   vB = tu(2,:);   vC = tu(3,:);
    eBC = vB-vC;  eBA = vB-vA;  eCA = vC-vA;
    num = cotA*dot(eBC,eBC) + cotB*dot(eBA,eBA) + cotC*dot(eCA,eCA);
    areau = 0.5 * cross2(eBA,eCA);

    Tarea(i) = areav;
    Eangle(i) = num / (2 * areau);
    Earea(i) = areav/areau + areau/areav;
end


end
