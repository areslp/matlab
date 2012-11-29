function [ area ] = meshArea( verts, tris )
%MESHAREA compute area of triangles defined by vertices
%  This function works for 2D and 3D embeddings.
%
%  verts is NxD list of vertices
%  tris is Nx3 list of vertex indices
%

[nv,nd] = size(verts);
nt = size(tris,1);

area = 0;
for i = 1:nt
    tv = verts( tris(i,:), : );
    
    if nd == 2
        vA = tv(1,:);   vB = tv(2,:);   vC = tv(3,:);
        eBA = vB-vA;  eCA = vC-vA;
        area = area + 0.5 * cross2(eBA,eCA);        
    elseif nd == 3
        vA = tv(1,:);   vB = tv(2,:);   vC = tv(3,:);
        area = area + triarea( vA, vB, vC );
    end
    
end

end
