function [ ] = writeMesh( mesh, filename )
% [ ] = writeMesh( mesh, filename )
%   write out mesh in .OBJ format


fid = fopen(filename, 'w');

has_tex = ( size(mesh.u,1) == size(mesh.v,1) );
has_normal = ( size(mesh.n,1) == size(mesh.v,1) );

for i = 1:numel(mesh.vidx)
    fprintf(fid, 'v %f %f %f\n', mesh.v(i,1), mesh.v(i,2), mesh.v(i,3));
end

if has_normal
    for i = 1:numel(mesh.vidx)
        fprintf(fid, 'vn %f %f %f\n', mesh.n(i,1), mesh.n(i,2), mesh.n(i,3));
    end
end

if has_tex
    for i = 1:numel(mesh.vidx)
        fprintf(fid, 'vt %f %f\n', mesh.u(i,1), mesh.u(i,2));
    end
end

for i = 1:numel(mesh.fidx)
    if has_tex & has_normal
        fprintf(fid, 'f %d/%d/%d %d/%d/%d %d/%d/%d\n', mesh.f(i,1),mesh.f(i,1),mesh.f(i,1),  mesh.f(i,2),mesh.f(i,2),mesh.f(i,2),  mesh.f(i,3),mesh.f(i,3),mesh.f(i,3)  );
    elseif has_normal
        fprintf(fid, 'f %d//%d %d//%d %d//%d\n', mesh.f(i,1),mesh.f(i,1),  mesh.f(i,2),mesh.f(i,2),  mesh.f(i,3),mesh.f(i,3)  );
    elseif has_tex
        fprintf(fid, 'f %d/%d %d/%d %d/%d\n', mesh.f(i,1),mesh.f(i,1),  mesh.f(i,2),mesh.f(i,2),  mesh.f(i,3),mesh.f(i,3)  );
    else
        fprintf(fid, 'f %d %d %d\n', mesh.f(i,1),mesh.f(i,2),mesh.f(i,3)  );
    end
end

fclose(fid);


end
