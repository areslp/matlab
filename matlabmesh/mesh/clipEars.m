function [ mesh_out ] = clipEars( mesh )
%CLIPEARS clip 'ear' triangles from mesh
%   (an ear triangle is a boundary triangle only connected on one edge)

mesh_out = mesh;

MAX_ITER = 10;
for iter = 1:MAX_ITER

    % 2 out of 3 edges of an ear triangle do not have a neighbour
    [vi,vj] = find(mesh_out.e == 1);
    
    isear = zeros(size(mesh_out.f,1), 1);
    for ii = 1:numel(vi)
        i = vi(ii);
        j = vj(ii);

        t1 = mesh_out.te(i,j);
        t2 = mesh_out.te(j,i);
        if t1 ~=0 && t2 ~= 0
            error('bug');
        end
        t = max(t1,t2);

        f = mesh_out.f(t,:);
        count = 0;
        if mesh_out.e(f(1),f(2)) == 1
            count = count+1;
        end
        if mesh_out.e(f(2),f(3)) == 1
            count = count+1;
        end
        if mesh_out.e(f(3),f(1)) == 1
            count = count+1;
        end
        
        if count == 2
            isear(t) = 1;
        end
    end

    if sum(isear) == 0
        break;
    end
    
    face_roi = mesh_out.fidx(isear==0);
    
    mesh_out = subMesh( mesh_out, face_roi );
end

if iter == MAX_ITER
    fprintf('[clipEars WARNING] Ear triangles still exist after %d iterations\n', MAX_ITER)
end


end
