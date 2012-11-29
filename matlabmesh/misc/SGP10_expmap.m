
%% compute expmap at north pole on set of hemisphere meshes

meshes = {'hemisphere_graphite_1000.obj', 'hemisphere_graphite_2000.obj', 'hemisphere_graphite_4000.obj', 'hemisphere_graphite_16000.obj'};
%meshes = {'hemisphere_graphite_1000.obj', 'hemisphere_graphite_2000.obj', 'hemisphere_graphite_4000.obj'};
northpole = [0,0,1];
radius = 1;
[tan1,tan2] = tangentFrame(normalize(northpole));

mesh_G = {};
expmap_G_err = {};
upwind_G_err = {};
expmap_T_err = {};
upwind_T_err = {};
expmap_L2_err = {};
upwind_L2_err = {};
    
for mi = 1:numel(meshes)
    meshname = meshes{mi};
    hemimesh = readMesh(meshname, 'n');
    hemimesh.v = radius * normalize(hemimesh.v);

    meshuv = zeros(numel(hemimesh.vidx),2);
    for vi = 1:numel(hemimesh.vidx)
        meshuv(vi,:) = sphereNormCoords(radius, northpole, hemimesh.v(vi,:), tan1,tan2);
    end
    hemimesh.u = meshuv;
    meshG = vmag(meshuv);   meshT = atan2(meshuv(:,2),meshuv(:,1));
    [tmp,idx] = sort(meshG);


    hemips = [];
    hemips.v = [ northpole; hemimesh.v ];
    hemips.n = normalize(hemips.v);

    options = [];
    options.graphnbrs = 8;
    options.stopcrit = 'none';
    options.stopparam = inf;
    options.upwindavg = 0;
    expmapuv = expmap(hemips.v, hemips.n, 1, options);
    expmapuv = expmapuv(2:end,:);
    rotangle = meshT(1)-atan2(expmapuv(1,2),expmapuv(1,1));
    rotmat = [cos(rotangle),sin(rotangle); -sin(rotangle),cos(rotangle)];
    expmapuv = mvmul(expmapuv,rotmat)';
    expmapG = vmag(expmapuv);  expmapT = atan2(expmapuv(:,2),expmapuv(:,1));

    options.upwindavg = 1;
    upwinduv = expmap(hemips.v, hemips.n, 1, options);
    upwinduv = upwinduv(2:end,:);
    rotangle = meshT(1)-atan2(upwinduv(1,2),upwinduv(1,1));
    rotmat = [cos(rotangle),sin(rotangle); -sin(rotangle),cos(rotangle)];
    upwinduv = mvmul(upwinduv,rotmat)';
    upwindG = vmag(upwinduv);  upwindT = atan2(upwinduv(:,2),upwinduv(:,1));
    
    mesh_G{mi} = meshG;
    
    % absolute errors |u-u'|
    expmap_G_err{mi} = abs(meshG-expmapG);
    upwind_G_err{mi} = abs(meshG-upwindG);
    expmap_T_err{mi} = abs(angledist(meshT,expmapT));
    upwind_T_err{mi} = abs(angledist(meshT,upwindT));
    expmap_L2_err{mi} = vmag(meshuv-expmapuv);
    upwind_L2_err{mi} = vmag(meshuv-upwinduv));
    
    % scaled errors |u-u'| / |u|
%     expmap_G_err{mi} = abs(meshG-expmapG) ./ meshG;
%     upwind_G_err{mi} = abs(meshG-upwindG) ./ meshG;
%     expmap_T_err{mi} = abs(angledist(meshT,expmapT));   % // how to scale
%     upwind_T_err{mi} = abs(angledist(meshT,upwindT));
%     expmap_L2_err{mi} = vmag(meshuv-expmapuv) ./ vmag(meshuv);
%     upwind_L2_err{mi} = vmag(meshuv-upwinduv) ./ vmag(meshuv);    
end

M = numel(meshes);

use_idx = {};
for mi = 1:M
    [tmp,idx] = sort(mesh_G{mi});
    idx=idx(find(mesh_G{mi}(idx)<pi/4));    % cull to points within pi/4 cone
    use_idx{mi} = idx;
end

use_upwind = 1;



%% scatterplots  x=true geodesic dist, y=error

figure;hold all;
for mi = 1:M
    idx=use_idx{mi};
    if use_upwind
        scatter(mesh_G{mi}(idx),upwind_G_err{mi}(idx),3,'filled');
    else
        scatter(mesh_G{mi}(idx),expmap_G_err{mi}(idx),3,'filled');
    end
end
title('Geodist Error');
hold off; drawnow;


figure;hold all;
for mi = 1:M
    idx=use_idx{mi};
    if use_upwind
        scatter(mesh_G{mi}(idx),upwind_T_err{mi}(idx),3,'filled');
    else
        scatter(mesh_G{mi}(idx),expmap_T_err{mi}(idx),3,'filled');
    end
end
title('Angle Error');
hold off; drawnow;


figure;hold all;
for mi = 1:M
    idx=use_idx{mi};
    if use_upwind
        scatter(mesh_G{mi}(idx),upwind_L2_err{mi}(idx),3,'filled');
    else
        scatter(mesh_G{mi}(idx),expmap_L2_err{mi}(idx),3,'filled');
    end
end
title('L2 Error');
hold off; drawnow;



%% plots sorted by magnitude of error values  (gives smooth curves)

figure;hold all;
for mi = 1:M
    idx=use_idx{mi};
    if use_upwind
        y = sort(upwind_G_err{mi}(idx));
    else
        y = sort(expmap_G_err{mi}(idx));
    end
    xi = (1:numel(idx)) / numel(idx);
    plot(xi,y);
end
title('Sorted Geodist Error');
hold off; drawnow;


figure;hold all;
for mi = 1:M
    idx=use_idx{mi};
    if use_upwind
        y = sort(upwind_T_err{mi}(idx));
    else
        y = sort(expmap_T_err{mi}(idx));
    end
    xi = (1:numel(idx)) / numel(idx);
    plot(xi,y);
end
title('Sorted Angle Error');
hold off; drawnow;


figure;hold all;
for mi = 1:M
    idx=use_idx{mi};
    if use_upwind
        y = sort(upwind_L2_err{mi}(idx));
    else
        y = sort(expmap_L2_err{mi}(idx));
    end
    xi = (1:numel(idx)) / numel(idx);
    plot(xi,y);
end
title('Sorted L2 Error');
hold off; drawnow;



%% print mean/max errors


fprintf('[expmap avg]\n');
for mi = 1:M
    idx=use_idx{mi};
    fprintf('[%6d verts] G: %8.6f    T: %8.6f     L2: %8.6f\n', numel(idx), mean(expmap_G_err{mi}(idx)), mean(expmap_T_err{mi}(idx)), mean(expmap_L2_err{mi}(idx)) );
end
fprintf('[upwind avg]\n');
for mi = 1:M
    idx=use_idx{mi};
    fprintf('[%6d verts] G: %8.6f    T: %8.6f     L2: %8.6f\n', numel(idx), mean(upwind_G_err{mi}(idx)), mean(upwind_T_err{mi}(idx)), mean(upwind_L2_err{mi}(idx)) );
end




fprintf('[expmap max]\n');
for mi = 1:M
    idx=use_idx{mi};
    fprintf('[%6d verts] G: %8.6f    T: %8.6f     L2: %8.6f\n', numel(idx), max(expmap_G_err{mi}(idx)), max(expmap_T_err{mi}(idx)), max(expmap_L2_err{mi}(idx)) );
end
fprintf('[upwind max]\n');
for mi = 1:M
    idx=use_idx{mi};
    fprintf('[%6d verts] G: %8.6f    T: %8.6f     L2: %8.6f\n', numel(idx), max(upwind_G_err{mi}(idx)), max(upwind_T_err{mi}(idx)), max(upwind_L2_err{mi}(idx)) );
end



