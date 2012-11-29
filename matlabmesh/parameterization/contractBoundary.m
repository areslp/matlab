function [ updateu ] = contractBoundary( mesh, boundaryUV, weights )
%CONTRACTBOUNDARY Summary of this function goes here
%   Detailed explanation goes here

do_reconstruct = 1;
smooth_iterations = 0;
smooth_alpha = 0.5;
do_offset = 0;
offset_size = 0.1;

% reconstruct boundary from weights
%  (current boundary is used...)
updateu = mesh.u;

if do_reconstruct
    for i = 1:size(mesh.v,1);
       if mesh.isboundaryv(i)
           nbrs = find(weights(i,:) ~= 0);
           newu = vdot( mesh.u(nbrs,:)', weights(i,nbrs) )';
           updateu(i,:) = newu;
%           updateu(i,:) = 0.5 * (newu + mesh.u(i,:));
       end
    end
end


% smooth new boundary
for SMOOTHITER = 1:smooth_iterations
    last_u = updateu;
    for ii = 1:size(mesh.loops)
       i = mesh.loops(ii);

       j = find(boundaryUV(:,1) == i);
       if j == size(boundaryUV,1)
           jn = 1; jp = j-1;
       elseif j == 1
           jn = j+1; jp = size(boundaryUV,1);
       else
           jn = j+1; jp = j-1;
       end
       j = boundaryUV(j,1);
       jn = boundaryUV(jn,1);
       jp = boundaryUV(jp,1);
       newu = 0.5 * ( last_u(jn,:) + last_u(jp,:) );
       newu = (1-smooth_alpha) * last_u(j,:) + smooth_alpha * newu;
       updateu(i,:) = newu;
    end
end


if do_offset
    n = size(mesh.loops,1);
    k1 = [1:n];
    k2 = [2:n,1];
    vdiff = mesh.v(mesh.loops(k1),:) - mesh.v(mesh.loops(k2),:);
    udiff = updateu(mesh.loops(k1),:) - updateu(mesh.loops(k2),:);
    perimv = sum( sqrt( sum( vdiff .* vdiff, 2) ) );
    perimu = sum( sqrt( sum( udiff .* udiff, 2) ) );
    bscale = perimv / perimu;

    % apply normal offset
    prev_u = updateu;
    for ii = 1:size(mesh.loops)
        i = mesh.loops(ii);

        j = find(boundaryUV(:,1) == i);
        if j == size(boundaryUV,1)
           jn = 1; jp = j-1;
        elseif j == 1
           jn = j+1; jp = size(boundaryUV,1);
        else
           jn = j+1; jp = j-1;
        end
        j = boundaryUV(j,1);
        jn = boundaryUV(jn,1);
        jp = boundaryUV(jp,1);
        un = prev_u( jn,: ) - prev_u( jp,: );
        udist = vmag(prev_u(jn,:)-prev_u(j,:)) + vmag(prev_u(jp,:)-prev_u(j,:));
        vdist = vmag(mesh.v(jn,:)-mesh.v(j,:)) + vmag(mesh.v(jp,:)-mesh.v(j,:));
        ratio = vdist / (udist * bscale);
        ratio = max( ratio-1, 0 );
        un = -normalize( perpdot(un) );

        updateu(i,:) = updateu(i,:) + offset_size * ratio * un;
    end
end