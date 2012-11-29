function [ H, data ] = meanCurv( mesh, vertices, mode )
%MEANCURV estimate mean curvature at vertices
%   [Ryan Schmidt  rms@dgp.toronto.edu  09/2008]
%   - if vertices is undefined/empty, compute for entire mesh
%   - returned data is mode specific
%   - valid modes are:
%     'normal' - Schneider & Kobbelt 01 version of Moreton & Sequin 92
%        ( data is list of [ qj_x, qj_y, qj_z, tj_x, tj_y ] for last vertex )
%     'dmc' - discrete mean-curvature normal operator  [Meyer02]
%        ( data is list of mean-curvature normals )

if ~ exist('mode', 'var')
    mode = 'normal';
end

if ~ exist('vertices', 'var') || numel(vertices) == 0
    vertices = 1:size(mesh.v,1);
end
n = numel(vertices);
   
H = zeros(n,1);
if strcmp(mode, 'normal')
    for i = 1:n
        [H(i), data] = meanCurv_SK01( mesh, vertices(i) );
    end
elseif strcmp(mode, 'dmc')
    data = zeros(n,3);
    for i = 1:n
        [H(i), tmpdata] = meanCurv_Meyer02( mesh, vertices(i) );
        data(i,:) = tmpdata;
    end
end

end


function [ H, data ] = meanCurv_Meyer02(mesh, v)
    W = cotanWeights(mesh, v);
    W = W(v,:);
    Amixed = vertexArea(mesh, v, 'mixed');
    
    Hn = zeros(1,3);
    nbrs = find(W ~= 0);
    for ni = 1:numel(nbrs)
        j = nbrs(ni);
        Hn = Hn + W(j)*( mesh.v(v,:) - mesh.v(j,:) );
    end
    Hn = Hn / (2*Amixed);
    
    H = vmag(Hn);
    data = Hn;
end


function [ H, data ] = meanCurv_SK01( mesh, v )
    [i,j] = find( mesh.e(v,:) );
    
    qi = mesh.v(v,:);
    ni = mesh.n(v,:);
    Qj = vadd( mesh.v(j,:), -qi );
    
    if mesh.valence(v) < 5
        newQ = [];
        [ot,ov] = onering(mesh,v);
        for k = 1:numel(ot)
            % figure out which two nbr verts to use
            f = mesh.f(ot(k),:);
            nbrs = f(f~=v);
            q1 = nbrs(1);  q2 = nbrs(2);
            f1 = mesh.te(q1,q2);   f2 = mesh.te(q2,q1);
            if f1 == 0 || f2 == 0
                continue;    % edge triangle
            end
            
            % construct E plane
            v1 = mesh.v(q1,:);  v2 = mesh.v(q2,:);
            Et1 = normalize(v2 - v1);
            Et2 = normalize( mesh.fn(f1,:) + mesh.fn(f2,:) );
            En = ncross( Et1, Et2 );
            
            % construct tangent vectors
            t1 = ncross( mesh.n(q1,:), En );
            if acos( vdot(-t1, Et1) ) < acos( vdot( t1, Et1) )
                t1 = -t1;
            end
            t2 = ncross( mesh.n(q2,:), En );
            if acos( vdot(-t2, Et1) ) < acos( vdot( t2, Et1) )
                t2 = -t2;
            end
            
            % construct basis for E plane and project points into it
            [bx,by] = tangentFrame( En );
            Eo = v1;
            v1 = [vdot((v1-Eo),bx), vdot((v1-Eo),by)];
            v2 = [vdot((v2-Eo),bx), vdot((v2-Eo),by)];
            t1 = normalize( [t1*bx', t1*by'] );
            t2 = normalize( [t2*bx', t2*by'] );
            
            % compute circles
            r1 = 2*vmag2(v2-v1) / ((v2-v1)*perpdot(t1)');
            c1 = circle2_2pr(v1, v2, r1);
            if normalize(v1-c1(1,:))*t1' > 0.1
                c1 = c1(2,:);
            else
                c1 = c1(1,:);
            end
            r2 = 2*vmag2(v1-v2) / ((v1-v2)*perpdot(t2)');
            c2 = circle2_2pr(v1, v2, r2);
            if normalize(v2-c2(1,:))*t2' > 0.1
                c2 = c2(2,:);
            else
                c2 = c2(1,:);
            end
            
            % [RMS] not sure why this happens...maybe triangle collapsed?
            if ~isfinite(r1) || ~isfinite(r2) 
                continue;
            end
            
            % intersect perp bisector with circles
            % [RMS] note that selection between two intersection
            %   points here is slightly different than specified
            %   in the paper - they use the nearest point to the
            %   line (v1->v2), while I'm using the nearest point
            %   to the midpoint 0.5(v1+v2). Seems like this should be
            %   equivalent where this hack is going to work anyway...
            pO = 0.5 * (v1 + v2);
            pD = perpdot(v2-v1);
            hit1 = isect_line2_circle2( pO, pO+pD, c1, r1 );
            if vmag2(hit1(1,:)-pO) < vmag2(hit1(2,:)-pO)
                hit1 = hit1(1,:);
            else
                hit1 = hit1(2,:);
            end
            hit2 = isect_line2_circle2( pO, pO+pD, c2, r2 );
            if vmag2(hit2(1,:)-pO) < vmag2(hit2(2,:)-pO)
                hit2 = hit2(1,:);
            else
                hit2 = hit2(2,:);
            end
            
            % project new point back to 3D
            p = 0.5 * (hit1 + hit2);
            P = Eo + p(1)*bx + p(2)*by;
            newQ = [newQ ; P - qi];
        end
        
%         newplot;    % draw points to see if they are reasonable...
%         hold all;
%         scatter3(Qj(:,1), Qj(:,2), Qj(:,3));
%         scatter3(newQ(:,1), newQ(:,2), newQ(:,3));
%         hold off;

        Qj = [Qj;newQ];
    end
    
    [bx, by] = tangentFrame( ni );
    Tj = [ vdot(Qj,bx), vdot(Qj,by) ];
    Tj = normalize(Tj);
    Tx = Tj(:,1);
    Ty = Tj(:,2);

    A = [Tx.^2, Tx.*Ty, Ty.^2];

    % [RMS] reverse normal here because our normals point 
    % outwards, while formulas assume pointing inwards (right?)
    K = 2 * vdot(Qj, -ni) ./ vdot(Qj,Qj);

    x = A \ K;

    H = (x(1) + x(3))/2;
    data = [Qj, Tj];
end
