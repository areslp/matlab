function [ faceValues ] = faceDistortion( mesh, type )
%FACEDISTORTION Summary of this function goes here
%   Detailed explanation goes here

N = numel(mesh.fidx);
faceValues = zeros(N,1);

if strcmp(type, 'dirichlet')   % dirichlet energy 
    for fi = 1:N
        f = mesh.f(fi,:);
        fsum2 = 0;
        A2 = triarea( mesh.u(f(1),:), mesh.u(f(2),:), mesh.u(f(3),:) );
        for vi = 1:3
            i = f(vi);
            j = f( mod(vi,3) + 1 );
            k = f( f~=i & f~=j );
            e1 = normalize(mesh.v(i,:) - mesh.v(k,:));
            e2 = normalize(mesh.v(j,:) - mesh.v(k,:));
            dist2 = vmag2(mesh.u(i,:) - mesh.u(j,:));
            fsum2 = fsum2 + (1/4)*vcot(e1,e2)*dist2;
        end
        faceValues(fi) = fsum2/A2;
    end
    
elseif strcmp(type,'qc')   % quasi-conformal
    t = mesh.f;
    u = mesh.u;
    v = mesh.v;
    for i = 1:N
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
        faceValues(i) = sMax / sMin;
    end    
    
end

end
