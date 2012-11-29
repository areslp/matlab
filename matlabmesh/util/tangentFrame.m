function [ tan1, tan2 ] = tangentFrame( normal )
%TANGENTFRAME create orthogonal tangent vectors to normal

n = normalize(normal);

if ( abs(n(1)) >= abs(n(2)) && abs(n(1)) >= abs(n(3)) ) 
    tan1 = [ -n(2), n(1), 0 ];
else
    tan1 = [ 0, n(3), -n(2) ];
end
tan1 = normalize(tan1);
tan2 = cross( tan1, n );
    
end
