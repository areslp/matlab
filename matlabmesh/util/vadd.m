function [ Mout ] = vadd( M, v )
%VADD adds v to each row of M

Mout = M + repmat(v, size(M,1), 1);

end
