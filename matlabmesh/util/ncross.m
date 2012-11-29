function [ nc ] = ncross( v1, v2 )
%NCROSS compute normalized cross-product
    nc = normalize( cross( v1, v2 ) );
end
