function [ Mv ] = mvmul( M, v )
%[ Mv ] = mvmul( M, v ) 
%   multiply row-vectors v by matrix M and return row-vector

Mv = (M * v')';

end
