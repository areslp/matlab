function [ Mout ] = vmul( M, v )
%[ Mout ] = vmul( M, v )
%   Detailed explanation goes here

if size(v,2) == size(M,2)
    Mout = M .* repmat(v, size(M,1), 1);
else
    Mout = M .* repmat(v, 1, size(M,2));
end


end
