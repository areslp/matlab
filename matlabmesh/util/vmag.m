function [ Mout ] = vmag( M )
%VMAG compute magnitude of rows in M

Mout = sqrt( sum( M .* M, 2 ) );

end
