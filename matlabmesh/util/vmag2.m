function [ Mout ] = vmag2( M )
%VMAG compute sqiared magnitude of rows in M

Mout = sum( M .* M, 2 );

end