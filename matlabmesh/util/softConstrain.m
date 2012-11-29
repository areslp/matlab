function [ Mcons, RHScons ] = softConstrain( M, RHS, consi, consv, consw )
%[ Mcons, RHScons ] = softConstrain( M, RHS, consi, consv )
%   add soft constraints to linear system
%   consi = row indices of values to constrain
%   consv = constrained values
%   consw = constrained row weights (default=1)


if ~ issparse(M)
    fprintf('softConstrain(): matrix is non-sparse, should handle this case\n');
end

[nrows,ncols] = size(M);
nrhs = size(RHS,2);
ncons = numel(consi);

if nrows ~= ncols
    error('softConstrain(): can only add soft constraints for square systems');
end
if ~exist('consw','var')
    consw = ones(ncons,1);
end


Mcons = sparse(consi,consi,consw.^2,nrows,ncols);
RHScons = zeros(size(RHS));
for ni = 1:ncons
    RHScons(consi(ni),:) = consw(ni).^2 * consv(ni,:);
end

Mcons = M + Mcons;
RHScons = RHS + RHScons;

end
