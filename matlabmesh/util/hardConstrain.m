function [ Mcons, RHScons ] = hardConstrain( M, RHS, consi, consv, rewrite )
%[ Mcons, RHScons ] = hardConstrain( M, RHS, consi, consv, rewrite )
%   add hard constraints by rewriting linear system 
%   consi = row indices of values to constrain
%   consv = constrained values
%   rewrite_system = 1 if variables in consi should be removed from system 
%      (this can improve the system condition number a large amount)

if ~ issparse(M)
    fprintf('hardConstrain(): matrix is non-sparse, should handle this case\n');
end
nrhs = size(RHS,2);
[nrows,ncols] = size(M);
if nrows ~= ncols
    error('hardConsrain(): can only add hard constraints for square systems');
end

if ~exist('rewrite','var')
    rewrite = 0;
end


Mcons = M;
RHScons = RHS;
ncons = numel(consi);
for ni = 1:ncons
    row = consi(ni);
    
    % set row to 1=value, leave other rows
    Mcons(row,:) = 0;
    Mcons(row,row) = 1;
    RHScons(row,:) = consv(ni,:);

    % move values over to RHS so we can make system smaller
    if rewrite
        rows = find(Mcons(:,row));
        rows = rows(rows~=row);
        Mcons(row,row) = 1;
        col = Mcons(rows,row);
        col = repmat(col,1,nrhs);
        val = repmat(consv(ni,:),numel(rows),1);
        RHScons(rows,:) = RHScons(rows,:) - (col.*val);
        Mcons(rows,row) = 0;
    end
end

% remove empy rows and columns from system
if rewrite
    keep = setdiff(1:nrows,consi);
    Mcons = Mcons(keep,keep);
    RHScons = RHScons(keep,:);
end


end
