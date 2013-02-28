function [ Y, svp, sv ] = singular_value_shrinkage_acc( X, tau, sv )
n=length(X);
if choosvd(n, sv) == 1
    % fprintf(1,'partial svd use lansvd\n');
    [U S V] = lansvd(X, sv, 'L');
else
    % fprintf(1,'matlab svd\n');
    [U S V] = svd(X, 'econ');
end
diagS = diag(S);
svp = length(find(diagS > tau));
if svp < sv
    sv = min(svp + 1, n);
else
    sv = min(svp + round(0.05*n), n);
end
Y = U(:, 1:svp) * diag(diagS(1:svp) - tau) * V(:, 1:svp)';
