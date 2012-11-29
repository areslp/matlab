function [ G ] = findKNN( points, k )
% [ G ] = findKNN( points, k )
%   G is sparse matrix, row G(i,:) is k-nearest-nbr dists for pts(i)

[dim,N] = size(points');
Kgraph = k;
G = spalloc(N, N, Kgraph * N);


if N < 2000
   
    X = points';
    X2 = sum(X.^2,1);
    distance = sqrt( repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X );
    [sorted,index] = sort(distance);
    for jj = 1:N
        nbrs = index(2:(1+Kgraph),jj);
        nbrdists = sorted(2:(1+Kgraph),jj);
        G(jj,nbrs) = nbrdists;
    end
    
else
    
    for jj = 1:N
        dists = vmag( vadd(points,-points(jj,:)) );
        [sorted,index] = sort(dists);
        nbrs = index(2:(1+Kgraph));
        nbrdists = sorted(2:(1+Kgraph));
        G(jj,nbrs) = nbrdists;
    end
    
end

end
