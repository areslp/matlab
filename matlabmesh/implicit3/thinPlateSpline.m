function [ xi, yi, F ] = thinPlateSpline( X, V, fSmooth )
%THINPLATESPLINE Summary of this function goes here
%   Detailed explanation goes here

N = numel(V);

if ~ exist('fSmooth')
    fSmooth = 0;
end

%minxy = min( min(X(:,1)), min(X(:,2)) );
%maxxy = max( min(X(:,1)), max(X(:,2)) );
minxy = -2;
maxxy = 2;

M = zeros(N,N);
B = V;
for i = 1:N
    pt = X(i,:);
    dists = vmag(vadd(X,-pt));
    vals = dists .* dists .* log(dists);
    vals(isnan(vals)) = 0;
    M(i,:) = vals;
end

M = M + eye(N,N)*fSmooth;

W = M \ B;

nsteps = 100;
k = (maxxy-minxy)/nsteps;

[xi,yi] = meshgrid(minxy:k:maxxy);
F = xi;
Ng = size(xi,1);
for i = 1:Ng
    for j = 1:Ng
        pt = [ xi(i,j), yi(i,j) ];
        dists = vmag(vadd(X,-pt));
        vals = dists .* dists .* log(dists);
        vals(isnan(vals)) = 0;
        f = sum( W .* vals );

        F(i,j) = f;
    end
end

        

