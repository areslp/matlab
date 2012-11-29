function [ Mout, Mlen ] = normalize( M )
%NORMALIZE normalize rows of M. If magnitude is 0, skip vector
% should we set degenerate vectors to 0 ??

mag = sum(M .* M, 2);
z = find(mag < eps);
mag(z) = 1;
Mlen = sqrt(mag);
Mout = M ./ repmat(Mlen, 1, size(M,2));

