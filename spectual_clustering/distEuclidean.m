function [ d ] = distEuclidean( M, N )
%DISTEUCLIDEAN Calculates Euclidean distances
%   distEuclidean calculates the Euclidean distances between n
%   d-dimensional points, where M and N are d-by-n matrices, and
%   returns a 1-by-n vector d containing those distances.
%
%   Author: Ingo Buerk
%   Year  : 2011/2012
%   Bachelor Thesis

d = sqrt(sum((M - N) .^ 2, 1));

end

