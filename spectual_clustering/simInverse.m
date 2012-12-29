function [ M ] = simInverse( M )
%SIMINVERSE Calculates similarity as inverse distances
%   simInverse(M) returns a matrix of the same size as the
%   distance matrix M, which contains the inverse values, which
%   can be used as a measure for similarity.
%
%   Author: Ingo Buerk
%   Year  : 2011/2012
%   Bachelor Thesis

M = 1 ./ M;

end
