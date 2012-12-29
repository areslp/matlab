function [ M ] = simGaussian( M, sigma )
%SIMGAUSSIAN Calculates Gaussian similarity on matrix
%   simGaussian(M, sigma) returns a matrix of the same size as
%   the distance matrix M, which contains similarity values
%   that are computed by using a Gaussian similarity function
%   with parameter sigma.
%
%   Author: Ingo Buerk
%   Year  : 2011/2012
%   Bachelor Thesis

M = exp(-M.^2 ./ (2*sigma^2));

end
