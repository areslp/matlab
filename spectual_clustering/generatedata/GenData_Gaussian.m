function [ T ] = GenData_Gaussian( num, mu, sigma )
%GENDATA_GAUSSIAN Calculates a Gaussian sample data set
%   GenData_Gaussian returns a matrix containing data points
%   as columns, which are Gaussian distributed around mu with
%   variance matrix sigma.
%
%   For example, GenData_Gaussian(500, [0 0], [1 0; 0 1]) will
%   return a set of 500 data points centered around the origin.
%
%   'num' - Number of data points to be created
%   'mu' - vector containing the center of the data set
%   'sigma' - variance matrix
%
%   Author: Ingo Buerk
%   Year  : 2011/2012
%   Bachelor Thesis

% make mu row vector
mu = mu(:);

% compute covariance matrix
Cov = chol(sigma);

% create Gaussian data set
T = repmat(mu, 1, num) + Cov * randn(size(sigma, 1), num);

end

