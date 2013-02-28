function [W A] = cotlp_matrix(filename)
%
% Compute the Laplace-Beltrami matrix from mesh by cot scheme
%
% INPUTS
%  filename:  off file of triangle mesh.
%
% OUTPUTS
%  W: weight matrix 
%  A: area weight per vertex, the Laplace-Beltrami matrix = diag(1./ D) * W


if nargin < 1
    error('Too few input arguments');	 
end;

[II JJ SS AA] = cotlpmatrix(filename);
W=sparse(II, JJ, SS);
A=AA;

