function [ output ] = errnorm( values, truth, p )
%ERRNORM Summary of this function goes here
%   Detailed explanation goes here

if ~ exist('p', 'var')
	p = 2;
end

output = norm(values-truth,p) / norm(truth,p);

end
