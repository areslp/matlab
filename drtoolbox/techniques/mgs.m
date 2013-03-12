function [Q, R] = mgs(A)
%MGS Performs a modified Gram-Schmidt orthogonalization
%
%   [Q, R] = mgs(A)
%
% Performs a modified Gram-Schmidt orthogonalization. This is a more stable
% way to compute the Gram-Schmidt orthogonalization than Matlab's one.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology

   
    % Perform Gram-Schmidt orthogonalization
    [m, n] = size(A);
    V = A;
    R = zeros(n, n);
    for i=1:n
        R(i,i) = norm(V(:,i));
        V(:,i) = V(:,i) / R(i, i);
        if (i < n)
            for j=i+1:n
                R(i,j) = V(:,i)' * V(:,j);
                V(:,j) = V(:,j) - R(i, j) * V(:,i);
            end
        end
    end
    Q = V;
