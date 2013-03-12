%SPARSE_NN
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology

function [edgesrow, edgescol,edgesdist] = sparse_nn(snn)
    % turn into sparse nearest neighbor graph snn into edgesrow and edgescol index
    N = size(snn,1);
    edgescol = zeros(N+1,1);
    nnzer = nnz(snn);
    edgesrow = zeros(nnzer,1);
    edgesdist = zeros(nnzer,1);

    edgescol(1) = 0;
    for jdx=1:N
        lst = find(snn(:, jdx)>0);
        %lst = lst(find(lst>jdx));
        edgescol(jdx+1) = edgescol(jdx)+length(lst);
        edgesrow(edgescol(jdx)+1:edgescol(jdx+1)) = lst-1;
        edgesdist(edgescol(jdx)+1:edgescol(jdx+1))=snn(lst,jdx);
    end