% Matlab Toolbox for Dimensionality Reduction
% Version 0.8 18-APR-2012
%
%
%Main functions
%----------------------
% compute_mapping    Performs dimension reduction using the specified technique
% out_of_sample      Performs out-of-sample extension using trained DR technique
% out_of_sample_est  Performs out-of-sample extension by approximation
% reconstruct_data   Reconstruct data from low-dimensional representation
% intrinsic_dim      Estimates intrinsic dimensionality of the data
%
%Helper functions
%----------------------
% drgui              Shows a GUI that facilitates access to all functions
% generate_data      Generates artificial data manifolds
% prewhiten          Performs whitening of a data set
%
%Installation functions
%----------------------
% mexall             Compiles all C++ code in the toolbox
% test_toolbox       Tests all functionalities of the toolbox
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology