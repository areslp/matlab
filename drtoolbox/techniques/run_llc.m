function mappedX = run_llc(X, no_dims, k, no_analyzers, max_iterations, eig_impl)
%RUN_LLC Performs the LLC algorithm for dimensionality reduction
%
%   mappedX = run_llc(X, no_dims, k, no_analyzers, max_iterations)
%
% Performs the Locally Linear Coordination (LLC) algorithm to reduce the
% dimensionality of dataset X to no_dims dimensions. The variable k
% indicates the number of neighbors that is used in the nieghborhood graph.
% The variable no_analyzers indicates the number of factor analyzers that
% is used, and max_iterations indicates the maximum number of iterations of
% the EM-algorithm that estimates the factor analyzers.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    if ~exist('no_dims', 'var')
        no_dims = 2;
    end
    if ~exist('k', 'var')
        k = 12;
    end
    if ~exist('no_analyzers', 'var')
        no_analyzers = 40;
    end
    if ~exist('max_iterations', 'var')
        max_iterations = 200;
    end

    % Initialize some parameters
    tol = 1e-10;                        % regularization parameter
    min_std = 1e-3;                     % minimum STD of Gaussians

    % Computes mixture of factor analyzers (using EM and PCA)
    disp('Running EM algorithm and compute local factor analyzers...');
    [LX, MX, PX] = mppca(X, no_dims, no_analyzers, tol, max_iterations, min_std);
    % Variables contain respectively:
    %  - LX     Lowdimensional representations of X of all factor analyzers
    %  - MX     Means of factor analyzers of X
    %  - PX     Noise covariance of X
    
    % Construct Mixture of Factor Analyzers based on results
    disp('Constructing mixture of factor analyzers...');
    [R, Z] = infermfa(X, LX, MX, PX);
    % Variables contain respectively:
    %   - R    Responsibilities of components of MFA
    %   - Z    Mean of posteriors over latent variables of MFA
 
    % Run the Local Linear Coordination algorithm on the MoFA
    mappedX = llc(X, k, no_dims, R, Z, eig_impl);
    