function [mappedX, mapping] = charting(X, no_dims, no_analyzers, max_iterations, eig_impl)
%CHARTING Performs manifold charting on dataset X 
%
%   [mappedX, mapping] = charting(X, no_dims, no_analyzers, max_iterations, eig_impl)
%
% Performs manifold charting on dataset X to reduce its dimensionality to
% no_dims dimensions. The variable no_analyzers determines the number of
% local factor analyzers that is used in the mixture of factor analyzers
% (default = 40). The variable max_iterations sets the maximum number of
% iterations that is employed in the EM algorithm (default = 200).
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
    if ~exist('no_analyzers', 'var')
        no_analyzers = 40;
    end
    if ~exist('max_iterations', 'var')
        max_iterations = 200;
    end
    
    % Initialize some parameters
    tol = 1e-10;                        % regularization parameter
    min_std = 1e-3;                     % minimum STD of Gaussians
    kf = no_analyzers * (no_dims + 1);

    % Construct MFA model on the data
    disp('Running EM algorithm and compute local factor analyzers...');
    [LX, MX, PX] = mppca(X', no_dims, no_analyzers, tol, max_iterations, min_std);
    [R, Z] = infermfa(X', LX, MX, PX);
    
    % Adds last entry = 1 in posterior mean to handle means of factor analyzers
    Z(no_dims + 1,:,:) = 1; 
    Z = permute(Z, [1 3 2]);
    
    % Construct blockdiagonal matrix D
    disp('Performing manifold charting...');
    D = zeros((no_dims + 1) * no_analyzers, (no_dims + 1) * no_analyzers);
    for i=1:no_analyzers
        Ds = zeros(no_dims + 1, no_dims + 1);
        for j=1:size(X, 1)
            Ds = Ds + R(i, j) .* (Z(:,i,j) * Z(:,i,j)');
        end
        D((i - 1) * (no_dims + 1) + 1:i * (no_dims + 1), (i - 1) * (no_dims + 1) + 1:i * (no_dims + 1)) = Ds;
    end
    
    % Construct responsibility weighted local representation matrix U
    R = reshape(R, [1 no_analyzers size(X, 1)]);
    U = reshape(bsxfun(@times, R, Z), [kf size(X, 1)])';
    
    % Solve generalized eigenproblem
    if strcmp(eig_impl, 'Matlab')
        options.disp = 0;
        options.isreal = 1;
        options.issym = 1;
        [V, lambda] = eigs(D - U' * U, U' * U, no_dims + 1, 'SM', options);
    else
        options.Disp = 0;
        options.LSolver = 'bicgstab';
        [V, lambda] = jdqz(D - U' * U, U' * U, no_dims + 1, 'SM', options);
    end
    [lambda, ind] = sort(diag(lambda));
    V = V(:,ind(2:end));
    
    % Compute final lowdimensional data representation
    mappedX = U * V;
    
    % Store mapping data for out-of-sample extension
    if nargout > 1
        mapping.LX = LX;
        mapping.MX = MX;
        mapping.PX = PX;
        mapping.V = V;
    end
    