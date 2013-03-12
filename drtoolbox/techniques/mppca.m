function [LX, MX, PX] = mppca(X, no_dims, no_analyzers, tol, maxiter, minstd)
%MPPCA Runs EM algorithm and computes local factor analyzers
%
%   [LX, MX, PX] = mppca(X, no_dims, no_analyzers, tol, maxiter, minstd)
%
% Runs EM algorithm to determine coordinates of factor analyzers. The data
% is given in the DxN matrix X. no_dims indicates the number of dimensions that
% the local factor analyzers compute their embedding in. The number of 
% factor analyzers that is used is given by no_analyzers. The variable tol indicates 
% the tolreance in considering the EM as converged, whereas maxiter
% indicates the maximum number of iterations for the EM algorithm. The
% variable minstd sets the minimum standard deviation of the Gaussians that
% the EM algorithm fits.
% The function returns in LX the lowdimensional representations of X of all
% factor analyzers, in MX the means of the factors analyzers, and in PX the
% noise covariances.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    % Initialize some variables
    [D N] = size(X);                % size of data
    epsilon = 1e-9;                 % regularization parameter
    
    % Estimate minimum variance allowed
    minvar = minstd ^ 2;

    % Randomly initialize factor analyzers
    mm = mean(X, 2);
    ss = cov(X');
    try                                                 % for small problems        
        cc = chol(ss);
        MX = bsxfun(@plus, cc' * randn(D, no_analyzers), mm);
        LX = minstd * randn(D, no_dims, no_analyzers);
        PX = 2 * mean(diag(cc)) * ones(D, 1);
    catch                                               % for large problems (or nearly singular covariance matrices)
        cc = std(X, [], 2);
        MX = bsxfun(@plus, bsxfun(@times, cc, randn(D, no_analyzers)), mm);
        LX = minstd * randn(D, no_dims, no_analyzers);
        PX = 2 * mean(cc) * ones(D, 1);
    end
    clear mm ss cc

    % Compute squared data
    X2 = X .^ 2;
    
    % Initialize some variables
    const = -D / 2 * log(2 * pi);
    lik = -Inf;
    R  = zeros(no_analyzers, N);
    czz = zeros(no_dims, no_dims, no_analyzers);
    zz  = zeros(no_dims, N, no_analyzers);

    % Run for maxiter iterations at max
    for i=1:maxiter

        % Progress bar
        if rem(i, 10) == 0
            fprintf('.');
        end

        % E step
        pii = 1 ./ PX;
        for k=1:no_analyzers
            l        = LX(:,:,k);                                                                       % select k-th local representation
            ltpi     = bsxfun(@times, pii, l)';
            ltpil    = ltpi * l;
            iltpil   = eye(no_dims) + ltpil;
            cc       = chol(iltpil);
            cci      = inv(cc);
            covz     = cci * cci';                                                                      % compute covariance
            delta    = X - MX(:,k * ones(1, N));
            meanz    = ((eye(no_dims) - ltpil * covz) * ltpi) * delta;
            czz(:,:,k) = covz;
            zz(:,:,k)  = meanz;                                                                         % update local representation
            R(k,:)    = -.5 * (pii' * (delta .* delta) - sum(meanz .* (iltpil * meanz), 1)) - ...       % update reponsibilie
                          sum(log(diag(cc)));
        end

        % Compute responsibilities of datapoints to the clusters
        R = R + const + .5 * sum(log(pi));
        R = exp(bsxfun(@minus, R, max(R, [], 1)));
        R = bsxfun(@rdivide, R, sum(R, 1));

        % Update likelihood of estimation
        oldlik = lik;
        lik = sum(max(R, [], 1) + log(sum(R, 1)));
        
        % Stop EM after convergence
        if abs(oldlik - lik) < tol
            break;
        end

        % M step
        PX = 0;
        for k=1:no_analyzers                         % Update all factor analyzers
            r    = R(k,:);
            z    = zz(:,:,k);
            rz   = bsxfun(@times, z, r);
            sr   = sum(r);
            srz  = sum(rz, 2);
            srxz = X * rz';
            srx  = X * r';
            m1   = [srxz srx];
            m2   = [sr * czz(:,:,k) + z * rz' srz; srz' sr];
            m1   = m1 / (m2 + (rand(size(m2)) * epsilon));     % random regularization to make sure that INV(m2) does not contain Infs
            LX(:,:,k) = m1(:,1:no_dims);
            MX(:,k) = m1(:,no_dims + 1);
            PX = PX + X2 * r' - sum(LX(:,:,k) .* srxz, 2) - MX(:,k) .* srx;
        end
        PX = max(minvar, PX / N);
        PX(:) = mean(PX);
    end

    % Done
    disp(' ');
