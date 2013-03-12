function [R, Z] = infermfa(X, LX, MX, PX)
%INFERMFA Infer MoFA using information from EM algorithm in MPPCA
% 
%   [R, Z] = infermfa(X, LX, MX, PX)
%
% Computes local data representations and responsibilities of the datapoints 
% to the clusters specified by LX, MX, and PX. Basically, this method performs
% an additional E-step of the EM-algorithm executed in MPPCA.
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
    [D N] = size(X);
    [D no_dims no_analyzers] = size(LX);
    const = -D / 2 * log(2 * pi);
    R  = zeros(no_analyzers, N);
    Z  = zeros(no_dims, N, no_analyzers);

    % Estimate the cluster centers based on input (= E-step)
    pii = 1 ./ PX;
    for kk=1:no_analyzers
        l        = LX(:,:,kk);
        ltpi     = (repmat(pii, [1 no_dims]) .* l)';
        ltpil    = ltpi * l;
        iltpil   = eye(no_dims) + ltpil;
        cc       = chol(iltpil);
        cci      = inv(cc);
        covz     = cci * cci';
        delta    = X - MX(:,kk * ones(1, N));
        meanz    = ((eye(no_dims) - ltpil * covz) * ltpi) * delta;

        Z(:,:,kk)  = meanz;
        R(kk,:)    = -.5 * (pii' * (delta .* delta) - sum(meanz .* (iltpil * meanz), 1)) - ...
                      sum(log(diag(cc)));
    end

    % Compute responsibilities of clusters to points
    R = R + const + .5 * sum(log(pii));
    R = exp(R - repmat(max(R, [], 1), [no_analyzers 1]));
    R = R ./ repmat(sum(R, 1), [no_analyzers 1]);

