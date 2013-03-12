function mappedX = cfa(X, no_dims, no_analyzers, no_iterations)
%CFA Performs manifold charting on dataset X 
%
%   mappedX = cfa(X, no_dims, no_analyzers, no_iterations)
%
% Performs manifold charting on dataset X to reduce its dimensionality to
% no_dims dimensions. The variable no_analyzers determines the number of
% local factor analyzers that is used in the mixture of factor analyzers
% (default = 40). The variable no_iterations sets the number of
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
    if ~exist('no_iterations', 'var')
        no_iterations = 200;
    end
    
    % Make sure data is zero-mean, unit-variance
    X = X - repmat(mean(X, 1), [size(X, 1) 1]);
    X = X ./ repmat(var(X, 1), [size(X, 1) 1]);
    
    % Initialize some parameters
    min_var = 1e-5;                     % minimum STD of Gaussians
    X = X';
    [D n] = size(X);
    c = no_analyzers;
    d = no_dims;
   
    % Randomly initialize the parameters of the CFA model
    SigmaN = repmat(eye(d), [1 1 n]);
    Z = rand(d, n) * .1;
    Pi = repmat(1 / c, [1 c]);
    Kappa = rand(d, c) * .1;
    Mu = randn(D, c) * .1;
    SigmaC = repmat(eye(d), [1 1 c]);
    Lambda = rand(D, d, c) * .1;
    Psi = zeros(D, D, c);
    for j=1:c
        tmp = zeros(D, D);
        tmp(1:size(tmp, 2) + 1:end) = rand(D, 1) * .05;
        Psi(:,:,j) = tmp + tmp';
    end
    
    % Perform the EM algorithm for the optimization
    iter = 0;
    while iter < no_iterations

        % E-step
        % ====================================================
        
        % Do some precomputations for speed
        invPsi = zeros(size(Psi));
        detSigmaC = zeros(1, c);
        detPsi = zeros(1, c);
        logPi = zeros(1, c);
        for j=1:c
            invPsi(:,:,j) = inv(Psi(:,:,j) + eye(D));
            detSigmaC(j) = det(SigmaC(:,:,j));
            detPsi(j) = det(Psi(:,:,j));
            logPi(j) = log(Pi(j));
        end
        
        % Compute matrices Epsilon, Vc, and m
        Eps = zeros(n, c);
        Vc = zeros(d, d, c);
        m = zeros(d, c);
        const = ((D + d) / 2) * log(2 * pi);
        for j=1:c
            
            % Precomputations
            Xnc = X - repmat(Mu(:,j), [1 n]);
            Znc = Z - repmat(Kappa(:,j), [1 n]);
            tmpProd = Lambda(:,:,j)' * invPsi(:,:,j) * Lambda(:,:,j);
            
            % Compute Epsilon
            for i=1:n            
                normX = (Xnc(:,i) - Lambda(:,:,j) * Znc(:,i));
                Eps(i, j) = -logPi(j) + const ...
                            + (.5 * log(detSigmaC(j))) + (.5 * detPsi(j)) ...
                            + (.5 * trace(SigmaC(:,:,j) * (SigmaN(:,:,i) + Znc(:,i) * Znc(:,i)'))) ...
                            + (.5 * trace(SigmaN(:,:,i) * tmpProd)) ...
                            + (.5 * normX' * invPsi(:,:,j) * normX);
            end
            
            % Compute Vc and m
            Vc(:,:,j) = inv(SigmaC(:,:,j) + eps * eye(d)) + tmpProd;
            m(:,j) = Kappa(:,j) + (Vc(:,:,j) \ Lambda(:,:,j)') * invPsi(:,:,j) * mean(Xnc, 2);
        end
        
        % Update estimate of Q
        Q = (repmat(sum(exp(-Eps), 2), [1 c]) .^ -1) .* exp(-Eps);
        
        % Update estimate of SigmaN
        for i=1:n
            tmp = zeros(d, d);
            for j=1:c
                tmp = tmp + Q(i, j) * Vc(:,:,j);
            end
            
            % Compute covariance matrix
            SigmaN(:,:,i) = inv(tmp + eps * eye(d));           % code above gave us inv(SigmaN)
        end
        
        % Update estimate of Z (= mappedX)
        for i=1:n
            tmp = zeros(d, 1);
            for j=1:c
                tmp = tmp + (Q(i, j) * Vc(:,:,j) * m(:,j));
            end
            Z(:,i) = SigmaN(:,:,i) * tmp;
        end
        
        % M-step
        % ====================================================
        
        % Update estimate of Pi
        Pi = sum(Q, 1) ./ n;
        
        % Update estimate of Mu and Kappa
        for j=1:c
            tmpQ = Q(:,c) ./ sum(Q(:,c));
            Mu(:,c) = sum(repmat(tmpQ', [D 1]) .* X, 2);
            Kappa(:,c) = sum(repmat(tmpQ', [d 1]) .* Z, 2);
        end
        
        % Update estimate of SigmaC
        for j=1:c
            tmp = 0;
            tmpQ = Q(:,c) ./ sum(Q(:,c));
            for i=1:n
                Znc = Z(:,i) - Kappa(:,j);
                tmp = tmp + (tmpQ(i) * (SigmaN(:,:,i) + Znc * Znc'));
            end
            
            % Enforce some variance
            tmp(1:size(tmp, 1) + 1:end) = max(min_var, tmp(1:size(tmp, 1) + 1:end));
            SigmaC(:,:,j) = tmp;
        end
        
        % Update estimate of Lambda
        for j=1:c
            Sc = zeros(D, d);
            tmpQ = Q(:,j) ./ sum(Q(:,j));
            Xnc = X - repmat(Mu(:,j), [1 n]);
            Znc = Z - repmat(Kappa(:,j), [1 n]);
            for i=1:n
                Sc = Sc + tmpQ(i) * (Xnc(:,i) * Znc(:,i)');
            end
            Lambda(:,:,j) = Sc / (SigmaC(:,:,j) + eps * eye(d));
        end
        
        % Update estimate of Psi
        for j=1:c
            tmpPsi = zeros(D, D);
            tmpQ = Q(:,j) ./ sum(Q(:,j));
            Xnc = X - repmat(Mu(:,j), [1 n]);
            Znc = Z - repmat(Kappa(:,j), [1 n]);
            tmpProd = Lambda(:,:,j) * SigmaN(:,:,i) * Lambda(:,:,j)';
            for i=1:n
                tmpPsi(1:size(tmpPsi, 2) + 1:end) = tmpPsi(1:size(tmpPsi, 2) + 1:end) + ...
                    tmpQ(i) * (((Xnc(:,i) - Lambda(:,:,j) * Znc(:,i)) .^ 2)' + tmpProd(1:size(tmpProd, 2) + 1:end));
            end
            Psi(:,:,c) = tmpPsi;
        end
        
        % Update number of iterations
        iter = iter + 1;
        if rem(iter, 5) == 0
            fprintf('.');
        end
    end
    
    % Transpose to get lowdimensional data representation
    mappedX = Z';
    