function [C, dC] = mvu_x_grad(Y, N, no_dims, ni, K, DX, lambda)
%MVU_X_GRAD This file is not currently used by the toolbox

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology



    % Decode solution
    Y = reshape(Y, [N no_dims]);
    
    % Compute slack
    DY = zeros(N, K);
    sum_Y = sum(Y .^ 2, 2);
    for k=1:K
        DY(:,k) = sum_Y + sum_Y(ni(:,k)) - 2 * sum(Y .* Y(ni(:,k),:), 2);
    end
    slack = DX - DY;
    
    % Compute cost function
    Yn = bsxfun(@minus, Y, mean(Y, 1));
    C = -lambda .* sum(Yn(:) .^ 2) + (1 - lambda) .* sum(abs(slack(:)));
    
    % Compute gradient
    if nargout > 1
    
        % Compute gradient w.r.t. objective
        dC = 2 .* lambda .* Yn;

        % Add gradient w.r.t. slack
        for k=1:K
            grad = bsxfun(@times, 2 .* (1 - lambda) .* sign(slack(:,k)), Y - Y(ni(:,k),:));
            dC = dC - grad;
            for d=1:no_dims
                dC(:,d) = dC(:,d) + accumarray(ni(:,k), grad(:,d), [N 1]);
            end
        end
        dC = dC(:);
    end
    