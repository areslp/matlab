function Y = mvu_x(X, no_dims, K, lambda)
%MVU_X This file is not currently used by the toolbox

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    if ~exist('lambda', 'var') || isempty(lambda)
        lambda = .5;
    end

    % Build nearest neighborhood graph and compute distances
    disp('Computing nearest-neighbor graph...');
    N = size(X, 1);
    [D, ni] = find_nn(X, K);
    DX = zeros(N, K);
    for n=1:N
        DX(n,:) = full(D(n, ni(n,:)));
    end
    DX = DX .^ 2;
        
    % Initialize solution
    init_Y = compute_mapping(X, 'PCA', no_dims);
    Y = cell(length(lambda), 1);
    
    % Loop over all values of lambda
%     checkgrad('mvu_x_grad', init_Y(:), 1e-7, N, no_dims, ni, K, DX, lambda)
    for i=1:length(lambda)
        
        % Perform learning
        Y{i} = init_Y;
        options.Method = 'lbfgs';
        options.MaxIter = 5000;
        options.MaxFunEvals = 10000;
        options.Display = 'on';
        Y{i} = reshape(minFunc(@mvu_x_grad, Y{i}(:), options, N, no_dims, ni, K, DX, lambda(i)), [N no_dims]);

        % Make sure solution is zero-mean
        Y{i} = bsxfun(@minus, Y{i}, mean(Y{i}, 1));
    end
end

