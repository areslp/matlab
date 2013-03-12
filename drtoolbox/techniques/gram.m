function G = gram(X1, X2, kernel, param1, param2)
%GRAM Computes the Gram-matrix of data points X using a kernel function
%
%   G = gram(X1, X2, kernel, param1, param2)
%
% Computes the Gram-matrix of data points X1 and X2 using the specified kernel
% function. If no kernel is specified, no kernel function is applied. The
% function GRAM is than equal to X1*X2'. The use of the function is different
% depending on the specified kernel function (because different kernel
% functions require different parameters. The possibilities are listed
% below.
% Linear kernel: G = gram(X1, X2, 'linear')
%           which is parameterless
% Gaussian kernel: G = gram(X1, X2, 'gauss', s)
%           where s is the variance of the used Gaussian function (default = 1).
% Polynomial kernel: G = gram(X1, X2, 'poly', R, d)
%           where R is the addition value and d the power number (default = 0 and 3)
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology

    % Check inputs
    if size(X1, 2) ~= size(X2, 2)
        error('Dimensionality of both datasets should be equal');
    end

    % If no kernel function is specified
    if nargin == 2 || strcmp(kernel, 'none')
        kernel = 'linear';
    end
    
    switch kernel
        
        % Linear kernel
        case 'linear'
            G = X1 * X2';
        
        % Gaussian kernel
        case 'gauss'
            if ~exist('param1', 'var'), param1 = 1; end
            G = L2_distance(X1', X2');
            G = exp(-(G.^2 / (2 * param1.^2)));
                        
        % Polynomial kernel
        case 'poly'
            if ~exist('param1', 'var'), param1 = 1; param2 = 3; end
            G = ((X1 * X2') + param1) .^ param2;
            
        otherwise
            error('Unknown kernel function.');
    end
    