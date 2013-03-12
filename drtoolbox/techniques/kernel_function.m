function y = kernel_function(v, X, center, kernel, param1, param2, type)
%KERNEL_FUNCTION Computes sum of (K * X) where X is a possible eigenvector
%
%   y = kernel_function(v, X, center, kernel, param1, param2)
%
% The function computes the sum of the elements of (K * v), where v is a
% possible eigenvector of K. This function is used to enable the use of
% EIGS in Kernel PCA. The other parameters of the function are the dataset 
% X, the name of the kernel function (default = 'gauss'), and its 
% corresponding parameters in param1 and param2.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    if ~exist('center', 'var')
        center = 0;
    end
    if ~exist('type', 'var')
        type = 'Normal';
    end
    if ~strcmp(type, 'ColumnSums'), fprintf('.'); end    
        
    % If no kernel function is specified
    if nargin == 2 || strcmp(kernel, 'none')
        kernel = 'linear';
    end
    
    % Construct result vector
    y = zeros(1, size(X, 1));
    n = size(X, 2);
    
    switch kernel
        
        case 'linear'
            % Retrieve information for centering of K
            if center || strcmp(type, 'ColumnSums')
                column_sum = zeros(1, n);
                for i=1:n
                    % Compute single row of the kernel matrix
                    K = X(:,i)' * X;
                    column_sum = column_sum + K;
                end
                % Compute centering constant over entire kernel
                total_sum = ((1 / n^2) * sum(column_sum));
            end
            
            if ~strcmp(type, 'ColumnSums')
                % Compute product K*v
                for i=1:n
                    % Compute single row of the kernel matrix
                    K = X(:,i)' * X;

                    % Center row of the kernel matrix
                    if center
                        K = K - ((1 / n) .* column_sum) - ((1 / n) .* column_sum(i)) + total_sum;
                    end

                    % Compute sum of products
                    y(i) = K * v;
                end
            else
                % Return column sums
                y = column_sum;
            end
            
        case 'poly'
            
            % Initialize some variables
            if ~exist('param1', 'var'), param1 = 1; param2 = 3; end            
                        
            % Retrieve information for centering of K
            if center || strcmp(type, 'ColumnSums')
                column_sum = zeros(1, n);
                for i=1:n
                    % Compute column sums of the kernel matrix
                    K = X(:,i)' * X;
                    K = (K + param1) .^ param2;
                    column_sum = column_sum + K;
                end
                % Compute centering constant over entire kernel
                total_sum = ((1 / n^2) * sum(column_sum));
            end       

            if ~strcmp(type, 'ColumnSums')
                % Compute product K*v
                for i=1:n
                    % Compute row of the kernel matrix
                    K = X(:,i)' * X;
                    K = (K + param1) .^ param2;

                    % Center row of the kernel matrix
                    if center
                        K = K - ((1 / n) .* column_sum) - ((1 / n) .* column_sum(i)) + total_sum;
                    end

                    % Compute sum of products
                    y(i) = K * v;
                end
            else
                % Return column sums
                y = column_sum;
            end
            
        case 'gauss'
            
            % Initialize some variables
            if ~exist('param1', 'var'), param1 = 1; end
            
            % Retrieve information for centering of K
            if center || strcmp(type, 'ColumnSums')
                column_sum = zeros(1, n);
                for i=1:n
                    % Compute row sums of the kernel matrix
                    K = L2_distance(X(:,i), X);
                    K = exp(-(K.^2 / (2 * param1.^2)));
                    column_sum = column_sum + K;
                end
                % Compute centering constant over entire kernel
                total_sum = ((1 / n^2) * sum(column_sum));
            end

            if ~strcmp(type, 'ColumnSums')
                % Compute product K*v
                for i=1:n
                    % Compute single row of the kernel matrix
                    K = L2_distance(X(:,i), X);
                    K = exp(-(K.^2 / (2 * param1.^2)));

                    % Center row of the kernel matrix                    
                    if center
                        K = K - ((1 / n) .* column_sum) - ((1 / n) .* column_sum(i)) + total_sum;
                    end
                    
                    % Compute sum of products
                    y(i) = K * v;
                end
            else
                % Return column sums
                y = column_sum;
            end
            
        otherwise
            error('Unknown kernel function.');
    end