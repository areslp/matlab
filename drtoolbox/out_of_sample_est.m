function t_points = out_of_sample_est(points, X, mappedX)
%TRANSFORM_SAMPLE_EST Performs out-of-sample extension using estimation technique
%
%   t_points = out_of_sample_est(points, X, mappedX)
%
% Performs out-of-sample extension using estimation technique on datapoints
% points. You also need to specify the original dataset in X, and the
% reduced dataset in mappedX (the two datasets may also be PRTools datasets).
% The function returns the coordinates of the transformed points in t_points.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    welcome;

    % Handle PRTools dataset
    if strcmp(class(points), 'dataset')
        prtools = 1;
        ppoints = points;
        points = points.data;
    else 
        prtools = 0;
    end

    % Handle PRTools datasets
    if strcmp(class(X), 'dataset')
        X = X.data;
    end
    if strcmp(class(mappedX), 'dataset')
        mappedX = mappedX.data;
    end
    
    % Remove duplicates from the dataset
    X = double(unique(X, 'rows'));

    % For all datapoints
    t_points = repmat(0, [size(points, 1) size(mappedX, 2)]);
    bb = sum(X' .* X');
    for i=1:size(points, 1)
        
        % Get current point
        point = points(i,:);
        
        % Find nearest neighbor for current point
        n = size(X, 1);
        aa = sum(point .* point);
        ab = point * X';
        d = sqrt(repmat(aa', [1 size(bb, 2)]) + repmat(bb, [size(aa, 2) 1]) - 2 * ab);
        [d, ind] = min(d);

        % Compute transformation matrix
        L = pinv(X(ind,:) - mean(X(ind,:))) * (mappedX(ind,:) - mean(mappedX(ind,:)));

        % Compute coordinates of transformed point
        t_points(i,:) = (mean(mappedX(ind,:)) + ((point - mean(X(ind,:))) * L))';
    end
    
    % Handle PRTools dataset
    if prtools == 1
        ppoints.data = t_points;
        t_points = ppoints;
    end
    