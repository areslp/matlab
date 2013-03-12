function t_point = out_of_sample(point, mapping)
%TRANSFORM_SAMPLE_EST Performs out-of-sample extension of new datapoints
%
%   t_points = out_of_sample(points, mapping)
%
% Performs out-of-sample extension of the new datapoints in points. The 
% information necessary for the out-of-sample extension is contained in 
% mapping (this struct can be obtained from COMPUTE_MAPPING).
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
    if strcmp(class(point), 'dataset')
        prtools = 1;
        ppoint = point;
        point = point.data;
    else 
        prtools = 0;
    end

    switch mapping.name
        
        % Linear mappings
        case {'PCA', 'LDA', 'LPP', 'NPE', 'LLTSA', 'SPCA', 'PPCA', 'FA', 'NCA', 'MCML', 'LMNN'}
            t_point = bsxfun(@minus, point, mapping.mean) * mapping.M;
            
        % Kernel PCA mapping
        case 'KernelPCA'
            
            % Compute and center kernel matrix
            K = gram(mapping.X, point, mapping.kernel, mapping.param1, mapping.param2);
            J = repmat(mapping.column_sums', [1 size(K, 2)]);
            K = K - repmat(sum(K, 1), [size(K, 1) 1]) - J + repmat(mapping.total_sum, [size(K, 1) size(K, 2)]);
            
            % Compute transformed points
            t_point = mapping.invsqrtL * mapping.V' * K;
            t_point = t_point';            
            
        case 'Autoencoder'
            [foo, t_point] = run_data_through_autoenc(mapping.network, point);
            
        case {'Isomap', 'LandmarkIsomap'}
            
            % Precomputations for speed
            if strcmp(mapping.name, 'Isomap')
                invVal = inv(diag(mapping.val));
                [val, index] = sort(mapping.val, 'descend');
                mapping.landmarks = 1:size(mapping.X, 1);
            else
                val = mapping.beta .^ (1 / 2);
                [val, index] = sort(real(diag(val)), 'descend');
            end
            val = val(1:mapping.no_dims);
            meanD1 = mean(mapping.DD .^ 2, 1);
            meanD2 = mean(mean(mapping.DD .^ 2));
            
            % Process all points (notice that in this implementation 
            % out-of-sample points are not used as landmark points)
            points = point;
            t_point = repmat(0, [size(point, 1) mapping.no_dims]);
            for i=1:size(points, 1)
                
                % Compute distance of new sample to training points
                point = points(i,:);
                tD = L2_distance(point', mapping.X');
                [tmp, ind] = sort(tD); 
                tD(ind(mapping.k + 2:end)) = 0;
                tD = sparse(tD);
                tD = dijkstra([0 tD; tD' mapping.D], 1);
                tD = tD(mapping.landmarks + 1) .^ 2;

                % Compute point embedding
                subB = -.5 * (bsxfun(@minus, tD, mean(tD, 2)) - meanD1 - meanD2);
                if strcmp(mapping.name, 'LandmarkIsomap')
                    vec = subB * mapping.alpha * mapping.invVal;
                    vec = vec(:,index(1:mapping.no_dims));
                else
                    vec = subB * mapping.vec * invVal;
                    vec = vec(:,index(1:mapping.no_dims));
                end
                t_point(i,:) = real(vec .* sqrt(val)');
            end
            
        case 'LLE'
            % Initialize some variables
            n = size(mapping.X, 1);
            t_point = repmat(0, [size(point, 1) numel(mapping.val)]);
            
            % Compute local Gram matrix
            D = (L2_distance(point', mapping.X') .^ 2);
            [foo, ind] = sort(D, 2, 'ascend');
            for i=1:size(point, 1)

                % Compute local Gram matrix
                C = (repmat(point(i,:), [mapping.k 1]) - mapping.X(ind(i, 2:mapping.k + 1),:)) * ...
                    (repmat(point(i,:), [mapping.k 1]) - mapping.X(ind(i, 2:mapping.k + 1),:))';

                % Compute reconstruction weights
                invC = inv(C);
                W = sum(invC, 2) ./ sum(sum(invC));

                % Compute kernel matrix
                K = repmat(0, [n 1]);
                K(ind(i, 2:mapping.k + 1)) = W;

                % Compute embedded point
                t_point(i,:) = sum(mapping.vec .* repmat(K, [1 size(mapping.vec, 2)]), 1);
            end
        
        case 'Laplacian'         
            % Initialize some other variables
            n = size(mapping.X, 1);
            
            % Compute embeddings
            t_point = repmat(0, [size(point, 1) numel(mapping.val)]);
            for i=1:size(point, 1)
                
                % Compute Gaussian kernel between test point and training points
                K = (L2_distance(point(i,:)', mapping.X') .^ 2)';
                [foo, ind] = sort(K, 'ascend');            
                K(ind(mapping.k+1:end)) = 0;
                K(K ~= 0) = exp(-K(K ~= 0) / (2 * mapping.sigma ^ 2));

                % Normalize kernel
                K = (1 ./ n) .* (K ./ sqrt(mean(K) .* mean(mapping.K, 2)));

                % Compute embedded point
                t_point(i,:) = sum(mapping.vec .* repmat(K, [1 size(mapping.vec, 2)]), 1);
            end
            
        case 'LandmarkMVU'
            % Initialize some variables
            n = size(point, 1);
            
            % Compute pairwise distances
            X = L2_distance(point', mapping.X); 
            neighbors = zeros(mapping.k2, n);
            
            % Compute reconstruction weights
            tol = 1e-7;
            Pia = sparse([], [], [], mapping.no_landmarks, n);
            for i=1:n
                
                % Identify nearest neighbors in distance matrix
                dist = L2_distance(X(:,i), mapping.D);
                [dist, ind] = sort(dist, 'ascend');
                neighbors(:,i) = ind(2:mapping.k2 + 1);
                
                % Compute reconstruction weights
                z = mapping.D(:,neighbors(:,i)) - repmat(X(:,i), 1, mapping.k2);
                C = z' * z;
                C = C + tol * trace(C) * eye(mapping.k2) / mapping.k2;
                invC = inv(C);
                Pia(neighbors(:,i), i) = sum(invC)' / sum(sum(invC));
            end

            % Fill sparse LLE weight matrix
            M = speye(n) + sparse([], [], [], n, n, n * mapping.k2 .^ 2);
            for i=1:n
                j = neighbors(:,i);
                w = Pia(j, i);
                M(i, j) = M(i, j) - w';
                M(j, i) = M(j, i) - w;
                M(j, j) = M(j, j) + w * w';
            end

            % Invert LLE weight matrix
            Pia = -M(mapping.no_landmarks + 1:end, mapping.no_landmarks + 1:end) \ ...
                   M(mapping.no_landmarks + 1:end, 1:mapping.no_landmarks);
            Pia = [eye(mapping.no_landmarks); Pia];
            
            % Apply mapping on the data
            t_point = mapping.Y * Pia';
            t_point = t_point(1:mapping.no_dims,:)';
            
        case 'FastMVU'
            
            if ~mapping.finetune
                % Initialize some other variables
                n = size(mapping.X, 1);   

                % Start with out-of-sample extension for Laplacian Eigenmaps
                Y = repmat(0, [size(point, 1) size(mapping.vec, 2)]);
                for i=1:size(point, 1)

                    % Compute adjecency matrix between test point and training points
                    K = L2_distance(point(i,:)', mapping.X(mapping.conn_comp,:)')' .^ 2;
                    [foo, ind] = sort(K, 'ascend');
                    K(ind(mapping.k + 1:end)) = 0;
                    K(ind(1:mapping.k)) = 1;

                    % Normalize kernel
                    K = (1 ./ n) .* (K ./ sqrt(mean(K) .* mean(mapping.D, 2)));

                    % Compute estimated eigenvectors of graph Laplacian
                    Y(i,:) = sum(mapping.vec .* repmat(K, [1 size(mapping.vec, 2)]), 1);
                end

                % Out-of sample extension to obtain initial solutions
                newY = mapping.L * Y';
                newY = mapping.newV' * newY;
                newY = newY(mapping.idx(end:-1:1),:)';

                % Apply the PCA mapping
                t_point = out_of_sample(newY, mapping.pca_map);
            else
                error('Out-of-sample extension for FastMVU is only available when finetuning was disabled.');
            end
            
        case 'ManifoldChart'
            
            % Set some variables
            no_dims = size(mapping.V, 2);
            no_analyzers = size(mapping.LX, 3);
            kf = no_analyzers * (no_dims + 1);
            
            % Infer locations and mixing proportions under MoPPCA model
            [R, Z] = infermfa(point', mapping.LX, mapping.MX, mapping.PX);
            
            % Construct matrix U
            Z(no_dims + 1,:,:) = 1; 
            Z = permute(Z, [1 3 2]);
            R = reshape(R, [1 no_analyzers size(point, 1)]);
            U = reshape(bsxfun(@times, R, Z), [kf size(point, 1)])';
            
            % Apply charting linear mapping
            t_point = U * mapping.V;
                        
        otherwise
            error(['An out-of-sample extension for ' mapping.name ' is not available in the toolbox. You might consider using OUT_OF_SAMPLE_EST instead.']);
    end
    
    % Handle PRTools dataset
    if prtools == 1
        ppoint.data = t_point;
        t_point = ppoint;
    end
    