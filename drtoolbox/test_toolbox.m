function test_toolbox
%TEST_TOOLBOX Tests all functionalities of the dimension reduction toolbox
%
%   test_toolbox
%
% Tests all functionalities of the dimension reduction toolbox.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    % Generate data
    disp('Testing data generation functions...');
    datasets = {'helix', 'twinpeaks', '3d_clusters', 'intersect', 'swiss'};
    for i=1:length(datasets)
        try        
            X = generate_data(datasets{i}, 500);
        catch e
            disp(e);
            warning(['Generation of data set ' datasets{i} ' failed! Press any key to continue tests...']);
            pause
        end
    end
    
    % Test prewhitening
    disp('Testing prewhitening...');
    try 
        X = prewhiten(X);
    catch e
        disp(e);
        warning('Prewhitening failed! Press any key to continue tests...');
        pause
    end
    unscaled_X = X;
    X = X - min(X(:));
    X = X / max(X(:));
    
    % Test all intrinsic dimensionality estimators
    disp('Testing intrinsic dimensionality estimators...');
    techniques = {'CorrDim', 'NearNbDim', 'GMST', 'PackingNumbers', 'EigValue', 'MLE'};
    for i=1:length(techniques)
        try
            intrinsic_dim(X, techniques{i});
        catch e
            disp(e);
            warning(['Intrinsic dimensionality estimation using ' techniques{i} ' failed! Press any key to continue tests...']);
            pause
        end
    end
    
    % Test all unsupervised dimension reduction techniques
    disp('Testing dimensionality reduction techniques...');
    techniques = {'PCA', 'MDS', 'ProbPCA', 'FactorAnalysis', 'GPLVM', 'Sammon', 'Isomap', ...
        'LandmarkIsomap', 'LLE', 'Laplacian', 'HessianLLE', 'LTSA', 'MVU', 'CCA', 'LandmarkMVU', ...
        'FastMVU', 'DiffusionMaps', 'KernelPCA', 'GDA', 'SNE', 'SymSNE', 'tSNE', 'LPP', 'NPE', ...
        'LLTSA', 'SPE', 'Autoencoder', 'LLC', 'ManifoldChart', 'CFA'};
    for i=1:length(techniques)
        
        % Test the dimension reduction technique
        try
            if any(strcmpi(techniques{i}, {'GPLVM', 'CFA'}))
                [mappedX, mapping] = compute_mapping(unscaled_X, techniques{i}, 2);
            else
                [mappedX, mapping] = compute_mapping(X, techniques{i}, 2);
            end
            if any(strcmpi(techniques{i}, {'Isomap', 'LandmarkIsomap', 'LLE', 'Laplacian', 'MVU', 'CCA', 'FastMVU', 'LPP', 'NPE', 'LLTSA'}))
                [mappedX, mapping] = compute_mapping(X, techniques{i}, 2, 'adaptive');
            end
        catch e
            disp(e);
            warning(['Technique ' techniques{i} ' failed! Press any key to continue tests...']);
            pause
        end
        
        % Test the out-of-sample extension code
        if any(strcmpi(techniques{i}, {'PCA', 'LPP', 'NPE', 'LLTSA', 'SPCA', 'PPCA', 'FA'}))
            try
                out_of_sample(X, mapping);
            catch e
                disp(e);
                warning(['Out-of-sample extension for technique ' techniques{i} ' failed! Press any key to continue tests...']);
                pause                
            end
        end
        
        % Test reconstruction code
        if any(strcmpi(techniques{i}, {'PCA', 'LPP', 'NPE', 'LLTSA', 'SPCA', 'PPCA', 'FA', 'Autoencoder'}))
            try
                reconstruct_data(mappedX, mapping);
            catch e
                disp(e);
                warning(['Reconstruction for technique ' techniques{i} ' failed! Press any key to continue tests...']);
                pause
            end
        end
    end
    
    % Test approximate out-of-sample function
    try
        out_of_sample_est(X, X, mappedX);
    catch e
        disp(e);
        warning(['Approximate out-of-sample extension failed! Press any key to continue tests...']);
        pause                
    end
    
    % Test all supervised dimension reduction techniques
    labels = double(X(:,1) > .5) + 1;
    X = [labels X];
    techniques = {'LDA', 'NCA', 'MCML', 'LMNN'};
    for i=1:length(techniques)
        
        % Test the actual technique
        try
            [mappedX, mapping] = compute_mapping(X, techniques{i}, 2);
        catch e
            disp(e);
            warning(['Technique ' techniques{i} ' failed! Press any key to continue tests...']);
            pause
        end
        
        % Test out-of-sample extension
        try
            out_of_sample(X(:,2:end), mapping);
        catch e
            disp(e);
            warning(['Out-of-sample extension for technique ' techniques{i} ' failed! Press any key to continue tests...']);
            pause
        end
        
        % Test reconstruction code
        try
            reconstruct_data(mappedX, mapping);
        catch e
            disp(e);
            warning(['Reconstruction for technique ' techniques{i} ' failed! Press any key to continue tests...']);
            pause
        end
        
    end
    disp('All tests completed!');    
    