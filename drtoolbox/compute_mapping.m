function [mappedA, mapping] = compute_mapping(A, type, no_dims, varargin)
%COMPUTE_MAPPING Performs dimensionality reduction on a dataset A
%
%   mappedA = compute_mapping(A, type)
%   mappedA = compute_mapping(A, type, no_dims)
%   mappedA = compute_mapping(A, type, no_dims, ...)
%
% Performs a technique for dimensionality reduction on the data specified 
% in A, reducing data with a lower dimensionality in mappedA.
% The data on which dimensionality reduction is performed is given in A
% (rows correspond to observations, columns to dimensions). A may also be a
% (labeled or unlabeled) PRTools dataset.
% The type of dimensionality reduction used is specified by type. Possible
% values are 'PCA', 'LDA', 'MDS', 'ProbPCA', 'FactorAnalysis', 'GPLVM', 
% 'Sammon', 'Isomap', 'LandmarkIsomap', 'LLE', 'Laplacian', 'HessianLLE', 
% 'LTSA', 'MVU', 'CCA', 'LandmarkMVU', 'FastMVU', 'DiffusionMaps', 
% 'KernelPCA', 'GDA', 'SNE', 'SymSNE', 'tSNE', 'LPP', 'NPE', 'LLTSA', 
% 'SPE', 'Autoencoder', 'LLC', 'ManifoldChart', 'CFA', 'NCA', 'MCML', and 'LMNN'. 
% The function returns the low-dimensional representation of the data in the 
% matrix mappedA. If A was a PRTools dataset, then mappedA is a PRTools 
% dataset as well. For some techniques, information on the mapping is 
% returned in the struct mapping.
% The variable no_dims specifies the number of dimensions in the embedded
% space (default = 2). For the supervised techniques ('LDA', 'GDA', 'NCA', 
% 'MCML', and 'LMNN'), the labels of the instances should be specified in 
% the first column of A (using numeric labels). 
%
%   mappedA = compute_mapping(A, type, no_dims, parameters)
%   mappedA = compute_mapping(A, type, no_dims, parameters, eig_impl)
%
% Free parameters of the techniques can be defined as well (on the place of
% the dots). These parameters differ per technique, and are listed below.
% For techniques that perform spectral analysis of a sparse matrix, one can 
% also specify in eig_impl the eigenanalysis implementation that is used. 
% Possible values are 'Matlab' and 'JDQR' (default = 'Matlab'). We advice
% to use the 'Matlab' for datasets of with 10,000 or less datapoints; 
% for larger problems the 'JDQR' might prove to be more fruitful. 
% The free parameters for the techniques are listed below (the parameters 
% should be provided in this order):
%
%   PCA:            - none
%   LDA:            - none
%   MDS:            - none
%   ProbPCA:        - <int> max_iterations -> default = 200
%   FactorAnalysis: - none
%   GPLVM:          - <double> sigma -> default = 1.0
%   Sammon:         - none
%   Isomap:         - <int> k -> default = 12
%   LandmarkIsomap: - <int> k -> default = 12
%                   - <double> percentage -> default = 0.2
%   LLE:            - <int> k -> default = 12
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   Laplacian:      - <int> k -> default = 12
%                   - <double> sigma -> default = 1.0
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   HessianLLE:     - <int> k -> default = 12
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   LTSA:           - <int> k -> default = 12
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   MVU:            - <int> k -> default = 12
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   CCA:            - <int> k -> default = 12
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   LandmarkMVU:    - <int> k -> default = 5
%   FastMVU:        - <int> k -> default = 5
%                   - <logical> finetune -> default = true
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   DiffusionMaps:  - <double> t -> default = 1.0
%                   - <double> sigma -> default = 1.0
%   KernelPCA:      - <char[]> kernel -> {'linear', 'poly', ['gauss']} 
%                   - kernel parameters: type HELP GRAM for info
%   GDA:            - <char[]> kernel -> {'linear', 'poly', ['gauss']} 
%                   - kernel parameters: type HELP GRAM for info
%   SNE:            - <double> perplexity -> default = 30
%   SymSNE:         - <double> perplexity -> default = 30
%   tSNE:           - <int> initial_dims -> default = 30
%                   - <double> perplexity -> default = 30
%   LPP:            - <int> k -> default = 12
%                   - <double> sigma -> default = 1.0
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   NPE:            - <int> k -> default = 12
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   LLTSA:          - <int> k -> default = 12
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   SPE:            - <char[]> type -> {['Global'], 'Local'}
%                   - if 'Local': <int> k -> default = 12
%   Autoencoder:    - <double> lambda -> default = 0
%   LLC:            - <int> k -> default = 12
%                   - <int> no_analyzers -> default = 20
%                   - <int> max_iterations -> default = 200
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   ManifoldChart:  - <int> no_analyzers -> default = 40
%                   - <int> max_iterations -> default = 200
%                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
%   CFA:            - <int> no_analyzers -> default = 2
%                   - <int> max_iterations -> default = 200
%   NCA:            - <double> lambda -> default = 0.0
%   MCML:           - none
%   LMNN:           - <int> k -> default = 3
%
%
% In the parameter list above, {.., ..} indicates a list of options, and []
% indicates the default setting. The variable k indicates the number of 
% nearest neighbors in a neighborhood graph. Alternatively, k may also have 
% the value 'adaptive', indicating the use of adaptive neighborhood selection
% in the construction of the neighborhood graph. Note that in LTSA and
% HessianLLE, the setting 'adaptive' might cause singularities. Using the
% JDQR-solver or a fixed setting of k might resolve this problem. SPE does
% not yet support adaptive neighborhood selection.
% 
% The variable sigma indicates the variance of a Gaussian kernel. The 
% parameters no_analyzers and max_iterations indicate repectively the number
% of factor analyzers that is used in an MFA model and the number of 
% iterations that is used in an EM algorithm. 
%
% The variable lambda represents an L2-regularization parameter.


% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    welcome;
    
    % Check inputs
    if nargin < 2
        error('Function requires at least two inputs.');
    end
    if ~exist('no_dims', 'var')
        no_dims = 2;
    end
    if ~isempty(varargin) && strcmp(varargin{length(varargin)}, 'JDQR')
        eig_impl = 'JDQR';
        varargin(length(varargin)) = [];
    elseif ~isempty(varargin) && strcmp(varargin{length(varargin)}, 'Matlab')
        eig_impl = 'Matlab';
        varargin(length(varargin)) = [];
    else
        eig_impl = 'Matlab';
    end        
    mapping = struct;
    
    % Handle PRTools dataset
    if strcmp(class(A), 'dataset')
        prtools = 1;
        AA = A;
        if ~strcmp(type, {'LDA', 'FDA', 'GDA', 'KernelLDA', 'KernelFDA', 'MCML', 'NCA', 'LMNN'})
            A = A.data;
        else
            A = [double(A.labels) A.data];
        end
    else 
        prtools = 0;
    end
    
    % Make sure there are no duplicates in the dataset
    A = double(A);
%     if size(A, 1) ~= size(unique(A, 'rows'), 1)
%         error('Please remove duplicates from the dataset first.');
%     end
    
    % Check whether value of no_dims is correct
    if ~isnumeric(no_dims) || no_dims > size(A, 2) || ((no_dims < 1 || round(no_dims) ~= no_dims) && ~any(strcmpi(type, {'PCA', 'KLM'})))
        error('Value of no_dims should be a positive integer smaller than the original data dimensionality.');
    end
    
    % Switch case
    switch type
        case 'Isomap'         
            % Compute Isomap mapping
			if isempty(varargin), [mappedA, mapping] = isomap(A, no_dims, 12);
            else [mappedA, mapping] = isomap(A, no_dims, varargin{1}); end
            mapping.name = 'Isomap';
			
		case 'LandmarkIsomap'
			% Compute Landmark Isomap mapping
            if isempty(varargin), [mappedA, mapping] = landmark_isomap(A, no_dims, 12, 0.2);
			elseif length(varargin) == 1, [mappedA, mapping] = landmark_isomap(A, no_dims, varargin{1}, 0.2);
            elseif length(varargin) >  1, [mappedA, mapping] = landmark_isomap(A, no_dims, varargin{1}, varargin{2}); end
            mapping.name = 'LandmarkIsomap';
            
        case {'Laplacian', 'LaplacianEig', 'LaplacianEigen' 'LaplacianEigenmaps'}
            % Compute Laplacian Eigenmaps-based mapping
            if isempty(varargin), [mappedA, mapping] = laplacian_eigen(A, no_dims, 12, 1, eig_impl);
			elseif length(varargin) == 1, [mappedA, mapping] = laplacian_eigen(A, no_dims, varargin{1}, 1, eig_impl);
            elseif length(varargin) > 1,  [mappedA, mapping] = laplacian_eigen(A, no_dims, varargin{1}, varargin{2}, eig_impl); end
            mapping.name = 'Laplacian';
            
        case {'HLLE', 'HessianLLE'}
            % Compute Hessian LLE mapping
			if isempty(varargin), mappedA = hlle(A, no_dims, 12, eig_impl);
            else mappedA = hlle(A, no_dims, varargin{1}, eig_impl); end
            mapping.name = 'HLLE';
            
        case 'LLE'
            % Compute LLE mapping
			if isempty(varargin), [mappedA, mapping] = lle(A, no_dims, 12, eig_impl);
            else [mappedA, mapping] = lle(A, no_dims, varargin{1}, eig_impl); end
            mapping.name = 'LLE';
            
        case 'GPLVM'
            % Compute GPLVM mapping
			if isempty(varargin), mappedA = gplvm(A, no_dims, 1);
            else mappedA = gplvm(A, no_dims, varargin{1}); end
            mapping.name = 'GPLVM';
            
        case 'LLC'
            % Compute LLC mapping
			if isempty(varargin), mappedA = run_llc(A', no_dims, 12, 20, 200, eig_impl);
            elseif length(varargin) == 1, mappedA = run_llc(A', no_dims, varargin{1}, 20, 200, eig_impl);
            elseif length(varargin) == 2, mappedA = run_llc(A', no_dims, varargin{1}, varargin{2}, 200, eig_impl);
            else mappedA = run_llc(A', no_dims, varargin{1}, varargin{2}, varargin{3}, eig_impl); end
            mappedA = mappedA';
            mapping.name = 'LLC';
                                    
        case {'ManifoldChart', 'ManifoldCharting', 'Charting', 'Chart'}
            % Compute mapping using manifold charting
            if isempty(varargin), [mappedA, mapping] = charting(A, no_dims, 40, 200, eig_impl);
            elseif length(varargin) == 1, [mappedA, mapping] = charting(A, no_dims, varargin{1}, 200, eig_impl);   
            else [mappedA, mapping] = charting(A, no_dims, varargin{1}, varargin{2}, eig_impl); end
            mapping.name = 'ManifoldChart';
            
        case 'CFA'
            % Compute mapping using Coordinated Factor Analysis
            if isempty(varargin), mappedA = cfa(A, no_dims, 2, 200);
            elseif length(varargin) == 1, mappedA = cfa(A, no_dims, varargin{1}, 200);
            else mappedA = cfa(A, no_dims, varargin{1}, varargin{2}); end
            mapping.name = 'CFA';
           
        case 'LTSA'
            % Compute LTSA mapping 
            if isempty(varargin), mappedA = ltsa(A, no_dims, 12, eig_impl);
            else mappedA = ltsa(A, no_dims, varargin{1}, eig_impl); end
            mapping.name = 'LTSA';
            
        case 'LLTSA'
            % Compute LLTSA mapping 
            if isempty(varargin), [mappedA, mapping] = lltsa(A, no_dims, 12, eig_impl);
            else [mappedA, mapping] = lltsa(A, no_dims, varargin{1}, eig_impl); end
            mapping.name = 'LLTSA';
            
        case {'LMVU', 'LandmarkMVU'}
            % Compute Landmark MVU mapping
            if isempty(varargin), [mappedA, mapping] = lmvu(A', no_dims, 5);
            else [mappedA, mapping] = lmvu(A', no_dims, varargin{1}); end
            mappedA = mappedA';
            mapping.name = 'LandmarkMVU';
            
        case 'FastMVU'
            % Compute MVU mapping
            if isempty(varargin), [mappedA, mapping] = fastmvu(A, no_dims, 12, eig_impl);
            elseif length(varargin) == 1, [mappedA, mapping] = fastmvu(A, no_dims, varargin{1}, true, eig_impl); 
            elseif length(varargin) == 2, [mappedA, mapping] = fastmvu(A, no_dims, varargin{1}, varargin{2}, eig_impl);end
            mapping.name = 'FastMVU';
            
        case {'Conformal', 'ConformalEig', 'ConformalEigen', 'ConformalEigenmaps', 'CCA', 'MVU'}
            % Perform initial LLE (to higher dimensionality)
            disp('Running normal LLE...')
            tmp_dims = min([size(A, 2) 4 * no_dims + 1]);
            if isempty(varargin), [mappedA, mapping] = lle(A, tmp_dims, 12, eig_impl);
            else [mappedA, mapping] = lle(A, tmp_dims, varargin{1}, eig_impl); end
            
            % Now perform the MVU / CCA optimalization            
            if strcmp(type, 'MVU'),
                disp('Running Maximum Variance Unfolding...');
                opts.method = 'MVU';
            else
                disp('Running Conformal Eigenmaps...');
                opts.method = 'CCA';
            end
            disp('CSDP OUTPUT =============================================================================');
            mappedA = cca(A(mapping.conn_comp,:)', mappedA', mapping.nbhd(mapping.conn_comp, mapping.conn_comp)', opts);
            disp('=========================================================================================');
            mappedA = mappedA(1:no_dims,:)';
            
        case {'DM', 'DiffusionMaps'}
            % Compute diffusion maps mapping
			if isempty(varargin), mappedA = diffusion_maps(A, no_dims, 1, 1);
            elseif length(varargin) == 1, mappedA = diffusion_maps(A, no_dims, varargin{1}, 1);
            else mappedA = diffusion_maps(A, no_dims, varargin{1}, varargin{2}); end
            mapping.name = 'DM';
            
        case 'SPE'
            % Compute SPE mapping
            if isempty(varargin), mappedA = spe(A, no_dims, 'Global');
            elseif length(varargin) == 1, mappedA = spe(A, no_dims, varargin{1}); 
            elseif length(varargin) == 2, mappedA = spe(A, no_dims, varargin{1}, varargin{2}); end
            mapping.name = 'SPE';
            
        case 'LPP'
            % Compute LPP mapping
            if isempty(varargin), [mappedA, mapping] = lpp(A, no_dims, 12, 1, eig_impl);
            elseif length(varargin) == 1, [mappedA, mapping] = lpp(A, no_dims, varargin{1}, 1, eig_impl); 
            else [mappedA, mapping] = lpp(A, no_dims, varargin{1}, varargin{2}, eig_impl); end
            mapping.name = 'LPP';
            
        case 'NPE'
            % Compute NPE mapping
            if isempty(varargin), [mappedA, mapping] = npe(A, no_dims, 12, eig_impl);
            else [mappedA, mapping] = npe(A, no_dims, varargin{1}, eig_impl); end
            mapping.name = 'NPE';
            
        case 'SNE'
            % Compute SNE mapping
			if isempty(varargin), mappedA = sne(A, no_dims);
            else mappedA = sne(A, no_dims, varargin{1}); end
            mapping.name = 'SNE';

        case {'SymSNE', 'SymmetricSNE'}
            % Compute Symmetric SNE mapping
			if isempty(varargin), mappedA = sym_sne(A, no_dims);
            elseif length(varargin) == 1, mappedA = sym_sne(A, no_dims, varargin{1});
            else mappedA = sym_sne(A, no_dims, varargin{1}, varargin{2}); end
            mapping.name = 'SymSNE';
            
        case {'tSNE', 't-SNE'}
            % Compute t-SNE mapping
			if isempty(varargin), mappedA = tsne(A, [], no_dims);
            else mappedA = tsne(A, [], no_dims, varargin{1}); end
            mapping.name = 't-SNE';
            
        case {'AutoEncoder', 'Autoencoder'}
            
            % Train deep autoencoder to map data
            layers = [ceil(size(A, 2) * 1.2) + 5 max([ceil(size(A, 2) / 4) no_dims + 2]) + 3 max([ceil(size(A, 2) / 10) no_dims + 1]) no_dims];
            disp(['Network size: ' num2str(layers)]);
%             [mappedA, network, binary_data, mean_X, var_X] = train_autoencoder(A, net_structure);
            if isempty(varargin), [network, mappedA] = train_deep_autoenc(A, layers, 0);
            else [network, mappedA] = train_deep_autoenc(A, layers, varargin{1}); end
            mapping.network = network;
            mapping.name = 'Autoencoder';
            
        case {'KPCA', 'KernelPCA'}
            % Apply kernel PCA with polynomial kernel
			if isempty(varargin), [mappedA, mapping] = kernel_pca(A, no_dims);
			else [mappedA, mapping] = kernel_pca(A, no_dims, varargin{:}); end
            mapping.name = 'KernelPCA';
			
		case {'KLDA', 'KFDA', 'KernelLDA', 'KernelFDA', 'GDA'}
			% Apply GDA with Gaussian kernel
            if isempty(varargin), mappedA = gda(A(:,2:end), uint8(A(:,1)), no_dims);
            else mappedA = gda(A(:,2:end), uint8(A(:,1)), no_dims, varargin{:}); end
            mapping.name = 'KernelLDA';
            
        case {'LDA', 'FDA'}
            % Run LDA on labeled dataset
            [mappedA, mapping] = lda(A(:,2:end), A(:,1), no_dims);
            mapping.name = 'LDA';
            
        case 'MCML'
            % Run MCML on labeled dataset
            mapping = mcml(A(:,2:end), A(:,1), no_dims);
            mappedA = bsxfun(@minus, A(:,2:end), mapping.mean) * mapping.M;
            mapping.name = 'MCML';
            
        case 'NCA'
            % Run NCA on labeled dataset
            if isempty(varargin), lambda = 0; else lambda = varargin{1}; end
            [mappedA, mapping] = nca(A(:,2:end), A(:,1), no_dims, lambda);
            mapping.name = 'NCA';
            
        case 'MDS'
            % Perform MDS
            mappedA = mds(A, no_dims);
            mapping.name = 'MDS';
            
        case 'Sammon'
            mappedA = sammon(A, no_dims);
            mapping.name = 'Sammon';
            
        case {'PCA', 'KLM'}
            % Compute PCA mapping
			[mappedA, mapping] = pca(A, no_dims);
            mapping.name = 'PCA';
            
        case {'SPCA', 'SimplePCA'}
            % Compute PCA mapping using Hebbian learning approach
			[mappedA, mapping] = spca(A, no_dims);
            mapping.name = 'SPCA';
            
        case {'PPCA', 'ProbPCA', 'EMPCA'}
            % Compute probabilistic PCA mapping using an EM algorithm
			if isempty(varargin), [mappedA, mapping] = em_pca(A, no_dims, 200);
            else [mappedA, mapping] = em_pca(A, no_dims, varargin{1}); end
            mapping.name = 'PPCA';
            
        case {'FactorAnalysis', 'FA'}
            % Compute factor analysis mapping (using an EM algorithm)
            [mappedA, mapping] = fa(A, no_dims);
            mapping.name = 'FA';
            
        case 'LMNN'
            % Perform large-margin nearest neighbor metric learning
            Y = A(:,1); A = A(:,2:end);
            mapping.mean = mean(A, 1);
            A = bsxfun(@minus, A, mapping.mean);
            [foo, mapping.M, mappedA] = lmnn(A, Y);
            mapping.name = 'LMNN';
            
        otherwise
            error('Unknown dimensionality reduction technique.');
    end
    
    % JDQR makes empty figure; close it
    if strcmp(eig_impl, 'JDQR')
        close(gcf);
    end
    
    % Handle PRTools dataset
    if prtools == 1
        if sum(strcmp(type, {'Isomap', 'LandmarkIsomap', 'FastMVU'}))
            AA = AA(mapping.conn_comp,:);
        end
        AA.data = mappedA;
        mappedA = AA;
    end
    