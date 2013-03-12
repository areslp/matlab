function mappedX = hlle(X, no_dims, k, eig_impl)
%HLLE Runs the Hessian LLE algorithm
%
%   mappedX = hlle(X, no_dims, k, eig_impl)
%
% Runs the Hessian LLE algorithm on dataset X to reduce its dimensionality
% to no_dims. The variable k specifies the number of nearest negihbors that
% is used.
%
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
    if ~exist('k', 'var')
        k = 12;
    end
    if ~exist('eig_impl', 'var')
        eig_impl = 'Matlab';
    end

    % Compute nearest neighbors
    if ischar(k)
        warning('Adaptive neighborhood selection often leads to problems in Hessian LLE.');
    end
    disp('Finding nearest neighbors...');
    [D, nind] = find_nn(X, k);
    max_k = size(nind, 2);
    
    % Size of original data
    n = size(X, 1);
    
    % Extra term count for quadratic form
    dp = no_dims * (no_dims + 1) / 2;
    W = sparse([], [], [], dp * n, n, dp * n * max_k);

    % For all datapoints
    disp('Building Hessian estimator for neighboring points...');
    for i=1:n
        % Center datapoints by substracting their mean
        tmp_ind = nind(i,:);
        tmp_ind = tmp_ind(tmp_ind ~= 0);
        kt = length(tmp_ind);
        thisx = X(tmp_ind,:);
        thisx = (thisx - repmat(mean(thisx, 1), kt, 1))';

        % Compute local coordinates (using SVD)
        [U, D, Vpr] = svd(thisx);
        if size(Vpr, 2) < no_dims
            no_dims = size(Vpr, 2);
            dp = no_dims * (no_dims + 1) / 2;
            warning(['Target dimensionality reduced to ' num2str(no_dims) '...']);
        end
        V = Vpr(:,1:no_dims);
		% Basically, the above is applying PCA to the neighborhood of Xi. The PCA mapping that is found
		% (and that is contained in V) is an approximation for the tangent space at Xi.

        % Build Hessian estimator
        clear Yi; clear Pii;
        ct = 0;
        for mm=1:no_dims
            startp = V(:,mm);
            for nn=1:length(mm:no_dims)
                indles = mm:no_dims;
                Yi(:,ct+nn) = startp .* (V(:,indles(nn)));
            end
            ct = ct + length(mm:no_dims);
        end
        Yi = [repmat(1, kt, 1) V Yi];
		        
        % Gram-Schmidt orthogonalization (works different from Matlab QR function)
        [Yt, Orig] = mgs(Yi);
        Pii = Yt(:,no_dims + 2:end)';
        
        % Double check weights sum to 1
        for j=1:dp
            if sum(Pii(j,:)) > 0.0001
                tpp = Pii(j,:) ./ sum(Pii(j,:)); 
            else
                tpp = Pii(j,:);
            end
            
            % Fill weight matrix
            W((i - 1) * dp + j, tmp_ind) = tpp;
        end
    end

    % The weight matrix W is now entirely filled, perform eigenanalysis of W
    disp('Computing HLLE embedding (eigenanalysis)...');

    % Make sparse matrix that is inproduct of W
    G = W' * W;
	G(isnan(G)) = 0;
    G = sparse(G);
    
    % Clear some memory
    clear X thisx W D nind U D Vpr;

    % Perform eigendecomposition
    tol = 0;
    if strcmp(eig_impl, 'JDQR')
        options.Disp = 0;
        options.LSolver = 'bicgstab';
        [mappedX, eigenvals] = jdqr(G, no_dims + 1, tol, options);
    else
        options.disp = 0;
        options.issym = 1;
        options.isreal = 1;
        [mappedX, eigenvals] = eigs(G, no_dims + 1, tol, options);
    end
    
    % Sort eigenvalues and eigenvectors
    [eigenvals, ind] = sort(diag(eigenvals), 'ascend');
    if size(mappedX, 2) < no_dims + 1, no_dims = size(mappedX, 2) - 1; end
    mappedX = mappedX(:,ind(2:no_dims + 1)); 
   
    % Extract nonzero coordinates
    mappedX = mappedX(:,1:no_dims)' * sqrt(n);
    mappedX = mappedX';

