function [mappedX, mapping] = lmvu(X, no_dims, K, LL)
%LMVU Performs Landmark MVU on dataset X
%
%   [mappedX, mapping] = lmvu(X, no_dims, k1, k2)
%
% The function performs Landmark MVU on the DxN dataset X. The value of k1
% represents the number of nearest neighbors that is employed in the MVU 
% constraints. The value of k2 represents the number of nearest neighbors
% that is employed to compute the reconstruction weights (for embedding the
% non-landmark points).
% 

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology



    if ~exist('K', 'var')
        K = 3;
    end
    if ~exist('LL', 'var')
        LL = 12;
    end
    
    % Save some data for out-of-sample extension
    if ischar(K)
        error('Adaptive neighborhood selection not supported for landmark MVU.');
    end
    mapping.k1 = K;
    mapping.k2 = LL;
    mapping.X = X;

    % Initialize some variables
    N = size(X, 2);     % number of datapoints        
    B = ceil(0.02 * N); % number of landmark points

    % Set some parameters parameters
    pars.ep = eps;
    pars.fastmode = 1;
    pars.warmup = K * B / N;
    pars.verify = 1;
    pars.angles = 1;
    pars.maxiter = 100;
    pars.noise = 0;
    pars.ignore = 0.1;
    pars.penalty = 1;
    pars.factor = 0.9999;

    % Identify nearest neighbors
    disp('Identifying nearest neighbors...');
    KK = max(LL, K);
    X = L2_distance(X, X);          % memory-intensive: O(n^{2}) !!!
    [foo, neighbors] = find_nn(X, KK);
    neighbors = neighbors';
    
    % Get inverse of reconstruction weight matrix
    Pia = getPia(B, LL, X, neighbors);
    
    % Generate SDP problem
    disp('Generating SDP problem...');
    neighbors = neighbors(1:K,:);
    clear temp index sorted
    nck = nchoosek(1:K + 1, 2);
    AA = zeros(N * K, 2);
    pos3 = 1;
    for i=1:N
        ne = neighbors(:,i);
        nne = [ne; i];
        pairs = nne(nck);
        js = pairs(:,1);
        ks = pairs(:,2);
        AA(pos3:pos3 + length(js) - 1,:) = sort([js ks], 2);
        pos3 = pos3 + length(js);
        if i == B
            AA = unique(AA, 'rows');
            ForceC = size(AA, 1);
        end
        if pos3 > size(AA, 1) && i < N
            AA = unique(AA, 'rows');
            pos3 = size(AA, 1) + 1;
            AA = [AA; zeros(round(N / (N - i) * pos3), 2)];
            fprintf('.');
        end
    end
    AA = unique(AA, 'rows');
    AA = AA(2:end,:);
    clear neighbors ne v2 v3 js ks
    bb = zeros(1, size(AA, 1));
    for i=1:size(AA, 1)
        bb(i) = sum((X(:,AA(i, 1)) - X(:,AA(i, 2))) .^ 2);
    end
    disp(' ');

    % Reduce the number of forced vectors
    ii = (1:ForceC)';
    jj = zeros(1, size(AA, 1));
    jj(ii) = 1; 
    jj = find(jj == 0);
    jj1 = jj(jj <= ForceC);
    jj2 = jj(jj >  ForceC);
    jj2 = jj2(randperm(length(jj2)));
    jj1 = jj1(randperm(length(jj1)));
    corder = [ii; jj1'; jj2'];
    AA = AA(corder,:);
    bb = bb(corder);
    ForceC = length(ii);
    clear temp jj1 jj2 jj ii
    Const = max(round(pars.warmup * size(AA, 1)), ForceC);
    [A, b, AA, bb] = getConstraints(AA, Pia, bb, B, Const, pars);
    Qt = sum(Pia, 1)' * sum(Pia, 1);
    A = [Qt(:)'; A];
    b = [0; b];
    clear K;
    solved = 0;

    % Start SDP iterations
    disp('Perform semi-definite programming...');
    disp('CSDP OUTPUT =============================================================================');
    while solved == 0
        
        % Initialize some variables
        c = -vec(Pia' * Pia);
        flags.s = B;
        flags.l = size(A, 1) - 1;
        A = [[zeros(1,flags.l); speye(flags.l)] A];

        % Set c (employ penalty)
        c = [ones(ForceC, 1) .* max(max(c)); zeros(flags.l - ForceC, 1); c];

        % Launch the CSDP solver
        options.maxiter=pars.maxiter;
        [x, d, z, info] = csdp(A, b, c, flags, options);
        K = mat(x(flags.l + 1:flags.l + flags.s ^ 2));

        % Check whether a solution is reached
        solved = isempty(AA);
        A = A(:,flags.l + 1:end);
        xx = K(:);
        if size(AA, 1)
            Aold = size(A,1);
            total = 0;
            while size(A, 1) - Aold < Const && ~isempty(AA)
                [newA, newb, AA, bb] = getConstraints(AA, Pia, bb, B, Const, pars);
                jj = find(newA * xx - newb > pars.ignore * abs(newb));
                if info == 2
                    jj = 1:size(newA, 1);
                end
                total = total + length(jj);
                A(size(A,1) + 1:size(A,1) + length(jj),:) = newA(jj,:);
                b(length(b) + 1:length(b) + length(jj)) = newb(jj);
            end
            if total == 0
                solved = 1;
            end
        else
            solved=1;
        end
        if solved == 1 && pars.maxiter < 100
            pars.maxiter = 100;
        end
    end
    disp('=========================================================================================');

    % Perform eigendecomposition of kernel matrix to compute Y
    disp('Perform eigendecomposition to obtain low-dimensional data representation...');
    [V, D] = eig(K);
    V = V * sqrt(D);
    Y = (V(:,end:-1:1))';
    mappedX = Y * Pia';

    % Reorder data in original order
    mappedX = mappedX(1:no_dims,:);
    
    % Set some information for the out-of-sample extension
    mapping.Y = Y;
    mapping.D = X;
    mapping.no_landmarks = B;
    mapping.no_dims = no_dims;
    
    

% Function that computes LLE weight matrix
function Q = getPia(B, LL, X, neighbors)

    % Initialize some variables
    N = size(X,2);

    % Compute reconstruction weights
    disp('Computing reconstruction weights...');
    tol = 1e-7;
    Pia = sparse([], [], [], B, N);
    for i=1:N
        z = X(:,neighbors(:,i)) - repmat(X(:,i), 1, LL);
        C = z' * z;
        C = C + tol * trace(C) * eye(LL) / LL;
        invC = inv(C);
        Pia(neighbors(:,i),i) = sum(invC)' / sum(sum(invC));
    end

    % Fill sparse LLE weight matrix
    M = speye(N) + sparse([], [], [], N, N, N * LL .^ 2);
    for i=1:N
        j = neighbors(:,i);
        w = Pia(j, i);
        M(i, j) = M(i, j) - w';
        M(j, i) = M(j, i) - w;
        M(j, j) = M(j, j) + w * w';
    end

    % Invert LLE weight matrix
    disp('Invert reconstruction weight matrix...');
    Q = -M(B + 1:end, B + 1:end) \ M(B + 1:end, 1:B);
    Q = [eye(B); Q];



% Functions that constructs the constraints for the SDP
function [A, b, AAomit, bbomit] = getConstraints(AA, Pia, bb, B, Const, pars)
    
    % Initialize some variables
    pos2 = 0;
    perm = 1:size(AA,1);
    if size(AA, 1) > Const
        AAomit = AA(perm(Const + 1:end),:);
        bbomit = bb(perm(Const + 1:end));
        AA = AA(perm(1:Const),:);
        bb = bb(perm(1:Const));
    else
        AAomit = [];
        bbomit = [];
    end

    % Allocate some memory
    persistent reqmem;
    if isempty(reqmem)
        A2 = zeros(size(AA, 1) * B, 3);
    else
        A2 = zeros(reqmem, 3);
    end

    % Set the constraints
    pos = 0;
    for j=1:size(AA, 1)

        % Evaluate for current row in AA
        ii = AA(j, 1);
        jj = AA(j, 2);
        Q = Pia(ii,:)' * Pia(ii,:) - 2 .* Pia(jj,:)' * Pia(ii,:) + Pia(jj,:)' * Pia(jj,:);
        Q = (Q + Q') ./ 2;
        it = find(abs(Q) > pars.ep .^ 2);

        % Constraint found
        if ~isempty(it)
            
            % Allocate more memory if needed
            pos = pos + 1;
            if pos2 + length(it) > size(A2, 1)             
                A2 = [A2; zeros(ceil((size(AA, 1) - j) / j * size(A2, 1)), 3)];
            end
            
            % Set constraint
            A2(1 + pos2:pos2 + length(it), 1) = ones(length(it), 1) .* pos;
            A2(1 + pos2:pos2 + length(it), 2) = it;
            A2(1 + pos2:pos2 + length(it), 3) = full(Q(it));
            pos2 = pos2 + length(it);
        end
    end

    % Construct sparse constraint matrix
    reqmem = pos2;
    A2 = A2(1:pos2,:);
    A = sparse(A2(:,1), A2(:,2), A2(:,3), size(AA,1), B .^ 2);
    b = bb';



% Function that vectorizes a matrix
function v = vec(M)
    v = M(:);


    
% Function that matrixizes a vector
function M = mat(C)
    r = round(sqrt(size(C, 1)));
    M = reshape(C, r, r);
    
