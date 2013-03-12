function [Z, ccaEigen, ccaDetails] = cca(X, Y, EDGES, OPTS)
%
% Function [Z, CCAEIGEN, CCADETAILS] = CCA(X, Y, EDGES, OPTS) computes a low
% dimensional embedding Z in R^d that maximally preserves angles among  input 
% data X that lives in R^D, with the algorithm Conformal Component Analysis.
%
% The embedding Z is constrained to be Z = L*Y where Y is a partial basis that 
% spans the space of R^d. Such Y can be computed from graph Laplacian (such as 
% the outputs of Laplacian eigenmap and Locally Linear Embedding, ie, LLE). 
% The parameterization matrix L is found by this function as to maximally
% prserve angles between edges coded in the sparse matrix EDGES.
%
% A basic usage of this function is given below:
%
% Inputs:
%   X: input data stored in matrix  (D x N) where D is the dimensionality
%
%   Y: partial basis stored in matrix (d x N)
%
%   EDGES: a sparse matrix of (N x N). In each column i, the row indices j to
%   nonzero entrices define data points that are in the nearest neighbors of
%   data point i.
%
%   OPTS:
%     OPTS.method: 'CCA' 
%
% Outputs:
%   Z: low dimensional embedding (d X N)
%   CCAEIGN: eigenspectra of the matrix P = L'*L. If P is low-rank (say d' < d),
%   then Z can be cutoff at d' dimension as dimensionality reduced further.
%
% The CCA() function is fairly versatile. For more details, consult the file
% README.
%
%  by feisha@cs.berkeley.edu Aug 18, 2006
%  Feel free to use it for educational and research purpose.

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    % sanity check
    if nargin ~= 4
        error('Incorrect number of inputs supplied to cca().');
    end
    N = size(X,2);
    if (N~=size(Y,2)) || (N ~= size(EDGES,1)) || (N~=size(EDGES,2))
        disp('Unmatched matrix dimensions in cca().');
        fprintf('# of data points: %d\n', N);
        fprintf('# of data points in Y: %d\n', size(Y,2));
        fprintf('Size of the sparse matrix for edges: %d x %d\n', size(EDGES,1), size(EDGES,2));
        error('All above 4 numbers should be the same.');
    end
    % check necessary programs
    if exist('mexCCACollectData') ~= 3
        error('Missing mexCCACollectData mex file on the path');
    end
    if exist('csdp') ~= 2
        error('You will need CSDP solver to run cca(). Please make sure csdp.m is in your path');
    end
    % check options
    OPTS = check_opt(OPTS);

    D = size(X, 1); 
    d = size(Y, 1);

    %disp('Step I.  collect data needed for SDP formulation');
    [tnn, vidx] = triangNN(EDGES, OPTS.CCA);
    [erow, ecol, evalue] = sparse_nn(tnn);

    irow = int32(erow); icol = int32(ecol);
    ividx = int32(vidx); ivalue = int32(evalue);
    [A,B, g] = mexCCACollectData(X,Y, irow, icol, int32(OPTS.relative), ivalue, ividx );
    clear erow ecol irow icol tnn ividx ivalue evalue vidx;
    lst = find(g~=0);
    g = g(lst); B = B(:, lst);
    if OPTS.CCA == 1
        BG = B*spdiags(1./sqrt(g),0, length(g),length(g));
        Q = A - BG*BG';
        BIAS = OPTS.regularizer*reshape(eye(d), d^2,1);
    else
        Q = A; BIAS = 2*sum(B,2)+OPTS.regularizer*reshape(eye(d), d^2,1);
    end
    [V, E] = eig(Q+eye(size(Q))); % adding an identity matrix to Q for numerical
    E = E-eye(size(Q));           % stability
    E(E<0) = 0;
    if ~isreal(diag(E))
        warning('\tThe quadratic matrix is not positive definite..forced to be positive definite...\n');
        E=real(E);
        V = real(V);
        S = sqrt(E)*V';
    else
        S = sqrt(E)*V';
    end

    % Formulate the SDP problem
    [AA, bb, cc] = formulateSDP(S, d, BIAS, (OPTS.CCA==1));
    sizeSDP = d^2+1 + d + 2*(OPTS.CCA==1);
    csdppars.s = sizeSDP;
    csdpopts.printlevel = 0;
    
    % Solve it using CSDP
    [xx, yy, zz, info] = csdp(AA, bb, cc, csdppars,csdpopts);
    ccaDetails.sdpflag = info;
    
    % The negate of yy is our solution
    yy = -yy;
    idx = 0;
    P = zeros(d);
    for col=1:d
        for row = col:d
            idx=idx+1;
            P(row, col) = yy(idx);
        end
    end
    
    % Convert P to a positive definite matrix
    P = P + P' - diag(diag(P));

    % Transform the original projection to the new projection
    [V, E] = eig(P);
    E(E < 0) = 0;
    L = diag(sqrt(diag(E))) * V';
    newY = L * Y;

    % Eigenvalue of the new projection, doing PCA using covariance matrix
    [newV, newE] = eig(newY * newY');
    newE = diag(newE);
    [dummy, idx] = sort(newE);
    newE = newE(idx(end:-1:1));
    newY = newV' * newY;
    Z = newY(idx(end:-1:1),:);

    ccaEigen = newE;
    ccaDetails.cost = P(:)'*Q*P(:) - BIAS'*P(:) + sum(g(:))*(OPTS.MVU==1);
    if OPTS.CCA == 1
        ccaDetails.c = spdiags(1./sqrt(g),0, length(g),length(g))*B'*P(:);
    else
        ccaDetails.c = [];
    end
    ccaDetails.P = P;
    ccaDetails.opts = OPTS;


%%%%%%%%%%%%%%%%%%%% FOLLOWING IS SUPPORTING MATLAB FUNCTIONS
function [A, b, c] = formulateSDP(S, D, bb, TRACE)
    [F0, FI, c] = localformulateSDP(S, D, bb, TRACE);
    [A, b, c] = sdpToSeDuMi(F0, FI, c);
    

function [F0, FI, c] = localformulateSDP(S, D, b, TRACE)
% formulate SDP problem
% each FI that corresponds to the LMI for the quadratic cost function has
% precisely 2*D^2 nonzero elements. But we need only D^2 storage for
% indexing these elements since the FI are symmetric
        tempFidx = zeros(D^2, 3);
        dimF = (D^2+1) + D + 2*TRACE;
        idx= 0;
        tracearray = ones(TRACE,1);
        for col=1:D
            for row=col:D
                idx = idx+1;
                lindx1 = sub2ind([D D], row, col);
                lindx2 = sub2ind([D D], col, row);
                tempFidx(:,1) = [1:D^2]';
                tempFidx(:,2) = D^2+1;
                if col==row
                    tempFidx(:,3) = S(:, lindx1) ;
                    FI{idx} = sparse([tempFidx(:,1); ...  % for cost function
                                        tempFidx(:,2); ... % symmetric
                                        row+D^2+1; ... % for P being p.s.d
                                        tracearray*(D^2+1+D+1); % for trace
                                        tracearray*(D^2+1+D+2); % for negate trace
                                    ], ...
                                    [tempFidx(:,2); ...  % for cost function
                                        tempFidx(:,1); ... % symmetric
                                        row+D^2+1; ... % for P being p.s.d
                                        tracearray*(D^2+1+D+1); % for trace
                                        tracearray*(D^2+1+D+2); % for negate trace
                                    ],...
                                    [tempFidx(:,3); ... % for cost function
                                        tempFidx(:,3); ... % symmetric
                                        1;                  % for P being p.s.d
                                        tracearray*1; % for trace
                                        tracearray*(-1); % for negate trace                               
                                    ], dimF, dimF);
                else

                    tempFidx(:,3) = S(:, lindx1) + S(:, lindx2);
                    FI{idx} = sparse([tempFidx(:,1); ...  % for cost function
                                        tempFidx(:,2); ... % symmetric
                                        row+D^2+1; ... % for P being p.s.d
                                        col+D^2+1; ... % symmetric
                                    ], ...
                                    [tempFidx(:,2); ...  % for cost function
                                        tempFidx(:,1); ... % symmetric
                                        col+D^2+1; ... % for P being p.s.d
                                        row+D^2+1; ... % being symmetric
                                    ],...
                                    [tempFidx(:,3); ... % for cost function
                                        tempFidx(:,3); ... % symmetric
                                        1;                  % for P being p.s.d
                                        1;                  % symmetric
                                    ], dimF, dimF);

                end
            end
        end
        idx=idx+1;
        % for the F matrix corresponding to t
        FI{idx} = sparse(D^2+1, D^2+1, 1, dimF, dimF);

        % now for F0
        if TRACE==1
            F0 = sparse( [[1:D^2] dimF-1 dimF], [[1:D^2] dimF-1 dimF], [ones(1, D^2) -1 1], dimF, dimF);
        else
            F0 = sparse( [[1:D^2]], [[1:D^2]], [ones(1, D^2)], dimF, dimF);
        end

        % now for c
        b = reshape(-b, D, D);
        b = b*2 - diag(diag(b)); 
        c = zeros(idx-1,1);
        kdx=0;
        %keyboard;
        for col=1:D
            for row=col:D
              kdx = kdx+1;
              c(kdx) = b(row, col);
            end
        end
        %keyboard;
        c = [c; 1]; % remember: we use only half of P
    return;


function [A, b, c] = sdpToSeDuMi(F0, FI, cc)
% convert the canonical SDP dual formulation:
% (see  Vandenberche and Boyd 1996, SIAM Review)
%  max -Tr(F0 Z)
% s.t. Tr(Fi Z) = cci and Z is positive definite
%
% in which cc = (cc1, cc2, cc3,..) and FI = {F1, F2, F3,...}
% 
% to SeDuMi format (formulated as vector decision variables ):
% min c'x
% s.t. Ax = b and x is positive definite (x is a vector, so SeDuMi
% really means that vec2mat(x) is positive definite)
%
% by feisha@cis.upenn.edu, June, 10, 2004

        if nargin < 3
            error('Cannot convert SDP formulation to SeDuMi formulation in sdpToSeDumi!');
        end

        [m, n] = size(F0);
        if m ~= n
            error('F0 matrix must be squared matrix in sdpToSeDumi(F0, FI, b)');
        end

        p = length(cc);
        if p ~= length(FI)
            error('FI matrix cellarray must have the same length as b in sdpToSeDumi(F0,FI,b)');
        end

        % should check every element in the cell array FI...later..

        % x = reshape(Z, n*n, 1);  % optimization variables from matrix to vector

        % converting objective function of the canonical SDP
        c = reshape(F0', n*n,1);

        % converting equality constraints of the canonical SDP
        zz= 0;
        for idx=1:length(FI)
            zz= zz + nnz(FI{idx});
        end
        A = spalloc( n*n, p, zz);
        for idx = 1:p
            temp = reshape(FI{idx}, n*n,1);
            lst = find(temp~=0);
            A(lst, idx) = temp(lst);
        end
        % The SeDuMi solver actually expects the transpose of A as in following
        % dual problem
        % max b'y
        % s.t. c - A'y is positive definite
        % Therefore, we transpose A
        % A = A';

        % b doesn't need to be changed
        b = cc;
    return;


    % Check OPTS that is passed into
    function OPTS = check_opt(OPTS)
        if isfield(OPTS,'method') == 0  
            OPTS.method = 'cca';
            disp('Options does''t have method field, so running CCA');
        end

        if strncmpi(OPTS.method, 'MVU',3)==1
            OPTS.CCA = 0; OPTS.MVU = 1;
        else
            OPTS.CCA = 1; OPTS.MVU = 0;
        end

        if isfield(OPTS, 'relative')==0
            OPTS.relative = 0;
        end

        if OPTS.CCA==1 && OPTS.relative ==1
            disp('Running CCA, so the .relative flag set to 0');
            OPTS.relative = 0;
        end

        if isfield(OPTS, 'regularizer')==0
            OPTS.regularizer = 0;
        end
    return



    function [tnn vidx]= triangNN(snn, TRI)
    % function [TNN VIDX]= triangNN(SNN) triangulates a sparse graph coded by spare matrix
    % SNN. TNN records the original edges in SNN as well as those that are
    % triangulated. Each edge is associated with a scaling factor that is specific
    % to a vertex. And VIDX records the id of the vertex.
    %
    % by feisha@cs.berkeley.edu  Aug. 15, 2006.

        N = size(snn,1);
        %fprintf('The graph has %d vertices\n', N);
        % figure out maximum degree a vertex has
        connectivs = sum(snn,1);
        maxDegree =  max(connectivs);
        tnn = spalloc(N, N, round(maxDegree*N)); % prealloc estimated storage for speedup 

        % triangulation
        for idx=1:N
            lst = find(snn(:, idx)>0);
            for jdx=1:length(lst)
                col = min (idx, lst(jdx));
                row = max(idx, lst(jdx));
                tnn(row, col) = tnn(row, col)+1;
                if TRI == 1
                    for kdx = jdx+1:length(lst)
                        col = min(lst(jdx), lst(kdx));
                        row = max(lst(jdx), lst(kdx));
                        tnn(row, col) = tnn(row, col)+1;
                    end
                end
            end
        end

        numVertexIdx = full(sum(tnn(:)));
        %fprintf('%d vertex entries are needed\n', numVertexIdx);

        rowIdx = zeros(numVertexIdx,1);
        colIdx = zeros(numVertexIdx,1);
        vidx = zeros(numVertexIdx,1);

        whichEdge = 0;
        for idx=1:N
            lst = find(snn(:, idx)>0);
            for jdx=1:length(lst)
                col = min(lst(jdx), idx);
                row  = max(lst(jdx), idx);
                whichEdge = whichEdge+1;
                rowIdx(whichEdge) = row;
                colIdx(whichEdge) = col;
                vidx(whichEdge)  = idx;
                if TRI==1
                    for kdx = jdx+1:length(lst)
                        col = min(lst(jdx), lst(kdx));
                        row = max(lst(jdx), lst(kdx));
                        whichEdge = whichEdge+1;
                        rowIdx(whichEdge) = row;
                        colIdx(whichEdge) = col;
                        vidx(whichEdge)  = idx;
                    end
                end
            end
        end
        linearIdx  = sub2ind([N N],rowIdx, colIdx);
        [sa, sIdx] = sort(linearIdx);
        vidx = vidx(sIdx);
    return


    % turn sparse graph snn into row and col indices
    function [edgesrow, edgescol, value] = sparse_nn(snn)
        N = size(snn,1);
        edgescol = zeros(N+1,1);
        nnzer = nnz(snn);
        edgesrow = zeros(nnzer,1);
        value = zeros(nnzer,1);

        edgescol(1) = 0;
        for jdx=1:N
            lst = find(snn(:, jdx)>0);
            %lst = lst(find(lst>jdx));
            edgescol(jdx+1) = edgescol(jdx)+length(lst);
            edgesrow(edgescol(jdx)+1:edgescol(jdx+1)) = lst-1;
            value(edgescol(jdx)+1:edgescol(jdx+1)) = snn(lst, jdx);
        end
    return