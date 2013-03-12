function  [P, newY, L, newV, idx]= sdecca2(Y, snn, regularizer, relative)
% doing semidefinitve embedding/MVU with output being parameterized by graph
% laplacian's eigenfunctions.. 
%
% the algorithm is same as conformal component analysis except that the scaling
% factor there is set as 1
%
%
% function [P, newY, Y] = CDR2(X, Y, NEIGHBORS)  implements the 
% CONFORMAL DIMENSIONALITY REDUCTION of data X. It finds a linear map
% of Y -> L*Y such that X and L*Y is related by a conformal mapping.
%
% No tehtat  The algorithm use the formulation of only distances.
%
% Input:
%   Y: matrix of d'xN, with each column is a point in R^d'
%   NEIGHBORS: matrix of KxN, each column is a list of indices (between 1
%   and N) to the nearest-neighbor of the corresponding column in X
% Output:
%   P: square of the linear map L, i.e., P = L'*L
%   newY: transformed data points, i.e., newY = L*Y;
%   Y: the linear map L itself,    i.e., L = L
%
% The algorithm finds L by solving a semidefinite programming problem. It
% calls csdp() SDP solver by default and assumes that it is on the path.
%
% written by feisha@cis.upenn.edu
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    % Collect the data
    [erow, ecol, edist] = sparse_nn(snn);
    irow = int32(erow); 
    icol = int32(ecol);
    [A, B, g] = mexCCACollectData2(Y, irow, icol, edist, int32(relative));
    BG = 2 * sum(B, 2);
    Q = A ;
    [V, E] = eig(Q + eye(size(Q)));
    E = E - eye(size(Q));
    E(E < 0) = 0;
    if ~isreal(diag(E))
        E = real(E);
        V = real(V);
        S = sqrt(E) * V';
    else
        S = sqrt(E) * V';
    end

    % Put the regularizer in there
    BG = BG + regularizer * reshape(eye(size(Y, 1)), size(Y, 1) ^ 2, 1);
    
    % Formulate the SDP problem
    [AA, bb, cc] = formulateSDP(S, size(Y, 1), BG);
    sizeSDP = size(Y, 1) ^ 2 + 1 + size(Y, 1);
    pars.s = sizeSDP;
    opts.printlevel = 1;

    % Solve it using CSDP
    [xx, yy] = csdp(AA, bb, cc, pars, opts);

    % The negate of yy is our solution
    yy = -yy;
    idx = 0;
    P = zeros(size(Y, 1));
    for col=1:size(Y, 1)
        for row = col:size(Y, 1)
            idx = idx + 1;
            P(row, col) = yy(idx);
        end
    end
    
    % Convert P to a positive definite matrix
    P = P + P' - diag(diag(P));

    % Transform the original projection to the new projection
    [V, E] = eig(P);
    E(E < 0) = 0;
    L = diag(sqrt(diag(E))) * V';
    newY = L * Y;               % multiply with Laplacian

    % Eigendecomposition of the new projection: doing PCA because the 
    % dimensionality of newY or Y is definitely less than the number of
    % points
    [newV, newE] = eig(newY * newY');
    newE = diag(newE);
    [dummy, idx] = sort(newE);
    newY = newV' * newY;
    newY = newY(idx(end:-1:1),:);
return


% Function that formulates the SDP problem
function [A, b, c]=formulateSDP(S, D, bb)
    [F0, FI, c] = localformulateSDP(S, D, bb);
    [A, b, c] = sdpToSeDuMi(F0, FI, c);
return


% Function that formulates the SDP problem
function [F0, FI, c] = localformulateSDP(S, D, b)

    % Each FI that corresponds to the LMI for the quadratic cost function has
    % precisely 2 * D^2 nonzero elements. But we need only D^2 storage for
    tempFidx = zeros(D ^ 2, 3);
    dimF = (D ^ 2 + 1) + D;
    idx = 0;
    for col=1:D
        for row=col:D
            idx = idx + 1;
            lindx1 = sub2ind([D D], row, col);
            lindx2 = sub2ind([D D], col, row);
            tempFidx(:,1) = [1:D ^ 2]';
            tempFidx(:,2) = D ^ 2 + 1;
            if col == row
                tempFidx(:,3) = S(:,lindx1) ;
                FI{idx} = sparse([tempFidx(:,1); ...    % for cost function
                    tempFidx(:,2); ...                  % symmetric
                    row + D^2 + 1 ...                   % for P being p.s.d

                    ], ...
                    [tempFidx(:,2); ...                 % for cost function
                    tempFidx(:,1); ...                  % symmetric
                    row + D^2 + 1; ...                  % for P being p.s.d

                    ],...
                    [tempFidx(:,3); ...                 % for cost function
                    tempFidx(:,3); ...                  % symmetric
                    1;                                  % for P being p.s.d

                    ], dimF, dimF);
            else

                tempFidx(:,3) = S(:, lindx1) + S(:,lindx2);
                FI{idx} = sparse([tempFidx(:,1); ...    % for cost function
                    tempFidx(:,2); ...                  % symmetric
                    row + D^2 + 1; ...                  % for P being p.s.d
                    col + D^2 + 1; ...                  % symmetric
                    ], ...                    
                    [tempFidx(:,2); ...                 % for cost function
                    tempFidx(:,1); ...                  % symmetric
                    col + D^2 + 1; ...                  % for P being p.s.d
                    row + D^2 + 1; ...                  % being symmetric
                    ], ...                    
                    [tempFidx(:,3); ...                 % for cost function
                    tempFidx(:,3); ...                  % symmetric
                    1;                                  % for P being p.s.d
                    1;                                  % symmetric
                    ], dimF, dimF);

            end
        end
    end
    idx = idx + 1;
    
    % For the F matrix corresponding to t
    FI{idx} = sparse(D^2 + 1, D^2 + 1, 1, dimF, dimF);

    % Now for F0
    F0 = sparse(1:D^2, 1:D^2, ones(1, D^2), dimF, dimF);

    % Now for c
    b = reshape(-b, D, D);
    b = b * 2 - diag(diag(b));
    c = zeros(idx - 1,1);
    kdx = 0;
    for col=1:D
        for row=col:D
            kdx = kdx + 1;
            c(kdx) = b(row, col);
        end
    end
    c = [c; 1];
return


% Function that convertsthe canonical SDP dual formulation to SeDuMi format
function [A, b, c] = sdpToSeDuMi(F0, FI, cc)

    % Check inputs
    if nargin < 3
        error('Cannot convert SDP formulation to SeDuMi formulation.');
    end
    [m, n] = size(F0);
    if m ~= n
        error('F0 matrix must be squared matrix.');
    end
    p = length(cc);
    if p ~= length(FI)
        error('FI matrix cellarray must have the same length as b.');
    end

    % Converting objective function of the canonical SDP
    c = reshape(F0', n * n, 1);

    % Converting equality constraints of the canonical SDP
    zz = 0;
    for idx=1:length(FI)
        zz= zz + nnz(FI{idx});
    end
    A = spalloc(n * n, p, zz);
    for idx=1:p
        temp = reshape(FI{idx}, n * n, 1);
        lst = find(temp ~= 0);
        A(lst, idx) = temp(lst);
    end

    % We do not need to convert b
    b = cc;
return
