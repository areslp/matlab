classdef pseudoinverse  % < factorization_generic
    % PIF = PSEUDOINVERSE(A)
    %
    % Implement the Moore-Penrose pseudo-inverse factorization object 'PIF'
    % on a matrix A so that later call
    %   x = PIF*b % returns the same result as
    %   x = pinv(A)*b % (beside floating point roundoff error)
    %
    % Similarly a product with a standard matrix on a left side
    %   y = d*PIF  % returns the same result as
    %   y = d*pinv(A)
    %
    % Also feature is PSEUDOINVERSE is output an iterative Tikhonov's
    % regularization object to solve similar kind of linear problem:
    %       A*x = b + noise (see 'LSQR' option)
    %
    % NOTES:
    %   - Product with the pseudo inverse is the convergent limit of
    %     Tikhonov's regularized solution when the regularization term goes
    %     to zero.
    %   - Work for full and sparse matrix (pinv does not supported sparse).
    %   - PSEUDOINVERSE(A, TOL) where TOl is a scalar in order to specify
    %     QR-tolerance for cutting range.
    %     Default tolerance used is max(size(A))*norm(A)*eps. The tolerance
    %     is not equivalent to the same parameter used by PINV.
    %   - Comment the attribute (Hidden = true) in line #131 for older
    %     versions of Matlab
    %   - For release 7.8 and anterior, sparse matrix is converted to full
    %     unless if SuiteSparseQR package from Tim Davis is installed
    %     (available from http://www.cise.ufl.edu/research/sparse/SPQR )
    %
    % PSEUDOINVERSE(A, TOL, 'LSQR') will return the least-square
    %     Tikhonov's regularization object PIF such that command
    %       x = PIB*b; returns x that minimizes the regularized least-square
    %           J(x) = |A*x - b|^2 + tol^2 |x|^2
    %     (The minimization algorithm is carried out using Matlab PCG)
    %
    % PSEUDOINVERSE(A, TOL, 'LSQR', 'Property1, Value1, ...) allows
    %    to customize the inversion. Valid properties are:
    %    'pcgtol', specified the iterative relative tolerance, default 1e-6
    %    'maxit', specified the maximum number of iterations, default
    %             min([m n 20])
    %    'M1', 'M2': preconditioner matrices, default [] (no
    %                preconditioning), see help PCG for more details
    %    'Tikhonov': specify other regularization. Possible forms are:
    %       - scalar mu: negative -> no reguarization
    %                    positive J(x) := |A*x - b|^2 + (mu*|x|)^2
    %       - matrix L:  J(x) := |A*x - b|^2 + |L*x|^2
    %       - function handle that computes R*x
    %         in the regularized objective: J(x) := |A*x - b|^2 + x'*R*x
    %         R is a user-supplied symmetric (semi-)norm function.
    %         Optional arguments can be passed as cell {Rfun, arg2, arg3, ...}
    %         Rfun will be called Rfun(x, arg2, arg3, ...) and should
    %         return the vector R*x
    %   'verbose': true/false -> turn on/off PCG message display. 'false'
    %              by default
    %
    %   All or some of the above properties can be grouped together in a
    %   single structure, named 'LSQROpt' with the Properties specified in
    %   the fields.
    %
    % NOTE: "Dual" Tikhonov's regularization form
    %       y = argmin J(y)
    %           J(y) := |A'*y - d|^2 + tol^2 |y|^2
    %  can be carried out by calling left product
    %       y = d'*pseudoinverse(A,[],'lsqr')
    %       y = y'
    %  A non-standard Tikhonov's regularization (using e.g. matrix) must
    %  specifically be configured for exlusively either primal or dual form 
    %
    % EXAMPLE:
    %
    %   % Matrix
    %   m1 = 7; m2=5;
    %   n = 5;
    %   rankA = 3;
    %   d = zeros(n,1);
    %   d(1:rankA) = rand(1,rankA);
    %   L1 = rand(n,m1);
    %   L2 = rand(n,m2);
    %   A = sparse(L1'*diag(d)*L2);
    %
    %   F=pseudoinverse(A);
    %
    %   % RHS
    %   b = rand(size(F,2),3);
    %
    %   % Checking
    %   x1 = F*b
    %   x2 = pinv(full(A))*b
    %   norm(x1-x2,'fro')/norm(x2,'fro')
    %
    %   % Regularized solution with coefficient relative to matrix norm 
    %   Fr = pseudoinverse(A,[],'lsqr',...
    %                     'tikhonov', {@(x,r) r*normest(A)*x, 1e-8});
    %   xregularized = Fr*b
    %
    % See also: PINV, PCG, FACTORIZE package by Tim Davis (on File
    %           Exchange)
    %
    % Author: Bruno Luong <brunoluong@yahoo.com>
    % History
    %   Original: 30-Sep-2009
    %   03-Oct-2009: LSQR with Tikhonov regularization
    %   08-Oct-2009: Put static methods to separate private functions
    %                for better compatibility
    %   18-Oct-2009: Newly methods supported: left/right multiplication,
    %                conjugate, transpose, and complex-transpose
    %   19-Oct-2009: call SuiteSparseQR when detected
    
    %%
    properties %(Hidden=true)
        % list of local variables that have been assigned through the
        % workspace
        type = 'Pseudo-inverse factorization';
        method = 'QR';
        is_cmplxtrans = false;
        A;
        Qt;
        Tt;
        ES;
        Msize;
        Mrank;
        tol;
        Misparse;
        LSQROpt;

    end % properties (Hidden=true)
    
    %%
    methods (Hidden=true) % Comment the attribute (Hidden=true) if
                          % your Matlab does not support it
        
        %%
        function PIF = pseudoinverseQR(PIF, A, tol, options) %#ok
        % Pseudo-inverse by QR factorization in both source and destination
        % vectorial spaces
        
            persistent SPARSEQR_WARN

            % A(:,p) = A*E = Q*R
            % E'*A*A*E = R'*R
            % Q'*Q = I
            % R triangular         
            [m n] = size(A);
            
            % Permutation QR for sparse exists only from Release 7.9 ;-(
            if issparse(A) && ...
               datenum(version('-date')) < datenum('August 12, 2009')
                if isempty(which('spqr'))
                    % issue a warning
                    if isempty(SPARSEQR_WARN)
                        warning('PSEUDOINVERSE:SPARSEQRWARN', ...
                                'PSEUDOINVERSE: convert to full matrix');
                        fprintf('Recommendation: install SuiteSparseQR');
                        SPARSEQR_WARN = 1;
                    end
                    % no choice but convert A to full matrix
                    [Q R p] = qr(full(A),0);
                else
                    % SuiteSparseQR package
                    [Q R p] = spqr(A,0);
                end
            else
                % Matlab built-in QR
                [Q R p] = qr(A,0);
            end
            
            % reverse the permutation
            ip = zeros(size(p),'double');
            ip(p) = 1:length(p);
            
            % Default cutting rank
            if nargin<2 || isempty(tol)
                wrn = warning('off','MATLAB:normest:notconverge');
                tol = max(m,n) * normest(A) * eps(class(A));
                warning(wrn.state,'MATLAB:normest:notconverge');
            end
            % estimate the rank
            % Need to convert to full because of a Bug in 2007B
            % that can't deal with find on SPARSE
            k = find(abs(full(diag(R)))<=tol,1,'first');
            
            if k==0 % not go inside if block for empty k
                % Matrix is all zeros
                % Pseudo inverse is zeros
                PIF.Qt = sparse([],[],[],n,m);
                PIF.Tt = speye(n);
                PIF.ES = speye(n);
            else
                if isempty(k)
                    k = min(m,n);
                else
                    k = k-1;
                end
                PIF.Qt = Q(:,1:k)';
                R = R(1:k,:);
                
                if k<n
                    % Factorization R';
                    % R' = S*T
                    % R = Tt * S'
                    % S'*S = I
                    % Tt triangular
                    [S T] = qr(R');
                    PIF.Tt = T';
                    
                    % S = E*S
                    PIF.ES = S(ip,:);
                else % full rank
                    PIF.Tt = R;
                    PIF.ES = sparse(1:length(ip),ip,1);
                end
            end
            PIF.Msize = [n m];
            PIF.Mrank = k;
            PIF.tol = tol;
            PIF.Misparse = issparse(A);
            PIF.method = 'QR';
        end % pseudoinverseQR
        
        function x = mtimesQR(PIF, b, dual)
            if ~dual
                % Solve A*x = b from QR factorization
                x = (PIF.ES * (PIF.Tt \ (PIF.Qt*b)));
            else
                % Solve x'*A = b from QR factorization
                x = (((b*PIF.ES) / (PIF.Tt)) * PIF.Qt);
            end % mtimesQR
        end
        
        %%
        function PIF = pseudoinverseLSQR(PIF, A, tol, options)
        % Least-square constructor
            if nargin<2 || isempty(tol)
                tol = 1e-6;
            end
            PIF.A = A;
            [m n] = size(A);
            PIF.Msize = [n m];
            PIF.Mrank = NaN;
            PIF.tol = tol;
            PIF.Misparse = issparse(A);         
            PIF.method = 'LSQR';
            
            % Default option for least square pcg
            DefLSQROpt = struct('tikhonov', tol, ...
                                'pcgtol', 1e-6, ... 
                                'maxit', min([m n 20]), ...
                                'm1', [], ...
                                'm2', [], ...
                                'verbose', false);
            % If LSQROpt is provided                
            DefLSQROpt = getoption(options, 'LSQROpt', DefLSQROpt);
            % parsing other options
            for field = fieldnames(DefLSQROpt)'
                fname = field{1};
                DefLSQROpt.(fname) = getoption(options,...
                                         fname, DefLSQROpt.(fname));
            end
            PIF.LSQROpt = DefLSQROpt;
        end % pseudoinverseLSQR
        
        function x = mtimesLSQR(PIF, b, dual)
        % Solve by least-square by conjugate gradient method
            
            Opt = PIF.LSQROpt;
            Amat = PIF.A;
            Lambda = Opt.tikhonov;
            
            if dual
                Amat = Amat';
                b = b';
            end
            
            x = zeros(size(Amat,2),size(b,2));
            
            % function that computes the Tikhonov's regularization
            % gradient
            function reg = regfun(x)  
                if isnumeric(Lambda)
                    if isscalar(Lambda)
                        % "Standard" form of Tikhonov: R = 0.5*Lambda^2*|x|^2
                        if Lambda > 0
                            reg = Lambda^2 * x;
                        else
                            reg = 0;
                        end
                    else isnumeric(Lambda) % seminorm Tikhonov matrix
                        % Compute R = 0.5*|L*x|^2
                        % Compute dR/dx = L'*L*x
                        Lx = Lambda*x;
                        reg = (Lx'*Lambda)';
                    end
                else % function handle that must compute H*x where
                     % H = L'*L
                    if isa(Lambda,'function_handle')
                        Lambda = {Lambda};
                    elseif ~iscell(Lambda)
                        error('Pseudo inverse requires matrix for first input');
                    end
                    reg = feval(Lambda{1},x,Lambda{2:end});
                end
            end % regfun

            function y = Anonsym(x)
            % Perform product y = (A'*A)*x
                Ax = Amat*x;
                y = (Ax'*Amat)';
                y = y + regfun(x);
            end % Anonsym     
               
            rhs = (b'*Amat)'; % A'*b
            Afun = @Anonsym;
            pcgtol = getoption(Opt, 'pcgtol', 1e-6);
            
            pcgmaxit = getoption(Opt, 'maxit', min([size(PIF) 20]));
            pcgM1 = getoption(Opt, 'm1', []);
            pcgM2 = getoption(Opt, 'm2', []);
            % Call iterative solver PCG for each rhs
            if ~getoption(Opt, 'verbose', false)
                pcgout = {[] []}; % no verbose
            else
                pcgout = {[]}; % with verbose
            end
            % Loop over vectors of RHS
            for k=1:size(b,2)
                [pcgout{:}] = pcg(Afun, rhs(:,k), pcgtol, pcgmaxit, ...
                                  pcgM1, pcgM2, x(:,k)); % pcgflag is needed
                                      % in order to disable verbose message
                x(:,k) = pcgout{1};            
            end
            
            if dual
                x = x';
            end
            
        end % mtimesLSQR
        
    end % hidden methods
    
    %%
    methods % public
       
        function PIF = pseudoinverse(A, tol, method, varargin)
            % Object constructor
            % PIF = pseudoinverse(A)
            % PIF = pseudoinverse(A, tol)
            % Factorization A for the purpose of performing pseudo inverse
            % 
            % PIF = pseudoinverse(A, tol, method, ...)
            % method: 'QR' 'LSQR'
            % Upcomming methods: 'SVD'
            
            if ndims(A)>2
                error('Pseudo inverse requires matrix for first input');
            end            
            
            if nargin<3 || isempty(method)
                method = 'QR';
            end
            if nargin<2 || isempty(tol)
                tol = [];
            end
            
            if nargin<4 || isempty(varargin{1})
                options = struct();
            else
                options = varargin{1};
                % Put Property1, val1, etc... into a structure
                if ~isstruct(options)
                    if mod(length(varargin),2)
                        error('Optional inputs must be Property/value pairs');
                    end
                    clear options
                    for k=1:2:length(varargin)
                        field = varargin{k};
                        value = varargin{k+1};
                        options.(field) = value;
                    end
                end
            end
            options = lowercase(options);

            switch lower(method)
                case 'qr',
                    PIF = PIF.pseudoinverseQR(A, tol, options);
                case 'lsqr',
                    PIF = PIF.pseudoinverseLSQR(A, tol, options);
                otherwise
                    error('Not implemented method %s', method);
            end
        end % constructor
        
        function s = size(PIF, n)
            % s = size(PIF)
            % Get the size of the pseudoinverse
            s = PIF.Msize;
            if PIF.is_cmplxtrans
                s = s([2 1]);
            end
            if nargin>=2
                if n<=length(s)
                    s = s(n);
                else % trailing dimension
                    s = 1;
                end
            end
        end % size
        
        function n = numel(PIF)
            n = prod(size(PIF)); %#ok
        end % numel
        
        function n = length(PIF)
            n = max(size(PIF)); %#ok
        end % length
        
        function x = mtimes(l, r)
        % Implement left/right product with regular matrix
            
            if isa(l,'pseudoinverse')
                % PIF*b
                [PIF b] = deal(l, r);
            else
                % b*PIF
                [PIF b] = deal(r, l);
            end
            
            if PIF.is_cmplxtrans
                b = b';
            end
            
            dual = ~xor(isa(l,'pseudoinverse'),PIF.is_cmplxtrans);            
            switch lower(PIF.method)
                case 'qr',
                    x = PIF.mtimesQR(b, dual);
                case 'lsqr',
                    x = PIF.mtimesLSQR(b, dual);
                otherwise
                    error('Not implemented method %s', PIF.method);
            end
            
            if PIF.is_cmplxtrans
                x = x';
            end

        end % mtimes
        
        function PIF = conj(PIF)
            % Implement Q = conj(P)
            PIF.A = conj(PIF.A);
            PIF.Qt = conj(PIF.Qt);
            PIF.Tt = conj(PIF.Tt);
            PIF.ES = conj(PIF.ES);
        end
        
        function PIF = ctranspose(PIF)
            % implement Q = P'
            PIF.is_cmplxtrans = ~PIF.is_cmplxtrans;
        end
        
        function PIF = transpose(PIF)
            % implement Q = P.'
            PIF = conj(ctranspose(PIF));
        end
        
    end % public methods
    
end % pseudoinverse
