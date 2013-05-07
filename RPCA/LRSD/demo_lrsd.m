function demo_lrsd
% Demo to solve
% min_A tau*||A(:)||_1 + ||C - A||_*
% min_A t*||A(:)||_1 + (1 - t) * ||C - A||_*
%
% tau = t   / (1 - t);
% t   = tau / (1 + tau);
%
close all
%% generate problem
% problem size
m = 100;
n = m;

%% Change settings here to est different problems
%==================================================================
impulsive = 1; % 1 for Gaussian and 0 for impulsive Sparse matrices
spr = 0.05; % sparsity ratio: #nonzeros/m/n
rB = 10; % rank of Low-Rank matrix
%==================================================================

% Low-Rank matrix
B = randn(m,rB) * randn(rB,n);

% Sparse matrix
A = zeros(m,n);
p = randperm(m*n);
L = round(spr*m*n);
A(p(1:L)) = randn(L,1);
if impulsive; mgB = max(abs(B(:))); A(p(1:L)) = mgB * sign(randn(L,1)); end


% Sparse + Low-Rank
C = A + B;
% Prob1;

% regularization parameter
t = .1;

%% algorithm
opts = [];
opts.beta =   .25/mean(abs(C(:)));%0.10;
opts.tol = 1e-6;
opts.maxit = 1000;
opts.A0 = zeros(m,n);
opts.B0 = zeros(m,n);
opts.Lam0 = zeros(m,n);
opts.Sparse = A;
opts.LowRank = B;
opts.print = 1;
out = lrsd(C, t/(1-t), opts);

errSP = norm(out.Sparse - A, 'fro') / (1 + norm(A,'fro'));
errLR = norm(out.LowRank - B, 'fro') / (1 + norm(B,'fro'));
errTotal = norm([out.Sparse,out.LowRank] - [A,B],'fro') / (1 + norm([A,B],'fro'));

fprintf('Image size %d x %d, Rank %d, Sparse Ratio %4.1f%%\n',m,n,rB,spr*100);
fprintf('Iter: %d, Sparse Error: %4.2e, Low-Rank Error: %4.2e, Total Error: %4.2e\n',out.iter,errSP,errLR,errTotal);

tx = 13;  
figure(1);  set(1,'position',[0,100,400,400]);
subplot(121); imshow(A,[]);     title('True Sparse','fontsize',tx);
subplot(122); imshow(out.Sparse,[]);  title('Recovered Sparse','fontsize',tx);
colormap(bone(5));

figure(2); set(2,'position',[420,100,400,400]);
subplot(121); imshow(B,[]); title('True Low-Rank','fontsize',tx);
subplot(122); imshow(out.LowRank,[]); title('Recovered Low-Rank','fontsize',tx);
colormap(bone(5));

figure(3); set(3,'position',[840,100,400,400]);
subplot(121); semilogy(out.errsSP,'-ro','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',3); ylabel('RelErr in Sparse Matrix');
subplot(122); semilogy(out.errsLR,'-bo','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',3); ylabel('RelErr in Low-Rank Matrix');
end
