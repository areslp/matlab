clear;
load synth3views_2clusters;
clear sigma;
num_views = 3;
numClust = 2;
X{1} = X1;
X{2} = X2;
X{3} = X3;
sigma1 = optSigma(X1);
sigma2 = optSigma(X2);
sigma3 = optSigma(X3);
sigma(1) = sigma1;
sigma(2) = sigma2;
sigma(3) = sigma3;
lambda = 0.5; 
numiter = 5;

%% single best view
fprintf('Running with view1\n');
[E F P R nmi avgent] = baseline_spectral(X1,numClust,sigma1,truth);
fprintf('Running with view2\n');
[E F P R nmi avgent] = baseline_spectral(X2,numClust,sigma2,truth);
fprintf('Running with view3\n');
[E F P R nmi avgent] = baseline_spectral(X3,numClust,sigma3,truth);


%% two views
% feature concat
fprintf('Running with the feature concatenation of two views\n');
[E F P R nmi avgent] = baseline_spectral([X1 X2],numClust,optSigma([X1 X2]),truth);
%de sa mindis
fprintf('Running with the De Sa MinDis\n');
[F P R nmi avgent AR] = spectral_mindis(X1, X2, numClust,sigma1,sigma2, truth);
fprintf('Our approach with 2 views (pairwise)\n');
[U U1 F P R nmi avgent AR] = spectral_pairwise(X1,X2,numClust,sigma1,sigma2,lambda,truth,numiter);
fprintf('Our approach with 2 views (centroid)\n');
lambda1 = 0.5; lambda2 = 0.5;
[U U1 F P R nmi avgent AR] = spectral_centroid(X1,X2,numClust,sigma1,sigma2,lambda1,lambda2,truth,numiter);


%% three views
fprintf('Running with the feature concatenation of three views\n');
[E F P R nmi avgent] = baseline_spectral([X1 X2 X3],numClust,optSigma([X1 X2 X3]),truth);

% multiview spectral (pairwise): more than 2 views
fprintf('Multiview spectral with 3 views\n');
[F P R nmi avgent AR] = spectral_pairwise_multview(X,num_views,numClust,sigma,lambda,truth,numiter);

% multiview spectral (centroid): more than 2 views
fprintf('Multiview spectral with 3 views\n');
lambda = [0.5 0.5 0.5];
[F P R nmi avgent AR] = spectral_centroid_multiview(X,num_views,numClust,sigma,lambda,truth,numiter);
