clear;
load synth3views_2clusters;
clear sigma;
num_views = 3;
numClust = 2;
data{1} = X1;
data{2} = X2;
data{3} = X3;
sigma1 = optSigma(X1);
sigma2 = optSigma(X2);
sigma3 = optSigma(X3);
sigma(1) = sigma1;
sigma(2) = sigma2;
sigma(3) = sigma3;
projev = 1.5; 

%% single best view
for j=1:num_views
    fprintf('Running with view %d\n', j);
    options.KernelType = 'Gaussian';
    options.t = sigma(j);
    K = constructKernel(data{j},data{j},options);
    [V Eval F P R nmi avgent AR C] = baseline_spectral_onkernel(K,numClust,truth,projev);
end

%% two views
% feature concat
fprintf('Running with the feature concatenation of two views\n');
[E F P R nmi avgent] = baseline_spectral([X1 X2],numClust,optSigma([X1 X2]),truth);


%% three views
fprintf('Running with the feature concatenation of three views\n');
[E F P R nmi avgent] = baseline_spectral([X1 X2 X3],numClust,optSigma([X1 X2 X3]),truth);

%% co-trained multiview spectral
projev = 1.5; 
numiter = 5;
fprintf('Co-trained multiview spectral with 3 views\n');
[nmi_max] = spectral_cotraining(data,num_views,numClust,sigma,truth,projev,numiter);

