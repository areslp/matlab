function [C] = cotraining(views,num_views,numClust,projev,numiter,prior)
% co-trained spectral clustering algorithm (ICML, 2011)

% INPUTS:
% DATA: cell array of length 'num_views' of data matrices of size N x dim
% numClust: desired number of clusters
% SIGMA: array of width parameters for Gaussian kernel in each view
% TRUTH: ground truth clustering
% PROJEV: determines number of top eigenvectors of Laplacian onto which the projection
% is done in the algorithm (the paper considers projev=1 but it was later observed
% that considering more eigenvectors in the projection step helps),
% although the final k-mean clustering is run only on top-'numClust'
% eigenvectors of the graph Laplacian.
% NUMITER: number of iterations 

% OUTPUTS: 
% nmi_max: maximum nmi value obtained

kmeans_avg_iter = 10;
opts.disp = 0;

N=length(views{1});

numEV = numClust;
numVects = numClust;
for i=1:num_views    
    fprintf('computing kernel for view %d\n',i);
    K(:,:,i) = views{i};
    %K1 = X1*X1';
    D = diag(sum(K(:,:,i),1));
    %L1 = D1 - K1; 
    L(:,:,i) = sqrt(inv(D))*K(:,:,i)*sqrt(inv(D));  
    L(:,:,i)=(L(:,:,i)+L(:,:,i)')/2;
end

V = cell(num_views,1);
for i=1:num_views
    [V{i}] = mybaseline_spectral_onkernel(K(:,:,i),numClust,projev);
end

X = V;
Y = K; Y_norm = Y;
for i=1:numiter
    fprintf ('iteration %d...\n', i);
    Sall = zeros(N,N);
    for j=1:num_views
        Sall = Sall + X{j}*X{j}';
    end
    for j=1:num_views
        Y(:,:,j) = K(:,:,j)*(Sall - X{j}*X{j}');
        Y(:,:,j) = (Y(:,:,j) + Y(:,:,j)')/2; % + 1*eye(N);
        %Y_norm(:,:,j) = renormalize(Y(:,:,j));
        Y_norm(:,:,j) = Y(:,:,j);
        opts.disp = 0;
        [X{j}] = mybaseline_spectral_onkernel(Y_norm(:,:,j),numClust,projev);    
    end                
end 

if prior==0
    % no prior view
    AU=[]; % all U from views, column-wise concate
    for i=1:num_views
        AU=[AU,X{i}(:,1:ceil(numClust*1))];
    end
else
    AU=X{prior}(:,1:ceil(numClust*1));
end

size(AU)

if (1)
    norm_mat = repmat(sqrt(sum(AU.*AU,2)),1,size(AU,2));
    %%avoid divide by zero
    for i=1:size(norm_mat,1)
        if (norm_mat(i,1)==0)
            norm_mat(i,:) = 1;
        end
    end
    AU = AU./norm_mat;
end

C = kmeans(AU,numClust,'EmptyAction','drop');
    
% function  [K] = renormalize(K)
    
% mn = min(min(K));
% mx = max(max(K));
% if (mn < 0)
    % K = (K - mn) / (mx-mn);
    % K = (K+K')/2;
    % K = K - mn;
% end
    
