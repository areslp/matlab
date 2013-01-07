function [V] = baseline_spectral_onkernel(K,numClust,projev)
% spectral clustering using similarity matrix K. Clustering is done by
% running k-means on the top-'numClust' eigenvectors of the normalized Laplacian
% INPUT:
% K: N x N similarity matrix. 
% numClust: desired number of clusters
% truth: N x 1 vector of ground truth clustering
% projev: number of top eigenvectors to return in V
% OUTPUT:
% V: top-'projev' eigenvectors of the Laplacian
% E: top-'projev' eigenvalues of the Laplacian
% F, P, R, nmi, avgent, AR: F-score, Precision, Recall, normalized mutual
% information, average entropy, Avg Rand Index
% C: obtained cluster labels 

    numEV = numClust*projev;
    N = size(K,1);
    %options.KernelType = 'Gaussian';
    %options.t = sigma; % width parameter for Gaussian kernel
    %fprintf('constructing kernel...\n');
    %K = constructKernel(X,X,options);
    
    D = diag(sum(K,1));
    inv_sqrt_D = sqrt(inv(abs(D)));
    L = inv_sqrt_D*K*inv_sqrt_D;     
    %L = inv(D) * K;    
    L = (L+L')/2;    
    %sum(sum(L-L'));
    % now do an eigen-decomposition of L
    %fprintf('doing eigenvalue decomp...\n');
    opts.disp = 0;
    [V E] = eigs(L,ceil(numEV),'LA',opts);  
    U = V(:,1:ceil(numClust*1));
    % fprintf(1,'single view E:\n');
    % E
    
    %[U E] = eig(L);   
    %[E1 I] = sort(diag(E));  %sort in increasing order
    %U = U(:,I(end-numEV+1:end));
    if (1)
    norm_mat = repmat(sqrt(sum(U.*U,2)),1,size(U,2));
    %%avoid divide by zero
    for i=1:size(norm_mat,1)
        if (norm_mat(i,1)==0)
            norm_mat(i,:) = 1;
        end
    end
    U = U./norm_mat;
    end
    %fprintf('running k-means...\n');
    
    % C = kmeans(U,numClust,'EmptyAction','drop');
