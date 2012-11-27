function [V E F P R nmi avgent AR] = baseline_spectral(X,numClust,sigma,truth)
% INPUT:
% X: N x P data matrix. Each row is an example
% numClust: desired number of clusters
% truth: N x 1 vector of ground truth clusterings
% OUTPUT:
% C, U, F, P, R: clustering, U matrix, F-score, Precision, Recall

    if (min(truth)==0)
        truth = truth + 1;
    end
    numEV = numClust*1;
    N = size(X,1);
    options.KernelType = 'Gaussian';
    options.t = sigma; % width parameter for Gaussian kernel
    fprintf('constructing kernel...\n');
    K = constructKernel(X,X,options);
    K = K;
    D = diag(sum(K,1));
    L = sqrt(inv(D))*K*sqrt(inv(D));    
    L = (L+L')/2;
    %convert_libsvm(L,'handwritten_L1.vb');
    % now do an eigen-decomposition of L
    fprintf('doing eigenvalue decomp...\n');
    opts.disp = 0;
    [V E] = eigs(L,numEV,'LA',opts);  
    U = V(:,1:numClust);
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
    fprintf('running k-means...\n');
    
    for i=1:20
        C = kmeans(U,numClust,'EmptyAction','drop');
        [A nmii(i) avgenti(i)] = compute_nmi(truth,C);
        [Fi(i),Pi(i),Ri(i)] = compute_f(truth,C);
        [ARi(i),RIi(i),MIi(i),HIi(i)]=RandIndex(truth,C);
    end
    F = mean(Fi);
    P = mean(Pi);
    R = mean(Ri);
    nmi = mean(nmii);
    avgent = mean(avgenti);
    AR = mean(ARi);
    
%    fprintf('F: %f(%f)\n', F, std(Fi));
%    fprintf('P: %f(%f)\n', P, std(Pi));    
%    fprintf('R: %f(%f)\n', R, std(Ri));
    fprintf('nmi: %f(%f)\n', nmi, std(nmii));
%    fprintf('avgent: %f(%f)\n', avgent, std(avgenti));
%    fprintf('AR: %f(%f)\n', AR, std(ARi));
