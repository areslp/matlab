function [V E F P R nmi avgent AR C] = baseline_spectral_onkernel(K,numClust,truth)
% INPUT:
% X: N x P data matrix. Each row is an example
% numClust: desired number of clusters
% truth: N x 1 vector of ground truth clusterings
% OUTPUT:
% C, U, F, P, R: clustering, U matrix, F-score, Precision, Recsigall

    numEV = numClust*1.5;
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
    
    for i=1:10
        C = kmeans(U,numClust,'EmptyAction','drop');
        [Fi(i),Pi(i),Ri(i)] = compute_f(truth,C);
        [A nmii(i) avgenti(i)] = compute_nmi(truth,C);
        if (min(truth)==0)
            [ARi(i),RIi(i),MIi(i),HIi(i)]=RandIndex(truth+1,C);
        else
            [ARi(i),RIi(i),MIi(i),HIi(i)]=RandIndex(truth,C);
        end       
    end
    F(1) = mean(Fi); F(2) = std(Fi);
    P(1) = mean(Pi); P(2) = std(Pi);
    R(1) = mean(Ri); R(2) = std(Ri);
    nmi(1) = mean(nmii); nmi(2) = std(nmii);
    avgent(1) = mean(avgenti); avgent(2) = std(avgenti);
    AR(1) = mean(ARi); AR(2) = std(ARi);
    
    %fprintf('F: %f(%f)\n', F(1), std(Fi));
    %fprintf('P: %f(%f)\n', P(1), std(Pi));    
    %fprintf('R: %f(%f)\n', R(1), std(Ri));
    fprintf('nmi: %f(%f)\n', nmi(1), std(nmii));
    %fprintf('avgent: %f(%f)\n', avgent(1), std(avgenti));
    %fprintf('AR: %f(%f)\n', AR(1), std(ARi));
    
    