function [nmi_max] = spectral_cotraining(data,num_views,numClust,sigma,truth,projev,numiter)
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


    if (min(truth)==0)
        truth = truth + 1;
    end        
    [N M1] = size(data{1});
    %[N M2] = size(data);
    
    for i=1:num_views
        %options(i) = [];
        options(i).KernelType = 'Gaussian';
        options(i).t = sigma(i);
        options(i).d = 4;
    end
        
    kmeans_avg_iter = 10;
    opts.disp = 0;
    
    numEV = numClust;
    numVects = numClust;
    for i=1:num_views    
        fprintf('computing kernel for view %d\n',i);
        K(:,:,i) = constructKernel(data{i},data{i},options(i));
        %K1 = X1*X1';
        D = diag(sum(K(:,:,i),1));
        %L1 = D1 - K1; 
        L(:,:,i) = sqrt(inv(D))*K(:,:,i)*sqrt(inv(D));  
        L(:,:,i)=(L(:,:,i)+L(:,:,i)')/2;
    end
 
    V = cell(num_views,1);
    for i=1:num_views
        [V{i} Eval F P R nmii avgent AR C] = baseline_spectral_onkernel(K(:,:,i),numClust,truth,projev);
        nmi(i,1) = nmii(1);    nmi(i,2) = nmii(2);
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
            [X{j} Eval F P R nmii avgent AR C] = baseline_spectral_onkernel(Y_norm(:,:,j),numClust,truth,projev);    
            nmi_K(i,1,j) = nmii(1); F_K(i,1,j) = F(1); P_K(i,1,j) = P(1); R_K(i,1,j) = R(1); 
            avgent_K(i,1,j) = avgent(1); AR_K(i,1,j) = AR(1);  
            nmi_K(i,2,j) = nmii(2); F_K(i,2,j) = F(2); P_K(i,2,j) = P(2); R_K(i,2,j) = R(2); 
            avgent_K(i,2,j) = avgent(2); AR_K(i,2,j) = AR(2);                      
        end                
    end 
    
    for j=1:num_views
        fprintf ('\nview %d: \n', j);
        fprintf ('\nF(%d): ', j);
        for i=1:numiter
            fprintf ('%f(%f) ', F_K(i,1,j), F_K(i,2,j));
        end
        fprintf ('\nP(%d): ', j);
        for i=1:numiter
            fprintf ('%f(%f) ', P_K(i,1,j), P_K(i,2,j));
        end
        fprintf ('\nR(%d): ', j);
        for i=1:numiter
            fprintf ('%f(%f) ', R_K(i,1,j), R_K(i,2,j));
        end    
        fprintf ('\nnmi(%d): ', j);
        for i=1:numiter
            fprintf ('%f(%f) ', nmi_K(i,1,j), nmi_K(i,2,j));
        end
        fprintf ('\navgent(%d): ', j);
        for i=1:numiter
            fprintf ('%f(%f) ', avgent_K(i,1,j), avgent_K(i,2,j));
        end
        fprintf ('\nAR(%d): ', j);
        for i=1:numiter
            fprintf ('%f(%f) ', AR_K(i,1,j), AR_K(i,2,j));
        end
        fprintf ('\n');
    end
    fprintf ('\n');
    
    fprintf ('Clustering performance of individual views: \n');
    for j=1:num_views
        fprintf ('NMI (view %d): %f (%f)\n', j, nmi(j,1), nmi(j,2));
    end
    fprintf ('\n\nClustering performance with Co-trained spectral algorithm: \n');
    for j=1:num_views
        fprintf ('maximum NMI (view %d): ', j);
        [nmi_max_view(j) ind(j)] = max(nmi_K(:,1,j));
        fprintf ('%f (%f) (obtained in iter %d)\n', nmi_max_view(j), nmi_K(ind(j),2,j), ind(j));        
    end
    [nmi_max r c] = maxMatrix(permute(nmi_K(:,1,:), [1 3 2]));
    fprintf ('maximum NMI: %f (%f) (obtained for view %d, iter %d)\n', nmi_max, nmi_K(r,2,c), c, r);    

    
function  [K] = renormalize(K)
    
    mn = min(min(K));
    mx = max(max(K));
    if (mn < 0)
        K = (K - mn) / (mx-mn);
        K = (K+K')/2;
        K = K - mn;
    end
    