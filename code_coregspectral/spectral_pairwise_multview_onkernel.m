function [F P R nmi avgent AR] = spectral_pairwise_multview_onkernel(S,num_views,numClust,lambda,truth,numiter)
% INPUT:
% OUTPUT:
    
    if (min(truth)==0)
        truth = truth + 1;
    end
    
    [N M1] = size(S{1});
    %[N M2] = size(X2);
    
        
    kmeans_avg_iter = 20;
    opts.disp = 0;

    numEV = numClust;
    numVects = numClust;
    for i=1:num_views
    % Laplacian for the first view of the data
        fprintf('computing kernel for X(%d)\n',i);
        K(:,:,i) = S{i}; %constructKernel(X{i},X{i},options(i));
        %K1 = X1*X1';
        D = diag(sum(K(:,:,i),1));
        %L1 = D1 - K1; 
        L(:,:,i) = sqrt(inv(D))*K(:,:,i)*sqrt(inv(D));  
        L(:,:,i)=(L(:,:,i)+L(:,:,i)')/2;
        [U(:,:,i) E] = eigs(L(:,:,i),numEV,'LA',opts);    
        objval(i,1) = sum(diag(E));
    end
    
    %%do clustering for first view
    U1 = U(:,:,1);
    normvect = sqrt(diag(U1*U1'));
    normvect(find(normvect==0.0)) = 1;
    U1 = inv(diag(normvect)) * U1;    
    for j=1:kmeans_avg_iter
        C = kmeans(U1(:,1:numVects),numClust,'EmptyAction','drop'); 
        [Fj(j),Pj(j),Rj(j)] = compute_f(truth,C); 
        [Aj nmi_j(j) avgent_j(j)] = compute_nmi(truth,C);
        [ARj(j),RIj(j),MIj(j),HIj(j)]=RandIndex(truth,C);
    end
    F(1) = mean(Fj); std_F(1) = std(Fj);
    P(1) = mean(Pj); std_P(1) = std(Pj);
    R(1) = mean(Rj); std_R(1) = std(Rj);
    nmi(1) = mean(nmi_j); std_nmi(1) = std(nmi_j);
    avgent(1) = mean(avgent_j); std_avgent(1) = std(avgent_j);
    AR(1) = mean(ARj); std_AR(1) = std(ARj);
    

    i = 2;
    % now iteratively solve for all U's
    while(i<=numiter+1)
        fprintf('Running iteration %d\n',i-1);
        
        for j=2:num_views            
            L_complement(1:N,1:N) = 0;
            for k=1:num_views
                if (k==j) 
                    continue;
                end
                L_complement =   L_complement + U(:,:,k)*U(:,:,k)';
            end
            L_complement = (L_complement+L_complement')/2;           
            [U(:,:,j) E] = eigs(L(:,:,j) + lambda*L_complement, numEV,'LA',opts);    
            objval(j,i) = sum(diag(E));
        end

        L_complement(1:N,1:N) = 0;
        for k=1:num_views
            if (k==1) 
                continue;
            end
            L_complement =   L_complement + U(:,:,k)*U(:,:,k)';
        end
        L_complement = (L_complement+L_complement')/2;
        [U(:,:,1) E] = eigs(L(:,:,1) + lambda*L_complement, numEV,'LA',opts);    
        objval(1,i) = sum(diag(E));
                
        if (1)  %use view 1 in actual clustering
            U1 = U(:,:,1);
            normvect = sqrt(diag(U1*U1'));    
            normvect(find(normvect==0.0)) = 1;
            U1 = inv(diag(normvect)) * U1;
            
            for j=1:kmeans_avg_iter
                C = kmeans(U1(:,1:numVects),numClust,'EmptyAction','drop'); 
                [Fj(j),Pj(j),Rj(j)] = compute_f(truth,C); 
                [Aj nmi_j(j) avgent_j(j)] = compute_nmi(truth,C);
                [ARj(j),RIj(j),MIj(j),HIj(j)]=RandIndex(truth,C);
            end
            F(i) = mean(Fj); std_F(i) = std(Fj);
            P(i) = mean(Pj); std_P(i) = std(Pj);
            R(i) = mean(Rj); std_R(i) = std(Rj); 
            nmi(i) = mean(nmi_j); std_nmi(i) = std(nmi_j);
            avgent(i) = mean(avgent_j); std_avgent(i) = std(avgent_j);
            AR(i) = mean(ARj); std_AR(i) = std(ARj);
        end
        i = i+1;
    end
    
    
    for j=1:num_views
        normvect = sqrt(diag(U(:,:,j)*U(:,:,j)'));
        normvect(find(normvect==0.0)) = 1;
        U_norm(:,:,j) = inv(diag(normvect)) * U(:,:,j);
    end

    V = U(:,:,1);
    for j=2:num_views
        V = [V U(:,:,j)];
    end
    normvect = sqrt(diag(V*V'));
    normvect(find(normvect==0.0)) = 1;
    V = inv(diag(normvect)) * V;
    %U = U./repmat(sqrt(sum(U.*U,2)),1,numClust*2); % normalize
    for j=1:kmeans_avg_iter
        C = kmeans(V(:,1:numVects),numClust,'EmptyAction','drop'); 
        [Fj(j),Pj(j),Rj(j)] = compute_f(truth,C); 
        [Aj nmi_j(j) avgent_j(j)] = compute_nmi(truth,C);
        [ARj(j),RIj(j),MIj(j),HIj(j)]=RandIndex(truth,C);
    end
            F(i) = mean(Fj); std_F(i) = std(Fj);
            P(i) = mean(Pj); std_P(i) = std(Pj);
            R(i) = mean(Rj); std_R(i) = std(Rj); 
            nmi(i) = mean(nmi_j); std_nmi(i) = std(nmi_j);
            avgent(i) = mean(avgent_j); std_avgent(i) = std(avgent_j);
            AR(i) = mean(ARj); std_AR(i) = std(ARj);
    
    %%%CCA on U1 and U2
    %i = i+1;
    %[feats1 feats2 F_c P_c R_c nmi_c avgent_c] = multiviewccacluster(U1_norm, U2_norm, numClust, sigma1, sigma2, truth);
    
    fprintf('F:   ');
    for i=1:numiter+2
        fprintf('%f(%f)  ', F(i), std_F(i));
    end
    fprintf('\n\n');
    fprintf('P:   ');    
    for i=1:numiter+2      
        fprintf('%f(%f)  ', P(i), std_P(i));
    end
    fprintf('\n\n');
    fprintf('R:   ');    
    for i=1:numiter+2      
        fprintf('%f(%f)  ', R(i), std_R(i));
    end
    fprintf('\n\n');
    fprintf('nmi:   ');    
    for i=1:numiter+2      
        fprintf('%f(%f)  ', nmi(i), std_nmi(i));
    end
    fprintf('\n\n');
    fprintf('avgent:   ');    
    for i=1:numiter+2      
        fprintf('%f(%f)  ', avgent(i), std_avgent(i));
    end
    fprintf('\n\n');
    fprintf('AR:   ');    
    for i=1:numiter+2      
        fprintf('%f(%f)  ', AR(i), std_AR(i));
    end
    fprintf('\n\n');
    for j=1:num_views
        fprintf('objval_u%d:   ', j);    
        for i=1:numiter+1
            fprintf('%f  ', objval(j,i));
        end
        fprintf('\n');
    end
        
    if (0)
    %%%%averaging of U1 and U2
    V = (U1_norm+U2_norm)/2;
    normvect = sqrt(diag(V*V'));
    normvect(find(normvect==0.0)) = 1;
    V = inv(diag(normvect)) * V;
    %U = U./repmat(sqrt(sum(U.*U,2)),1,numClust*2); % normalize
    for j=1:kmeans_avg_iter
        C = kmeans(V(:,1:numVects),numClust,'EmptyAction','drop'); 
        [Fj(j),Pj(j),Rj(j)] = compute_f(truth,C); 
        [Aj nmi_j(j) avgent_j(j)] = compute_nmi(truth,C);
        [ARj(j),RIj(j),MIj(j),HIj(j)]=RandIndex(truth+1,C);
    end
    i = i+1;
    F(i) = mean(Fj);
    P(i) = mean(Pj);
    R(i) = mean(Rj);
    nmi(i) = mean(nmi_j);
    avgent(i) = mean(avgent_j);
    AR(i) = mean(ARj);    
    
    %C = kmeans(U,numClust,'EmptyAction','drop');  
    %[F(i),P(i),R(i)] = compute_f(truth,C); 
    %[A nmi(i) avgent(i)] = compute_nmi(truth,C);
    
    end
    
function  [K] = renormalize(K)
    
    mn = min(min(K));
    mx = max(max(K));
    if (mn < 0)
        %K = (K - mn) / (mx-mn);
        %K = (K+K')/2;
        K = K - mn;
    end    