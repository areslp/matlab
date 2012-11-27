function [U U1 F P R nmi avgent AR] = spectral_pairwise(X1,X2,numClust,sigma1,sigma2,lambda,truth,numiter)
% INPUT:
% OUTPUT:
    
    if (min(truth)==0)
        truth = truth + 1;
    end
    
    [N M1] = size(X1);
    [N M2] = size(X2);
    options1 = [];
    options1.KernelType = 'Gaussian';
    options1.t = sigma1;
    options1.d = 4;
    
    options2 = [];
    options2.KernelType = 'Gaussian';
    options2.t = sigma2;
    options2.d = 4;
    
    kmeans_avg_iter = 20;
    opts.disp = 0;
    
    % Laplacian for the first view of the data
    fprintf('computing kernel for X1\n');
    K1 = constructKernel(X1,X1,options1);
    %K1 = X1*X1';
    D1 = diag(sum(K1,1));
    %L1 = D1 - K1; 
    L1 = sqrt(inv(D1))*K1*sqrt(inv(D1));  
    L1=(L1+L1')/2;
    
    % Laplacian for the second view of the data
    fprintf('computing kernel for X2\n');    
    K2 = constructKernel(X2,X2,options2);
    %K2 = X2*X2';
    D2 = diag(sum(K2,1));
    %L2 = D2 - K2; 
    L2 = sqrt(inv(D2))*K2*sqrt(inv(D2));    
    L2 = (L2+L2')/2;
    
    numEV = numClust;
    numVects = numClust;
    % initialize U1 (using only the 1st view)
    %[U1 E1] = eigs(L1, numClust, 'SM');
    [U1 E1] = eigs(L1,numEV,'LA',opts);   
    objval_u1(1) = sum(diag(E1));
    %[U1 E1] = eig(L1);
    %[E11 I] = sort(diag(E1));  %sort in increasing order    
    %U1 = U1(:,I(end-numClust+1:end));
    U = U1;
    %U = U./repmat(sqrt(sum(U.*U,2)),1,numClust); % normalize
    normvect = sqrt(diag(U*U'));
    normvect(find(normvect==0.0)) = 1;
    U = inv(diag(normvect)) * U;
    
    for j=1:kmeans_avg_iter
        C = kmeans(U(:,1:numVects),numClust,'EmptyAction','drop'); 
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
    
    [U2 E2] = eigs(L2,numEV,'LA',opts);    
    objval_u2(1) = sum(diag(E2));


    i = 2;
    % now iteratively solve for U1 and U2
    while(i<=numiter+1)
        fprintf('Running iteration %d\n',i-1);
        % update U2 (keeping U1 fixed)
        %normvect = sqrt(diag(U1*U1'));    
        %normvect(find(normvect==0.0)) = 1;
        %U1 = inv(diag(normvect)) * U1;    %use normalized U1        
        %K_u1 = U1*U1'; D_u1 = diag(sum(K_u1,1)); L_u1 = sqrt(inv(D_u1))*K_u1*sqrt(inv(D_u1));        
        L_u1 =   U1*U1';
        %K_u1 = constructKernel(U1,U1,options1); D_u1 = diag(sum(K_u1,1)); L_u1 = sqrt(inv(D_u1))*K_u1*sqrt(inv(D_u1));
        L_u1 = (L_u1+L_u1')/2;
        [U2 E2] = eigs(L2 + lambda*L_u1, numEV,'LA',opts);    
        objval_u2(i) = sum(diag(E2));
        %[U2 E2] = eig(L2 + lambda*L_u1);
        %[E22 I] = sort(diag(E2));  %sort in increasing order        
        %U2 = U2(:,I(end-numClust+1:end));
        
        if (0)  %use view 2 in actual clustering
            U = U2;
            normvect = sqrt(diag(U*U'));    
            normvect(find(normvect==0.0)) = 1;
            U = inv(diag(normvect)) * U;
            %U = U./repmat(sqrt(sum(U.*U,2)),1,numClust); % normalize
            C = kmeans(U,numClust,'EmptyAction','drop');  
            [F(i),P(i),R(i)] = compute_f(truth,C); 
            [A nmi(i) avgent(1)] = compute_nmi(truth,C);
        end
        %U2 = U2./repmat(sqrt(sum(U2.*U2,2)),1,numClust); % normalize
        % update U1 (keeping U2 fixed)
        %U2 = U2./repmat(sqrt(sum(U2.*U2,2)),1,numClust);   %use normalized U2        
        %K_u2 = U2*U2'; D_u2 = diag(sum(K_u2,1)); L_u2 = sqrt(inv(D_u2))*K_u2*sqrt(inv(D_u2));        
        L_u2 =   U2*U2';
        %K_u2 = constructKernel(U2,U2,options2); D_u2 = diag(sum(K_u2,1)); L_u2 = sqrt(inv(D_u2))*K_u2*sqrt(inv(D_u2));
        L_u2 = (L_u2+L_u2')/2;
        [U1 E1] = eigs(L1 + lambda*L_u2, numEV,'LA',opts);    
        objval_u1(i) = sum(diag(E1));
        %[U1 E1] = eig(L1 + lambda*L_u2);
        %[E11 I] = sort(diag(E1));  %sort in increasing order
        %U1 = U1(:,I(end-numClust+1:end));
        
        if (1)  %use view 1 in actual clustering
            U = U1;
            normvect = sqrt(diag(U*U'));    
            normvect(find(normvect==0.0)) = 1;
            U = inv(diag(normvect)) * U;
            
            for j=1:kmeans_avg_iter
                C = kmeans(U(:,1:numVects),numClust,'EmptyAction','drop'); 
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
    
%     U1 = U1./repmat(sqrt(sum(U1.*U1,2)),1,numClust); % normalize
%     C = kmeans(U1,numClust,'EmptyAction','drop'); 
%     %[F,P,R] = compute_f(truth,C);    
    
    normvect = sqrt(diag(U1*U1'));
    normvect(find(normvect==0.0)) = 1;
    U1_norm = inv(diag(normvect)) * U1;

    normvect = sqrt(diag(U2*U2'));
    normvect(find(normvect==0.0)) = 1;
    U2_norm = inv(diag(normvect)) * U2;
    
    U = [U1 U2];
    normvect = sqrt(diag(U*U'));
    normvect(find(normvect==0.0)) = 1;
    U = inv(diag(normvect)) * U;
    %U = U./repmat(sqrt(sum(U.*U,2)),1,numClust*2); % normalize
    for j=1:kmeans_avg_iter
        C = kmeans(U(:,1:numVects),numClust,'EmptyAction','drop'); 
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
    fprintf('objval_u1:   ');    
    for i=1:numiter+1
        fprintf('%f  ', objval_u1(i));
    end
    fprintf('\n');
    fprintf('objval_u2:   ');    
    for i=1:numiter+1
        fprintf('%f  ', objval_u2(i));
    end
    fprintf('\n');
        
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