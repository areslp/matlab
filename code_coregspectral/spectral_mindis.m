function [F P R nmi avgent AR] = spectral_mindis(X1, X2, numClust,sigma1,sigma2, truth)

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
        
    
    K1 = constructKernel(X1,X1,options1);
    K2 = constructKernel(X2,X2,options2);
    
    W = K1*K2;
    D_row = sum(W,2); %sum each row
    D_col = sum(W,1); %sum each col
    Lw = diag(D_row.^(-1/2))*W*diag(D_col.^(-1/2));
    [U S V] = svds(Lw,numClust);
    X = [U; V];
    normvect = sqrt(diag(X*X'));
    normvect(find(normvect==0.0)) = 1;
    Y = inv(diag(normvect)) * X;
    
    for i=1:20
        %view1
        C1 = kmeans(Y(1:N,:),numClust,'EmptyAction','drop');
        [F1i(i),P1i(i),R1i(i)] = compute_f(truth+1,C1);
        [A1 nmi1i(i) avgent1i(i)] = compute_nmi(truth,C1);
        [AR1i(i),RI1i(i),MI1i(i),HI1i(i)]=RandIndex(truth,C1);       
        %view2
        C2 = kmeans(Y(N+1:2*N,:),numClust,'EmptyAction','drop');
        [F2i(i),P2i(i),R2i(i)] = compute_f(truth+1,C2);
        [A2 nmi2i(i) avgent2i(i)] = compute_nmi(truth,C2);
        [AR2i(i),RI2i(i),MI2i(i),HI2i(i)]=RandIndex(truth,C2);        
        %joint
        C = kmeans((Y(1:N,:)+Y(N+1:2*N,:))/2, numClust, 'EmptyAction','drop');
        [Fi(i),Pi(i),Ri(i)] = compute_f(truth+1,C);
        [A nmii(i) avgenti(i)] = compute_nmi(truth,C);
        [ARi(i),RIi(i),MIi(i),HIi(i)]=RandIndex(truth,C);
    end
    F(1) = mean(F1i);
    P(1) = mean(P1i);
    R(1) = mean(R1i);
    nmi(1) = mean(nmi1i);
    avgent(1) = mean(avgent1i);
    AR(1) = mean(AR1i);
    F(2) = mean(F2i);
    P(2) = mean(P2i);
    R(2) = mean(R2i);
    nmi(2) = mean(nmi2i);
    avgent(2) = mean(avgent2i);
    AR(2) = mean(AR2i);
    F(3) = mean(Fi);
    P(3) = mean(Pi);
    R(3) = mean(Ri);
    nmi(3) = mean(nmii);
    avgent(3) = mean(avgenti);
    AR(3) = mean(ARi);    
    
    fprintf('nmi_1 = %f(%f), nmi_2 = %f(%f), nmi = %f(%f)\n', nmi(1), std(nmi1i),nmi(2), std(nmi2i),nmi(3), std(nmii));
    fprintf('F_1 = %f(%f), F_2 = %f(%f), F = %f(%f)\n', F(1), std(F1i),F(2), std(F2i),F(3), std(Fi));
    fprintf('P_1 = %f(%f), P_2 = %f(%f), P = %f(%f)\n', P(1), std(P1i),P(2), std(P2i),P(3), std(Pi));
    fprintf('R_1 = %f(%f), R_2 = %f(%f), R = %f(%f)\n', R(1), std(R1i),R(2), std(R2i),R(3), std(Ri));
    fprintf('Entropy_1 = %f(%f), Entropy_2 = %f(%f), Entropy = %f(%f)\n', avgent(1), std(avgent1i),avgent(2), std(avgent2i),avgent(3), std(avgenti));
    fprintf('AR_1 = %f(%f), AR_2 = %f(%f), AR = %f(%f)\n', AR(1), std(AR1i),AR(2), std(AR2i),AR(3), std(ARi));
    
    
    
function  [K] = renormalize(K)
    
    mn = min(min(K));
    mx = max(max(K));
    if (mn < 0)
        %K = (K - mn) / (mx-mn);
        %K = (K+K')/2;
        K = K - mn;
    end    
    
    