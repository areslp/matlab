%% Compare speeds of algorithms on Basis Pursuit problem
%  min ||x||_1 s.t. b = Ax

clear ;
clc ;

addpath L1Solvers;
% You need to download the TFOCS and SESOP_PACK packages from their respective authors
addpath \TFOCS;
addpath(genpath('\SESOP_PACK'));


%% Initialize variables

% l1Method = {'Homotopy'};
l1Method = {'Homotopy', 'PALM', 'DALM', 'PDIPA', 'L1LS', 'FISTA', 'SesopPCD','AMP', 'TFOCS'} ;
methodColor = {'r','g','b','m','c','k','r--','b--','k--'} ;

numTrials = 20 ; % no. of trials for averaging speed

numMethods = length(l1Method) ;

timeTaken = zeros(numTrials,numMethods) ;
errorEst = zeros(numTrials,numMethods) ;

avgTime = zeros(1,numMethods) ;

avgError = zeros(1,numMethods) ;

d = 800 ; % no. of observations
n = 1000 ; % length of unknown signal (x)
k = 100 ; % no. of non-zero entries in x

STOPPING_TIME = -2 ;
STOPPING_GROUND_TRUTH = -1;
STOPPING_DUALITY_GAP = 1;
STOPPING_SPARSE_SUPPORT = 2;
STOPPING_OBJECTIVE_VALUE = 3;
STOPPING_SUBGRADIENT = 4;

maxTime = 8;

%% Run various algorithms
timeEst = cell(numTrials,numMethods);
errEst = cell(numTrials,numMethods);

for trialIndex = 1 : numTrials
    
    disp(['After Trial ' num2str(trialIndex)]) ;
    
    % Generate random data and observations
    A = randn(d,n) ; % observation matrix
    
    for i = 1 : n
        A(:,i) = A(:,i)/norm(A(:,i)) ;
    end
    
    x = zeros(n,1) ;
    p = randperm(n) ;
    x(p(1:k)) = 20*(rand(k,1)-0.5) ; % signal to be estimated
    
    b = A * x ; % observation vector b
    
    for  methodIndex = 1 : numMethods
        
        methodName = ['Solve' l1Method{methodIndex}] ;
        
        methodHandle = str2func(methodName) ;       
        
        tic ;
        [xEst, nIter, curTimeEst, curErrEst] = methodHandle(A, b, 'stoppingCriterion', STOPPING_TIME,...
            'groundTruth', x, 'maxtime', maxTime, 'maxiteration', 1e6) ;
        tEst = toc ;       
        
        timeEst{trialIndex,methodIndex} = curTimeEst;
        errEst{trialIndex,methodIndex} = curErrEst /norm(x);
        
        estError = norm(x-xEst)/norm(x) ;
        timeTaken(trialIndex,methodIndex) = tEst ;
        errorEst(trialIndex,methodIndex) = estError ;
        
        avgTime(methodIndex) = avgTime(methodIndex) + tEst ;
        avgError(methodIndex) = avgError(methodIndex) + estError ;
        
        fprintf(1,'%10s\t\t%e\t\t%e\n',l1Method{methodIndex},avgTime(methodIndex)/trialIndex,...
            avgError(methodIndex)/trialIndex) ;
        pause(0.2) ;
        
    end
    
    disp('----------------') ;
    
    avgCurve = cell(1, numMethods);
    avgCount = cell(1, numMethods);
    
    figure(1) ; clf ;
    for methodIndex = 1 : numMethods
        avgCurve{methodIndex} = zeros(1,maxTime * 100);
        avgCount{methodIndex} = zeros(1,maxTime * 100);
        for tIdx = 1 : trialIndex
            curTimeEst = timeEst{tIdx,methodIndex};
            curErrEst = errEst{tIdx,methodIndex};
            
            idx_length = min(maxTime * 100, (100*max(curTimeEst)+1));
            cur_pt = 1;
            for idx0 = 1:idx_length
                while cur_pt <= length(curErrEst) && curTimeEst(cur_pt) < .01 * idx0
                    cur_pt = cur_pt + 1;
                end
                if cur_pt > 1
                    avgCurve{methodIndex}(idx0) = avgCurve{methodIndex}(idx0) + curErrEst(cur_pt-1);
                else
                    avgCurve{methodIndex}(idx0) = avgCurve{methodIndex}(idx0) + 1;
                end
                avgCount{methodIndex}(idx0) = avgCount{methodIndex}(idx0) + 1;
            end
            
        end
        maxCount = max(find(avgCount{methodIndex} ~= 0));
        avgCount{methodIndex} = avgCount{methodIndex}(1:maxCount);
        avgCurve{methodIndex} = avgCurve{methodIndex}(1:maxCount) ./ avgCount{methodIndex};
        
        semilogy(0.01*(1:maxCount),avgCurve{methodIndex},methodColor{methodIndex},'LineWidth',3) ;
        hold on;
    end
    xlim([0,maxTime]);
    ylim([1e-20 1]);
    set(gca,'fontsize',16);
    legend(l1Method,'fontsize',12);
    xlabel('Time(s)','fontsize',20);
    ylabel('Relative Error of x','fontsize',20);
    hold off;
    pause(.1);
end

avgTime = avgTime/numTrials ;
avgError = avgError/numTrials ;
