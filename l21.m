function [ M ] = l21( Q, lambda )
QN=sqrt(sum(Q.^2,1));
zero=find(QN==0);
coefficient=(QN-lambda)./QN; % may be NaN
coefficient(zero)=0;
invalid=find(coefficient<=0);
coefficient(invalid)=0;
M=bsxfun(@times,Q,coefficient);
