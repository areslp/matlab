function [ M ] = l21( Q, lambda )
QN=sqrt(sum(Q.^2,1));
% assert(length(find(QN<=0))==0);
coefficient=(QN-lambda)./QN;
invalid=find(coefficient<=0);
coefficient(invalid)=0;
M=bsxfun(@times,Q,coefficient);
