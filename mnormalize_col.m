function [M]=mnormalize_col(M)
% matrix normalize column
MN=sqrt(sum(M.^2,1));
coefficient=MN;
invalid=find(coefficient<=0);
coefficient=coefficient.^-1;
coefficient(invalid)=0;
M=bsxfun(@times,M,coefficient);

