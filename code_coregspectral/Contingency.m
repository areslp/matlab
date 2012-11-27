function Cont=Contingency(Mem1,Mem2)
%CONTINGENCY Form contigency matrix for two vectors
% C=Contingency(Mem1,Mem2) returns contingency matrix for two
% column vectors Mem1, Mem2. These define which cluster each entity 
% has been assigned to.
%
% See also RANDINDEX.
%

%(C) David Corney (2000)   		D.Corney@cs.ucl.ac.uk

if nargin < 2 | min(size(Mem1)) > 1 | min(size(Mem2)) > 1
   error('Contingency: Requires two vector arguments')
   return
end

Cont=zeros(max(Mem1),max(Mem2));

for i = 1:length(Mem1);
   Cont(Mem1(i),Mem2(i))=Cont(Mem1(i),Mem2(i))+1;
end
