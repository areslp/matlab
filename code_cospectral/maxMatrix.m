function [val r c] = maxMatrix(A)

   minval = min(min(abs(A)));
   A = A + minval*1e-10*abs(rand(size(A,1),size(A,2)));
   rowind = 1:size(A,1);
   colind = 1:size(A,2);
   val = max(max(A));
   indicator = (A==val);
   r = rowind(sum(indicator,2)==1);
   c = colind(sum(indicator,1)==1);
