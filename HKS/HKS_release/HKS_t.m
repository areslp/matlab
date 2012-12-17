function [hks_t] = HKS(evecs, evals, A, t, scale)

% INPUTS
%  evecs:  ith each column in this matrix is the ith eigenfunction of the Laplace-Beltrami operator
%  evals:  ith element in this vector is the ith eigenvalue of the Laplace-Beltrami operator
%  A:      ith element in this vector is the area associated with the ith vertex
%  t:      the time scale 
%  scale:  if scale = true, output the scaled hks
%          o.w. ouput the not scaled hks

% OUTPUTS
%  hks_t: heat kernel signature restricted to the time scale t


   %area = sum(A);
   %A = (1/area) * A;
   %evals = area * evals;
   %evecs = sqrt(area) * evecs;
                     
   if scale == true, 
      hks_t = abs( evecs(:, 2:end) ).^2 * exp( ( abs(evals(2)) - abs(evals(2:end)) )  * t);
      Am = sparse([1:length(A)], [1:length(A)], A);
      colsum = sum(Am*hks_t);
      scale = 1.0 / colsum; 
      scalem = sparse([1:length(scale)], [1:length(scale)], scale);
      hks_t = hks_t * scalem;
   else
      hks_t = abs( evecs(:, 2:end) ).^2 * exp( - abs(evals(2:end)) * t);

   end


