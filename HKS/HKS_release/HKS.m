function [hks] = HKS(evecs, evals, A, scale)

% INPUTS
%  evecs:  ith each column in this matrix is the ith eigenfunction of the Laplace-Beltrami operator
%  evals:  ith element in this vector is the ith eigenvalue of the Laplace-Beltrami operator
%  A:      ith element in this vector is the area associated with the ith vertex
%  scale:  if scale = true, output the scaled hks
%          o.w. ouput the hks that is not scaled

% OUTPUTS
%  hks: ith row in this matrix is the heat kernel signature of the ith vertex


   %area = sum(A);
   %A = (1/area) * A;
   %evals = area * evals;
   %evecs = sqrt(area) * evecs;
                     
   tmin = abs(4*log(10) / evals(end));
   tmax = abs(4*log(10) / evals(2));
   nstep = 100;

   stepsize = (log(tmax) - log(tmin)) / nstep;
   logts = log(tmin):stepsize:log(tmax);
   ts = exp(logts);

   if scale == true, 
      hks = abs( evecs(:, 2:end) ).^2 * exp( ( abs(evals(2)) - abs(evals(2:end)) )  * ts);
      Am = sparse([1:length(A)], [1:length(A)], A);
      colsum = sum(Am*hks);
      scale = 1.0./ colsum; 
      scalem = sparse([1:length(scale)], [1:length(scale)], scale);
      hks = hks * scalem;
   else
      hks = abs( evecs(:, 2:end) ).^2 * exp( - abs(evals(2:end)) * ts);

   end


