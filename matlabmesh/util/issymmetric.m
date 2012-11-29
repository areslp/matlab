function [ issym, err ] = issymmetric( M, use_epsilon )
% [ issym, err ] = issymmetric( M, use_epsilon )
%   returns 1 if M is symmetric (within tolerance use_epsilon)

if ~exist('use_epsilon','var')
    use_epsilon = sqrt(eps);
end

err = sum(sum(abs(M-M')));
if err > use_epsilon
    issym = 0;
else
    issym = 1;
end
    


end
