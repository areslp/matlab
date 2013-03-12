function [M,IND] = combn(V,N)
% COMBN - all combinations of elements
% M = COMBN(V,N) returns all combinations of N elements of the elements in
% vector V. M has the size (length(V).^N)-by-N.
%
% [M,I] = COMBN(V,N) also returns the index matrix I so that M = V(I).
%
% V can be an array of numbers, cells or strings.
%
% Example:
%   M = COMBN([0 1],3) returns the 8-by-3 matrix:
%     0     0     0
%     0     0     1
%     0     1     0
%     0     1     1
%     ...
%     1     1     1
%
% All elements in V are regarded as unique, so M = COMBN([2 2],3) returns 
% a 8-by-3 matrix with all elements equal to 2.
%
% NB Matrix sizes increases exponentially at rate (n^N)*N.
% 
% COMBN is very fast using a single matrix multiplication, without any
% explicit for-loops. 
%
% See also PERMS, NCHOOSEK
% and ALLCOMB, PERMPOS on the File Exchange

% for Matlab R13
% version 3.2 (may 2006)
% (c) Jos van der Geest
% email: jos@jasen.nl

% History
% 1.1 updated help text
% 2.0 new faster algorithm
% 3.0 (aug 2006) implemented very fast algorithm
% 3.1 (may 2007) Improved algorithm Roger Stafford pointed out that for some values, the floor
% operation on floating points, according to the IEEE 754 standard, could return 
% erroneous values. His excellent solution was to add (1/2) to the values
% of A. 
% 3.2 (may 2007) changed help and error messages slightly

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


if nargin ~=2,
    error('Two arguments required.') ;
end

if isempty(V) || N == 0,
    M = [] ;
    IND = [] ;
elseif fix(N) ~= N || N < 1 || numel(N) ~= 1 ;
    error('Second argument should be a positive integer.') ;
else
    nV = numel(V) ;
    % use a math trick
    A = [0:nV^N-1]+(1/2) ;
    B = [nV.^(1-N:0)] ;
    IND = rem(floor((A(:) * B(:)')),nV) + 1 ;
    M = V(IND) ;     
end

% Previous algorithms

% Version 2.0 
%     for i = N:-1:1
%         X = repmat(1:nV,nV^(N-i),nV^(i-1));
%         IND(:,i) = X(:);
%     end
%     M = V(IND) ;

% Version 1.0
%     nV = numel(V) ;
%     % don waste space, if only one output is requested
%     [IND{1:N}] = ndgrid(1:nV) ;
%     IND = fliplr(reshape(cat(ndims(IND{1}),IND{:}),[],N)) ;
%     M = V(IND) ;

% Combinations using for-loops
% can be implemented in C or VB
% nv = length(V) ;
% C = zeros(nv^N,N) ; % declaration
% for ii=1:N,     
%     cc = 1 ;
%     for jj=1:(nv^(ii-1)),
%         for kk=1:nv,
%             for mm=1:(nv^(N-ii)),
%                 C(cc,ii) = V(kk) ;
%                 cc = cc + 1 ;
%             end
%         end
%     end
% end  


