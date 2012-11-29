function [ cotAB ] = vcot( A, B )
% cotAB = vcot( A, B)
%  return cotangent of angle between n-D vectors A and B
%  Formula from page 19 of [Meyer02]
%  [TODO] make this work for arrays of vectors

tmp = dot(A,B);
cotAB = tmp / sqrt( dot(A,A)*dot(B,B) - tmp*tmp );

end
